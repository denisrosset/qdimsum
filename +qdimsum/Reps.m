classdef Reps
% Numerical algorithms to decompose a representation given by generalized permutations into 
% isotypic/irreducible components

    methods (Static) % Helpers
        
        function gdfo = groupDecompositionForOrbit(group, orbit)
            groupDecomposition = group.decomposition;
            n = size(groupDecomposition{1}, 2);
            backIndex = zeros(n, 1);
            backIndex(orbit) = 1:length(orbit);
            computeBackIndex = @(el) sign(el).*backIndex(abs(el));
            gdfo = cellfun(@(groupElements) computeBackIndex(groupElements(:, orbit)), groupDecomposition, 'UniformOutput', false);
        end
        
        function [U reps] = irreducibleDecomposition(group, settings)
        % Given a generalized permutation group provided by the given decomposition, 
        % returns a change of basis matrix U such that
        %
        % U' * matrix * U is block diagonal
        %
        % with the structure given by the representation dimensions reps(1,:) and multiplicities reps(2,:)
            import qdimsum.*
            groupDecomposition = group.decomposition;
            n = size(groupDecomposition{1}, 2);            
            sample1 = Reps.sampleSymmetricMatrix(group);
            sample2 = Reps.sampleSymmetricMatrix(group);
            sample3 = Reps.sampleGenericMatrix(group);
            if settings.checkLevel > 0
                sample4 = Reps.sampleSymmetricMatrix(group);
            else
                sample4 = [];
            end
            orbits = Reps.findOrbits(group, settings);
            nOrbits = length(orbits);
            gfo = cell(1, nOrbits);
            for i = 1:nOrbits
                gfo{i} = group.restrictToOrbit(orbits{i});
            end
            [U reps colOrbit] = Reps.isotypicComponentsFromSamples(orbits, sample1, sample2, settings, sample4);
            shift = 0;
            nReps = size(reps, 2);
            for r = 1:nReps
                d = reps(1, r); % dimension
                m = reps(2, r); % multiplicity
                repRange = shift+(1:d*m); % range of indices corresponding to the current isotypic component
                colOrbitRep = colOrbit(repRange); % orbit number for each column considered
                if settings.blockDiagRefine
                    Urep = zeros(n, d*m);
                    for i = 1:nOrbits
                        colRangeRep = find(colOrbitRep == i); % which columns for the current orbit
                        if ~isempty(colRangeRep)
                            rowRange = orbits{i}; % which rows for the current orbit
                            UrepOrbit = U(rowRange, repRange(colRangeRep));
                            UrepOrbit = Reps.refineIsotypicSubspace(gfo{i}, UrepOrbit, settings);
                            Urep(rowRange, colRangeRep) = UrepOrbit;
                        end
                    end
                else
                    Urep = U(:, repRange);
                end
                if m > 1
                    genSample = Urep'*sample3*Urep;
                    switch Reps.identifyRepresentationType(genSample, settings)
                      case 'R'
                        t1 = Urep'*sample1*Urep;
                        t2 = Urep'*sample2*Urep;
                        t1 = t1 + t1';
                        t2 = t2 + t2';
                        V = Reps.decomposeIsotypicRealType(t1, t2, d, m, colOrbitRep);
                        Urep = Urep * V;
                      case 'C'
                        error('Complex representations are not supported');
                      case 'H'
                        error('Quaternionic representations are not supported');
                    end
                end
                U(:, repRange) = Urep;
                shift = shift + d*m;
            end
            if settings.checkLevel > 0
                sample4p = U'*sample4*U;
                sws = Reps.swapSpaces(reps(1,:), reps(2,:));
                sample4p = sample4p(sws, sws);
                sample4pav = Reps.averageBlocksKeepCopies(sample4p, reps(1,:), reps(2,:));
                diff = (sample4p - sample4pav);
                diff = (diff + diff')/2;
                diffeig = abs(eig(diff));
                if ~isequal(settings.blockDiagMatHist, [])
                    settings.blockDiagMatHist.register(diffeig);
                end
                assert(all(diffeig < settings.blockDiagMatTol), 'Irreducible decomposition failed');
            end
        end
        
        function perm = swapSpaces(dimensions, multiplicities)
        % Returns the permutation that switches between the kron(A, id) and kron(id, B) forms
            shift = 0;
            nReps = length(dimensions);
            perm = [];
            for r = 1:nReps
                d = dimensions(r);
                m = multiplicities(r);
                s = reshape(reshape(1:d*m, [d m])', [1 d*m]);
                perm = [perm, shift + s];
                shift = shift + d*m;
            end
        end
        
        function red = averageBlocksKeepCopies(mat, dimensions, multiplicities)
            dRed = sum(multiplicities);
            nReps = length(dimensions);
            blocks = {};
            shift = 0;
            for i = 1:nReps
                d = dimensions(i);
                m = multiplicities(i);
                block = zeros(m, m);
                for j = 1:d
                    range = shift + (1:m);
                    block = block + mat(range, range);
                    shift = shift + m;
                end
                block = block / d;
                for j = 1:d
                    blocks = horzcat(blocks, {block});
                end
            end
            red = blkdiag(blocks{:});
        end

        function sample = sampleGenericMatrix(group)
        % Samples a real matrix that commutes with the given group action.
            import qdimsum.*
            groupDecomposition = group.decomposition;
            n = size(groupDecomposition{1}, 2);
            sample = GenPerm.symmetrize(randn(n), groupDecomposition);
        end
        
        function sample = sampleSymmetricMatrix(group)
        % Samples a real symmetric matrix (in the sense M = M') that additionally commutes with
        % the given group action.
            import qdimsum.*
            groupDecomposition = group.decomposition;
            n = size(groupDecomposition{1}, 2);
            sample = GenPerm.symmetrize(Random.symmetricGaussian(n), groupDecomposition);
            sample = sample + sample';
        end
        
        function [U reps fromOrbit] = isotypicComponents(group, settings)
        % Finds the isotypic components of the (generalized) permutation representation afforded
        % by the given group decomposition.
        %
        % Returns U, a change of basis matrix such that U'*sample*U is block diagonal
        % and reps = double(2,nBlocks) with reps(1,i) the dimension and reps(2,i) the multiplicity.
        %
        % In each isotypic component, the copies of the representations are ordered.
        %
        % Note: does not try to split the group action into orbits, as irreducibleDecomposition does.
            import qdimsum.*
            orbits = Reps.findOrbits(group, settings);
            sample1 = Reps.sampleSymmetricMatrix(group);
            sample2 = Reps.sampleSymmetricMatrix(group);
            if settings.checkLevel > 0
                sample3 = Reps.sampleSymmetricMatrix(group);
                [U reps fromOrbit] = Reps.isotypicComponentsFromSamples(orbits, sample1, sample2, settings, sample3);
            else
                [U reps fromOrbit] = Reps.isotypicComponentsFromSamples(orbits, sample1, sample2, settings);
            end
        end
        
        function [U reps fromOrbit] = isotypicComponentsFromSamples(orbits, sample1, sample2, settings, sample3)
        % Finds the isotypic components of a representation from two generic samples of
        % real, symmetric (M=M') matrices commuting with the group action.
        %
        % See Reps.isotypicComponents for the output format.
        %
        % If settings.checkLevel > 0, then sample3 must be provided.
            import qdimsum.*
            nOrbits = length(orbits);
            n = size(sample1, 1);
            vD = [];
            runs = {};
            U = zeros(n, n);
            fromOrbit = zeros(1, n);
            shift = 0;
            for o = 1:nOrbits
                range = orbits{o};
                orbitSize = length(range);
                basisRange = shift+(1:orbitSize);
                [Uo Do] = eig(sample1(range, range)); % eigenvalue decomposition of one orbit of the first sample
                Do = diag(Do);
                [~, I] = sort(Do); % sort eigenvalues just to be sure
                Do = Do(I); % reorder accordingly
                Uo = Uo(:, I); 
                U(range, basisRange) = Uo;
                fromOrbit(basisRange) = o;
                runso = Reps.findRuns(Do, settings); % find subspaces corresponding to repeated eigenvalues
                runs = horzcat(runs, cellfun(@(r) r + shift, runso, 'UniformOutput', false));
                shift = shift + orbitSize;
            end
            sample2p = U'*sample2*U; % compute the block structure using the second sample
            nRuns = length(runs);
            blockMask = zeros(nRuns, nRuns);
            for i = 1:nRuns
                for j = 1:nRuns
                    blockMask(i, j) = norm(sample2p(runs{i}, runs{j}));
                end
            end
            if ~isequal(settings.blockDiagMatHist, [])
                for i = 1:nRuns
                    for j = 1:nRuns
                        svds = svd(sample2p(runs{i}, runs{j}));
                        settings.blockDiagMatHist.register(svds(svds <= settings.blockDiagMatTol));
                    end
                end
            end
            blockMask = blockMask > settings.blockDiagMatTol;
            % find the subspaces corresponding to the same irreducible representation
            % by looking at the connected components of the graph defined by the adjacency
            % matrix of the block mask
            cc = Reps.findConnectedComponents(blockMask);
            Nc = length(cc);
            reps = zeros(2, Nc);
            for i = 1:Nc
                c = cc{i};
                reps(2, i) = length(c);
                dims = arrayfun(@(i) length(runs{i}), c);
                % verify that the dimensions are all the same for consistency
                assert(all(dims - dims(1) == 0));
                reps(1, i) = dims(1);
            end
            % sort the irreducible representations first by increasing dimension
            % then by increasing multiplicity
            [~, I] = sortrows(reps');
            reps = reps(:, I);
            cc = cc(I);
            reorder = cellfun(@(i) horzcat(runs{i}), cc, 'UniformOutput', false);
            reorder = [reorder{:}];
            U = U(:, reorder);
            fromOrbit = fromOrbit(reorder);
            if settings.checkLevel > 0
                % Checks the result on a third sample
                mask = abs(U'*sample3*U) > settings.blockDiagMatTol;
                shift = 0;
                for i = 1:size(reps, 2)
                    repSize = reps(1, i) * reps(2, i);
                    blockIndices = shift + (1:repSize);
                    restIndices = (shift+repSize+1):size(mask, 1);
                    assert(~any(any(mask(blockIndices, restIndices))), 'Invalid isotypic decomposition');
                    assert(~any(any(mask(restIndices, blockIndices))), 'Invalid isotypic decomposition');
                    shift = shift + repSize;
                end
            end
        end

        function refinedBasis = refineIsotypicSubspace(group, basis, settings)
        % Refines the basis of an isotypic subspace; the basis is given by the column vectors
        % of "basis"
            import qdimsum.*
            n = size(basis, 2);
            if size(basis, 1) == n
                refinedBasis = eye(n);
            else
                T = GenPerm.symmetrize(basis*Random.symmetricGaussian(n)*basis', group.decomposition);
                T = T + T';
                [refinedBasis, lambda] = eig(T);
                lambda = diag(lambda);
                [~, I] = sort(abs(lambda)); % sort eigenvalues (needed on Matlab 2015b)
                I = flipud(I(:)); % largest magnitude first
                lambda = lambda(I); % reorder accordingly
                refinedBasis = refinedBasis(:, I);
                refinedBasis = refinedBasis(:, 1:n); % cuts the additional eigenvector
                if settings.checkLevel > 0
                    assert(abs(lambda(n + 1)) < settings.blockDiagEigTol, ...
                           'The provided approximate basis does not match an isotypic component');                
                end
            end
        end

        function type = identifyRepresentationType(sample, settings)
        % Given a sample of a symmetrized block, returns its representation type
        % genSample = basis'*GenPerm.symmetrize(randn(N, N))*basis
        %
        % Note: the sample must not be a symmetric matrix
        %
        % TODO: support complex and quaternionic representations
            import qdimsum.*
            sampleSym = sample + sample';
            lambda = eig(sample);
            lambdaSym = eig(sampleSym);
            check = real(lambda) + rand*imag(lambda); % TODO: less hacky comparison
            checkSym = real(lambdaSym) + rand*imag(lambdaSym);
            r = rand;
            [~, I] = sort(check);
            [~, Isym] = sort(checkSym);
            runs = Reps.findRuns(check(I), settings);
            runsSym = Reps.findRuns(checkSym(Isym), settings);
            if length(runs) == length(runsSym)
                type = 'R';
            else
                error('No support for the complex or quaternionic case');
            end
            
        end
        
        function V = decomposeIsotypicRealType(sample1, sample2, dim, mult, colOrbit)
        % Given two (symmetric matrix) samples from an isotypic component of real type for a representation
        % of the given dimension and multiplicity, computes the change of basis matrix that puts all copies
        % of the irreducible representation in the same basis.
        %
        % Note that the representation copies should already be ordered in the isotypic basis.
            n = dim * mult;
            if nargin < 5
                colOrbit = ones(1, n);
            end
            Q = zeros(n, n);
            for i = unique(colOrbit)
                cols = find(colOrbit == i);
                [Qo, ~] = eig(sample1(cols, cols));
                Q(cols, cols) = Qo;
            end
            T = Q(:,1:dim)'*sample2*Q;
            P = cell(1, mult);
            P{1} = eye(dim);
            for j = 2:mult
                P{j} = T(:, (j-1)*dim+(1:dim))';
                P{j} = P{j} * sqrt(dim/trace(P{j}*P{j}'));
            end
            V = Q * blkdiag(P{:});
        end
        
        function blockMask = coarseGrainMask(mask, blockIndices)
        % Given a logical matrix of size nxn (mask) and a partition of 1:n (blockIndices)
        % returns the blockMask such that blockMask(r,c) is true if 
        % any of mask(blockIndices{r}, blockIndices{c})) is true
            nBlocks = length(blockIndices);
            blockMask = zeros(nBlocks, nBlocks); % mask is indexed using the eigenvalue "runs"
            for i = 1:nBlocks
                bi = blockIndices{i};
                for j = i:nBlocks
                    bj = blockIndices{j};
                    if any(any(mask(bi, bj)))
                        blockMask(i,j) = 1;
                        blockMask(j,i) = 1;
                    end
                end
            end
        end

        function orbits = findOrbits(group, settings)
        % Returns the orbits of the given group acting on the integers 1...n
        % as a partition (see DisjointSetForest.toPartition)
        %
        % If 'settings' is provided, and settings.blockDiagOrbits = false,
        % a single orbit {[1 ... n]} is returned and has the effect of disabling
        % orbit handling
            if nargin > 1 && ~settings.blockDiagOrbits
                orbits = {1:n}; % disable orbit lookup
                return
            end
            orbits = group.permOrbits.orbits;
        end
        
        function nonZeroRanges = findBlocks(mask)
        % Given a block-diagonal logical matrix (mask), returns the
        % ranges of nonzero blocks
        %
        % i.e. for [1 0 0; 0 0 0; 0 0 1] it returns {[1] [3]}
            d = size(mask, 1);
            from = 1;
            to = 1;
            nonZeroRanges = {};
            while from < d
                indices = sum(mask(from:to, :), 1) > 0;
                if ~any(indices)
                    from = to + 1;
                    to = to + 1;
                else
                    newTo = max(find(indices));
                    if newTo == to
                        nonZeroRanges{1, end+1} = from:to;
                        from = to + 1;
                        to = to + 1;
                    else
                        to = newTo;
                    end
                end
            end
        end

        function runs = findRuns(vector, settings)
        % Given a vector of ascending values, returns ranges of indices
        % [start(i) ... end(i)] such that vec(end(i))-vec(start(i)) <= settings.blockDiagEigTol
            runs = {};
            start = 1;
            for i = 1:length(vector) + 1
                if (i == length(vector) + 1) || vector(i) - vector(start) > settings.blockDiagEigTol
                    runs = horzcat(runs, {start:i-1});
                    start = i;
                end
            end
            if ~isequal(settings.blockDiagEigHist, [])
                for i = 1:length(runs)
                    run = runs{i};
                    diff = vector(run(2:end)) - vector(run(1:end-1));
                    settings.blockDiagEigHist.register(diff);
                end
            end
        end

        function connectedComponents = findConnectedComponents(adj)
        % Given an adjacency matrix adj, returns the sets of vertices corresponding
        % to connected components
        %
        % For adj = [0 0 1; 0 0 0; 1 0 0], it returns {[1 3] [2]}
            n = size(adj, 1);
            remaining = 1:n;
            cc = {};
            while ~isempty(remaining)
                comp = [];
                test = remaining(1);
                remaining = remaining(2:end);
                while ~isempty(test)
                    t = test(1);
                    test = test(2:end);
                    comp = [comp t];
                    newV = remaining(find(adj(t, remaining)));
                    test = [test newV];
                    remaining = setdiff(remaining, newV);
                end
                cc = horzcat(cc, {comp});
            end
            connectedComponents = cc;
        end

    end
end
