classdef IsoDec < handle
% Isotypic decomposition of the natural action of a generalized permutation group
%
% The change of basis matrix is adapted to the orbits of the group
    
    properties (GetAccess = public, SetAccess = protected)
        group;        % Generalized permutation group of which we decompose the natural representation
        fromOrbit;    % fromOrbit(i) is the index of orbit in group.orbit from which
                      % the basis vector U(:,i) comes from
        U;            % Orthonormal change of basis matrix
        ordered;      % Whether the representations inside the isotypic components have already been ordered
                      % so that an element of the group algebra is block-diagonal in the isotypic component basis
                      % with m copies of size d x d, where m is the multiplicity and d the dimension
                      % However, this ordering does not mean that equivalent representations are expressed in the
                      % same basis, or that real/complex/quaternionic representations are recognized: that is the
                      % job of IrrDec.
        compDims;     % Isotypic component dimensions
        repDims;      % Representation (real) dimensions
        repMuls;      % Representation (real) multiplicities
        settings;
    end
    
    methods
        
        function self = IsoDec(group, fromOrbit, U, ordered, repDims, repMuls, settings)
        % Constructs an IsoDec from full data
            assert(isreal(U));
            self.group = group;
            self.fromOrbit = fromOrbit;
            self.U = U;
            self.ordered = ordered;
            self.repDims = repDims(:)';
            self.repMuls = repMuls(:)';
            self.compDims = repDims(:)' .* repMuls(:)';
            self.settings = settings;
            if settings.checkLevel > 0
                self.check;
            end
        end

        function n = nComponents(self)
        % Number of isotypic components
            n = length(self.compDims);
        end
        
        function R = compRange(self, r)
        % Indices corresponding to the r-th isotypic component
        % Correspond to columns of U, and to indices in fromOrbit
            from = sum(self.compDims(1:r-1)) + 1;
            to = sum(self.compDims(1:r));
            R = from:to;
        end
        
        function Urep = compBasis(self, r)
        % Returns the basis of the r-th isotypic component
            Urep = self.U(:, self.compRange(r));
        end
        
        function refinedBasis = refinedBasis(self, r)
        % Returns the refined basis elements for the r-th representation
            import qdimsum.*
            range = self.compRange(r); % basis indices
            orbits = self.fromOrbit(range); % orbits present in that component
            n = self.group.n;
            refinedBasis = zeros(n, length(range));
            for o = unique(orbits) % refine orbits individually
                basisInd = range(orbits == o); % basis elements we refine
                n = length(basisInd);
                % the o-th orbit elements
                oOrbit = self.group.permOrbits.orbits{o};
                % find restriction of group to the o-th orbit
                resGroup = self.group.permOrbitRestriction(o);
                % basis for the r-th representation in the o-th orbit
                basis = self.U(oOrbit, basisInd);
                % compute a sample
                T = basis*Random.symmetricGaussian(n)*basis'; % range in the corresponding representation
                T = resGroup.phaseConfiguration.project(T); % project in the invariant subspace
                                                            % compute eigenvalues, the n largest eigenvalues correspond to the representation basis
                T = T + T'; % force symmetry
                [U, ~] = sortedEig(T, 'descend', true);
                assert(isreal(U));
                refinedBasis(oOrbit, orbits == o) = U(:, 1:n); % replace basis cutting the possible additional eigenvectors
            end
        end
        
        function I = refine(self)
        % Refine the change of basis by performing a second step of eigenvalue decomposition inside each
        % isotypic component. As a bonus, it orders the irreducible representations inside the isotypic
        % components, so that I.ordered = true.
            import qdimsum.*
            U = self.U;
            for r = 1:self.nComponents % refine each isotypic component
                range = self.compRange(r);
                U(:, range) = self.refinedBasis(r);
            end
            I = IsoDec(self.group, self.fromOrbit, U, true, self.repDims, self.repMuls, self.settings);
        end
        
        function o = smallestOrbitInRep(self, r)
        % Returns the smallest orbit present in the r-th representation
            range = self.compRange(r);
            % need only to consider one orbit
            orbits = self.fromOrbit(range);
            O = unique(orbits(:));
            % Compute number of elements in the orbit
            O = [arrayfun(@(o) sum(orbits == o), O) O];
            O = sortrows(O);
            o = O(1, 2);
        end
        
        function r = repIsReal(self, r)
        % Returns whether the r-th representation is real
            import qdimsum.*
            % Eigenvalue tolerance
            tol = self.settings.blockDiagEigTol;
            % Find smallest group orbit to perform the test in
            o = self.smallestOrbitInRep(r);
            range = self.compRange(r);
            % Full orbit, can include other representations, used to select the
            % phase configuration to be sampled
            fullOrbit = self.group.permOrbits.orbits{o};
            % Representation basis vectors corresponding to that orbit
            repOrbit = range(self.fromOrbit(range) == o);
            Urep = self.U(fullOrbit, repOrbit);
            % Group restriction to the selected orbit
            resGroup = self.group.permOrbitRestriction(o);
            % Compute sample, transform in the representation space
            sampleGen = Urep'*resGroup.phaseConfiguration.sampleRealGaussian*Urep;
            sampleSym = sampleGen + sampleGen';
            % Compute eigenvalues of both the nonsymmetric and the made-symmetric matrix
            lambdaGen = eig(sampleGen);
            lambdaSym = eig(sampleSym);
            % Compute eigenvalues that are close
            distGen = abs(bsxfun(@minus, lambdaGen, lambdaGen.'));
            distSym = abs(bsxfun(@minus, lambdaSym, lambdaSym.'));
            maskGen = distGen <= tol;
            maskSym = distSym <= tol;
            % Histogram update to check precision
            if ~isequal(self.settings.blockDiagEigHist, [])
                self.settings.blockDiagEigHist.register(distGen(maskGen));
                self.settings.blockDiagEigHist.register(distSym(maskSym));
            end
            % Connect close eigenvalues
            conGen = Reps.findConnectedComponents(maskGen);
            conSym = Reps.findConnectedComponents(maskSym);
            lenGen = cellfun(@(x) length(x), conGen);
            lenSym = cellfun(@(x) length(x), conSym);
            assert(all(lenGen == lenGen(1)));
            assert(all(lenSym == lenSym(1)));
            % Same number of distinct eigenvalues between nonsym and sym? Then the representation is real
            r = length(conGen) == length(conSym);
        end
        
        function blocks = projectInIsoBasis(self, M)
            blocks = cell(1, self.nComponents);
            for r = 1:self.nComponents
                range = self.compRange(r);
                blocks{r} = self.U(:, range)' * M * self.U(:, range);
            end
        end

        function check(self)
        % Checks the validity of this isotypic decomposition
            import qdimsum.*
            tol = self.settings.blockDiagMatTol;
            % Checks that the isotypic components are correct by considering
            % a sample from matrices that commute with the group
            sample = self.group.phaseConfiguration.sampleRealGaussian;
            sample = self.U'*sample*self.U;
            for i = 1:self.nComponents
                ir = self.compRange(i);
                for j = 1:self.nComponents
                    jr = self.compRange(j);
                    block = sample(ir, jr);
                    assert(isNonZeroMatrix(block, tol) == (i == j));
                end
            end
            % Second check by using sampling from the group algebra
            M1 = self.U'*GenPerm.orthogonalMatrix(self.group.randomElement)*self.U;
            M2 = self.U'*GenPerm.orthogonalMatrix(self.group.randomElement)*self.U;
            M = randn * M1 + randn * M2;
            for i = 1:self.nComponents
                ir = self.compRange(i);
                for j = 1:self.nComponents
                    jr = self.compRange(j);
                    % standard check
                    block = M(ir, jr);
                    assert(isNonZeroMatrix(block, tol) == (i == j));
                    if i == j && self.ordered
                        % verify that irreducible representations are grouped correctly inside the
                        % isotypic component
                        m = self.repMuls(i);
                        d = self.repDims(i);
                        for k = 1:m
                            kr = d*(k-1) + (1:d);
                            for l = 1:m
                                lr = d*(l-1) + (1:d);
                                assert(isNonZeroMatrix(block(kr, lr), tol) == (k == l));
                            end
                        end 
                    end
                end
            end
        end
        
    end

    methods (Static)
       
        function I = forGroup(group, settings)
        % Computes the isotypic decomposition of the given generalized permutation group
            import qdimsum.*
            if nargin < 2
                settings = NVSettings; % use default settings if not provided
            end
            % Get problem structure
            n = group.phaseConfiguration.n;
            orbits = group.permOrbits.orbits;
            nOrbits = length(orbits);
            % Get first sample
            sample1 = group.phaseConfiguration.sampleSymmetricGaussian;
            % Data to be prepared
            fromOrbit = zeros(1, n);
            runs = {}; % identify blocks of repeated eigenvalues
            U = zeros(n, n); % use dense matrix for now, but should switch to sparse
                             % if the savings due to orbits are worth it
            shift = 0;
            % Treat each orbit individually, to preserve some sparsity
            for o = 1:nOrbits
                orbit = orbits{o}; % indices in the current orbit
                orbitSize = length(orbit);
                basisIndices = shift + (1:orbitSize); % indices of basis elements to compute
                [Uo Do] = sortedEig(sample1(orbit, orbit), 'ascend', false);
                Do = diag(Do);
                U(orbit, basisIndices) = Uo;
                fromOrbit(basisIndices) = o;
                runso = findRuns(Do, settings); % find subspaces corresponding to repeated eigenvalues
                runso = cellfun(@(r) r + shift, runso, 'UniformOutput', false); % shift to cater to basis indices
                runs = horzcat(runs, runso); % concatenate runs
                shift = shift + orbitSize;
            end
            % Now, U provides a basis that splits isotypic components. Remains to group them
            % according to their equivalent representations.
            sample2 = group.phaseConfiguration.sampleSymmetricGaussian;
            sample2p = U'*sample2*U;
            % We are computing the block mask, where each block corresponds to a run of identical
            % eigenvalues.
            % Blocks corresponding to inequivalent representations should be zero; to cater
            % for numerical errors, we check whether the matrix 2-norm is above or below
            % the tolerance 'blockDiagMatTol'.
            nRuns = length(runs);
            blockMask = logical(zeros(nRuns, nRuns));
            v = zeros(nRuns, nRuns);
            tol = settings.blockDiagMatTol;
            registerHist = ~isequal(settings.blockDiagMatHist, []);
            for i = 1:nRuns
                for j = 1:nRuns
                    block = sample2p(runs{i}, runs{j});
                    blockMask(i, j) = isNonZeroMatrix(block, tol);
                    if ~blockMask(i, j) && registerHist
                        % register all singular values of the block
                        singvals = svd(block);
                        settings.blockDiagMatHist.register(singvals(singvals <= tol));
                    end
                end
            end
            % find the subspaces corresponding to the same irreducible representation
            % by looking at the connected components of the graph defined by the adjacency
            % matrix of the block mask
            cc = findConnectedComponents(blockMask);
            Nc = length(cc);
            reps = zeros(2, Nc);
            for i = 1:Nc
                c = cc{i};
                reps(2, i) = length(c);
                dims = arrayfun(@(i) length(runs{i}), c);
                % verify that the dimensions are all the same for consistency
                assert(all(dims - dims(1) == 0), 'Inconsistent representation dimensions');
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
            I = qdimsum.IsoDec(group, fromOrbit, U, true, reps(1, :), reps(2, :), settings);
        end
        
    end
    
end
