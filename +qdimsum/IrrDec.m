classdef IrrDec < handle
% Describes a decomposition of the group natural representation in irreducible representations
% over the reals.
%
% TODO: it identifies but does not handle quaternionic representations.
    
    properties (GetAccess = public, SetAccess = protected)
        group;        % Generalized permutation group of which we decompose the natural representation
        fromOrbit;    % fromOrbit(i) is the index of orbit in group.orbit from which
                      % the basis vector U(:,i) comes from
        U;            % Orthonormal change of basis matrix
        compDims;     % Isotypic component dimensions
        repDims;      % Representation dimensions
        repMuls;      % Representation multiplicities
                      % all types are real
        repTypes;     % 1 - real, 2 - complex, 3 - quaternionic
        settings;
    end
    
    methods
       
        function self = IrrDec(group, fromOrbit, U, repDims, repMuls, repTypes, settings)
            self.group = group;
            self.U = U;
            self.fromOrbit = fromOrbit;
            self.repDims = repDims(:)';
            self.repMuls = repMuls(:)';
            self.compDims = repDims(:)' .* repMuls(:)';
            self.repTypes = repTypes(:)';
            self.settings = settings;
            self.check;
        end
        
        function n = nComponents(self)
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

        function perm = swap(self)
        % Returns the permutation that switches between the kron(A, id) and kron(id, B) forms
            shift = 0;
            perm = [];
            for r = 1:self.nComponents
                d = self.repDims(r);
                m = self.repMuls(r);
                s = reshape(reshape(1:d*m, [d m])', [1 d*m]);
                perm = [perm, shift + s];
                shift = shift + d*m;
            end
        end
        
        function I = toIsoDec(self)
        % Returns the isotypic decomposition corresponding to this irreducible decomposition
            import qdimsum.*
            I = IsoDec(self.group, self.fromOrbit, self.U, true, self.repDims, self.repMuls, self.settings);
        end
        
        function M1 = project(self, M)
            M1 = self.U*self.projectInIrrepBasis(M, true)*self.U';
        end
        
        function M1 = projectInIrrepBasis(self, M, preserveCopies)
        % Projects the matrix M in the invariant subspace, changing the basis of M from the basis on
        % which the group acts to the basis that expresses the irreducible decomposition.
        %
        % preserveCopies determines whether we keep the copies of blocks when representations have dimension > 1
            if nargin < 3
                preserveCopies = true;
            end
            shift = 0;
            nC = self.nComponents;
            blocks = cell(1, nC);
            for r = 1:nC
                m = self.repMuls(r);
                d = self.repDims(r);
                switch self.repTypes(r)
                  case 1 % REAL TYPE
                    block = zeros(m, m);
                    for i = 1:d % Average over blocks
                        blockInd = shift + (i:d:m*d);
                        block = block + self.U(:,blockInd)'*M*self.U(:,blockInd);
                    end
                    block = block/d;
                    if preserveCopies
                        block = kron(block, eye(d));
                    end
                  case 2 % COMPLEX TYPE
                    % the block is made of 2x2 blocks of the form
                    % [r -i
                    %  i  r]
                    R = zeros(m, m);
                    I = zeros(m, m);
                    Rmask = eye(2);
                    Imask = [0 -1
                             1  0];
                    for i = 1:d/2 % Average over elements
                        blockInd = shift + (i:2*d:m*d);
                        R = R + self.U(:,blockInd)'*M*self.U(:,blockInd);
                        R = R + self.U(:,blockInd+1)'*M*self.U(:,blockInd+1);
                        I = I - self.U(:,blockInd)'*M*self.U(:,blockInd+1);
                        I = I + self.U(:,blockInd+1)'*M*self.U(:,blockInd);
                    end
                    R = R/d;
                    I = I/d;
                    if preserveCopies
                        R = kron(R, eye(d/2));
                        I = kron(I, eye(d/2));
                    end
                    block = kron(R, Rmask) + kron(I, Imask);
                  case 3 % QUATERNIONIC TYPE
                    % the block is made of 4x4 blocks of the form
                    % [a -b -c -d
                    %  b  a -d  c
                    %  c  d  a -b
                    %  d -c  b  a]
                    A = zeros(m, m);
                    B = zeros(m, m);
                    C = zeros(m, m);
                    D = zeros(m, m);
                    Amask = [1  0  0  0
                             0  1  0  0
                             0  0  1  0
                             0  0  0  1];
                    Bmask = [0 -1  0  0
                             1  0  0  0
                             0  0  0 -1
                             0  0  1  0];
                    Cmask = [0  0 -1  0
                             0  0  0  1
                             1  0  0  0
                             0 -1  0  0];
                    Dmask = [0  0  0 -1
                             0  0 -1  0
                             0  1  0  0
                             1  0  0  0];
                    for i = 1:d/4 % Average over elements
                        I = shift + (i:4*d:m*d);
                        A = A + self.U(:,I)'  *M*self.U(:,I);
                        A = A + self.U(:,I+1)'*M*self.U(:,I+1);
                        A = A + self.U(:,I+2)'*M*self.U(:,I+2);
                        A = A + self.U(:,I+3)'*M*self.U(:,I+3);
                        B = B - self.U(:,I)'  *M*self.U(:,I+1);
                        B = B + self.U(:,I+1)'*M*self.U(:,I);
                        B = B - self.U(:,I+2)'*M*self.U(:,I+3);
                        B = B + self.U(:,I+3)'*M*self.U(:,I+2);
                        C = C - self.U(:,I)'  *M*self.U(:,I+2);
                        C = C + self.U(:,I+1)'*M*self.U(:,I+3);
                        C = C + self.U(:,I+2)'*M*self.U(:,I);
                        C = C - self.U(:,I+3)'*M*self.U(:,I+1);
                        D = D - self.U(:,I)'  *M*self.U(:,I+3);
                        D = D - self.U(:,I+1)'*M*self.U(:,I+2);
                        D = D + self.U(:,I+2)'*M*self.U(:,I+1);
                        D = D + self.U(:,I+3)'*M*self.U(:,I);
                    end
                    A = A / d;
                    B = B / d;
                    C = C / d;
                    D = D / d;
                    if preserveCopies
                        A = kron(A, eye(d/4));
                        B = kron(B, eye(d/4));
                        C = kron(C, eye(d/4));
                        D = kron(D, eye(d/4));
                    end
                    block = kron(A, Amask) + kron(B, Bmask) + kron(C, Cmask) + kron(D, Dmask);
                end
                blocks{r} = block;
                shift = shift + m*d;
            end
            M1 = blkdiag(blocks{:});
        end
        
        function check(self)
        % Checks the validity of this irreducible decomposition
            import qdimsum.*
            % Perform isotypic checks
            self.toIsoDec.check;
            tol = self.settings.blockDiagMatTol;
            % Checks that the isotypic components are correct by considering
            % a sample from matrices that commute with the group
            sampleI = self.group.phaseConfiguration.sampleRealGaussian;
            sampleI = self.U'*sampleI*self.U;
            G1 = self.U'*GenPerm.orthogonalMatrix(self.group.randomElement)*self.U;
            G2 = self.U'*GenPerm.orthogonalMatrix(self.group.randomElement)*self.U;
            sampleG = randn * G1 + randn * G2;
            for r = 1:self.nComponents
                range = self.compRange(r);
                testI = sampleI(range, range);
                testG = sampleG(range, range);
                switch self.repTypes(r)
                  case 1
                    resI = IrrDec.projectRealBlocks(testI, self.repDims(r), false, true);
                    resG = IrrDec.projectRealBlocks(testG, self.repMuls(r), true, true);
                  case 2
                    resI = IrrDec.projectComplexBlocks(testI, self.repDims(r)/2, false, true);
                    resG = IrrDec.projectComplexBlocks(testG, self.repMuls(r), true, true);
                  case 3
                    resI = IrrDec.projectQuaternionicBlocks(testI, self.repDims(r)/4, false, true);
                    resG = IrrDec.projectQuaternionicBlocks(testG, self.repMuls(r), true, true);
                end
                assert(~isNonZeroMatrix(testI - resI, tol), 'Matrix that commutes with the group has not the proper form');
                assert(~isNonZeroMatrix(testG - resG, tol), 'Group algebra matrix has not the proper form');
            end
        end
    end

    methods (Static)
       
        function projected = projectRealBlocks(block, nCopies, collected, preserveCopies)
        % Projects a real type component of type I_nCopies (x) R(dR x dC)
        %
        % If collected = false, block is of the form kron(A, eye(nCopies))
        % If collected = true, block is of the form kron(eye(nCopies), A)
        %
        % preserveCopies = false  => keep a single copy of the block
        % preserveCopies = true   => restore the copies so that projected has the same size as block
            nR = size(block, 1);
            nC = size(block, 2);
            assert(mod(nR, nCopies) == 0);
            assert(mod(nC, nCopies) == 0);
            dR = nR/nCopies; % block dimensions
            dC = nC/nCopies;
            projected = zeros(dR, dC);
            for i = 1:nCopies
                if collected
                    IR = (i-1)*dR + (1:dR);
                    IC = (i-1)*dC + (1:dC);
                else
                    IR = i:nCopies:nR;
                    IC = i:nCopies:nC;
                end
                projected = projected + block(IR, IC);
            end
            projected = projected / nCopies;
            if preserveCopies
                if collected
                    projected = kron(eye(nCopies), projected);
                else
                    projected = kron(projected, eye(nCopies));
                end
            end
        end
        
        function projected = projectComplexBlocks(block, nCopies, collected, preserveCopies)
        % Project complex blocks
            nR = size(block, 1);
            nC = size(block, 2);
            assert(mod(nR, nCopies * 2) == 0);
            assert(mod(nC, nCopies * 2) == 0);
            dR = nR/nCopies/2; % complex dimension
            dC = nC/nCopies/2;
            R = zeros(dR, dC);
            I = zeros(dR, dC);
            for i = 1:nCopies
                if collected
                    IR = (i-1)*dR + (1:2:2*dR);
                    IC = (i-1)*dC + (1:2:2*dC);
                else
                    IR = i:2*nCopies:nR;
                    IC = i:2*nCopies:nC;
                end
                R = R + block(IR, IC) + block(IR+1, IC+1);
                I = I - block(IR, IC+1) + block(IR+1, IC);
            end
            R = R / (nCopies * 2);
            I = I / (nCopies * 2);
            if preserveCopies
                if collected
                    R = kron(eye(nCopies), R);
                    I = kron(eye(nCopies), I);
                else
                    R = kron(R, eye(nCopies));
                    I = kron(I, eye(nCopies));
                end
            end
            Rmask = [1  0
                     0  1];
            Imask = [0 -1
                     1  0];
            projected = kron(R, Rmask) + kron(I, Imask);
        end
        
        function projected = projectQuaternionicBlocks(block, nCopies, collected, preserveCopies)
            nR = size(block, 1);
            nC = size(block, 2);
            assert(mod(nR, nCopies * 4) == 0);
            assert(mod(nC, nCopies * 4) == 0);
            IR = 1:4:nR;
            IC = 1:4:nC;
            dR = nR/nCopies/4; % quaternionic dimension
            dC = nC/nCopies/4;
            A = zeros(dR, dC);
            B = zeros(dR, dC);
            C = zeros(dR, dC);
            D = zeros(dR, dC);
            for i = 1:nCopies
                if collected
                    IR = (i-1)*dR + (1:4:4*dR);
                    IC = (i-1)*dC + (1:4:4*dC);
                else
                    IR = i:4*nCopies:nR;
                    IC = i:4*nCopies:nC;
                end
                A =  block(IR   , IC    ) + block(IR + 1, IC + 1) + block(IR + 2, IC + 2) + block(IR + 3, IC + 3);
                B = -block(IR   , IC + 1) + block(IR + 1, IC    ) - block(IR + 2, IC + 3) + block(IR + 3, IC + 2);
                C = -block(IR   , IC + 2) + block(IR + 1, IC + 3) + block(IR + 2, IC    ) - block(IR + 3, IC + 1);
                D = -block(IR   , IC + 3) - block(IR + 1, IC + 2) + block(IR + 2, IC + 1) + block(IR + 3, IC    );
            end
            A = A/(nCopies * 4);
            B = B/(nCopies * 4);
            C = C/(nCopies * 4);
            D = D/(nCopies * 4);
            if preserveCopies
                if collected
                    A = kron(eye(nCopies), A);
                    B = kron(eye(nCopies), B);
                    C = kron(eye(nCopies), C);
                    D = kron(eye(nCopies), D);
                else
                    A = kron(A, eye(nCopies));
                    B = kron(B, eye(nCopies));
                    C = kron(C, eye(nCopies));
                    D = kron(D, eye(nCopies));
                end
            end
            Amask = [1  0  0  0
                     0  1  0  0
                     0  0  1  0
                     0  0  0  1];
            Bmask = [0 -1  0  0
                     1  0  0  0
                     0  0  0 -1
                     0  0  1  0];
            Cmask = [0  0 -1  0
                     0  0  0  1
                     1  0  0  0
                     0 -1  0  0];
            Dmask = [0  0  0 -1
                     0  0 -1  0
                     0  1  0  0
                     1  0  0  0];
            projected = kron(A, Amask) + kron(B, Bmask) + kron(C, Cmask) + kron(D, Dmask);
        end
        
    end
    methods (Static)
        
        function refinedBasis = orderedComplexBasis(iso, r)
        % Reorder the components of a complex representation so that the complex structure is visible
        % or returns [] if the representation is quaternionic
        % Optimized for precision
            import qdimsum.*
            tol = iso.settings.blockDiagEigTol;
            range = iso.compRange(r);
            d = iso.repDims(r);
            m = iso.repMuls(r);
            newU = iso.U;
            orbitsForRange = iso.fromOrbit(range);
            refinedBasis = zeros(iso.group.n, length(range));
            for o = unique(orbitsForRange) % iterate over all orbits present in this rep
                basisInd = range(orbitsForRange == o); % orbit indices for rep
                realRank = length(basisInd);
                % Euclidean coordinates of the o-th orbit elements (regardless of representation)
                oOrbit = iso.group.permOrbits.orbits{o};
                n = length(oOrbit);
                % find restriction of group to the o-th orbit
                resGroup = iso.group.permOrbitRestriction(o);
                % basis for the r-th representation in the o-th orbit
                basis = iso.U(oOrbit, basisInd);
                % compute a generic invariant sample (non-symmetric matrix), restricted to oOrbit x oOrbit
                sample = basis*Random.realGaussian(realRank)*basis';
                sample = resGroup.phaseConfiguration.project(sample); % project in the invariant subspace
                % Perform the Schur decomposition
                [U T] = schur(sample, 'real');
                % Get diagonal part = real part of eigenvalues
                D = diag(T);
                % Find non zero eigenvalues => part of the isotypic component
                nzInd = find(abs(D) > tol);
                % Find the zero eigenvalues => part of the null space
                zInd = find(abs(D) <= tol);
                % Cluster eigenvalues by reordering the Schur blocks
                clusters = zeros(1, n);
                % Put null space last
                clusters(zInd) = 1;
                % Find similar eigenvalues
                nzD = D(nzInd(1:2:end));
                distD = abs(bsxfun(@minus, nzD, nzD'));
                maskD = distD <= tol;
                conD = findConnectedComponents(maskD);
                % If the rank is not saturated, we have a quaternionic representation
                if length(conD) * 2 ~= realRank
                    % wrong rank found, it is a quaternionic representation
                    refinedBasis = [];
                    return
                end
                % Group similar eigenvalues
                for i = 1:length(conD)
                    compD1 = nzInd(conD{i}*2-1);
                    compD2 = nzInd(conD{i}*2);
                    assert(length(compD1) == d/2);
                    clusters(compD1) = i + 1;
                    clusters(compD2) = i + 1;
                end
                % Reorder
                [US TS] = ordschur(U, T, clusters);
                % Force blocks corresponding to the same complex eigenvalue to express it the same way
                % by removing the degeneracy due to conjugation
                start = 1;
                b = 1;
                for i = 1:m
                    b = b + 2;
                    for j = 2:n/m/2 % number of blocks
                        if sign(TS(b, b + 1)) == -sign(TS(start, start+1))
                            TS(:,[b b+1]) = TS(:,[b+1 b]);
                            TS([b b+1],:) = TS([b+1 b],:);
                            US(:,[b b+1]) = US(:,[b+1 b]);
                        end
                        b = b + 2;
                    end
                    start = b;
                end
                % Refined reordered basis is found
                refinedBasis(oOrbit, orbitsForRange == o) = US(:, 1:realRank);
            end
        end
        
    end
    
    methods (Static)
        
        function I = fromIsoDec(iso)
        % Constructs an irreducible decomposition from an isotypic decomposition
            import qdimsum.*
            sample = iso.group.phaseConfiguration.sampleRealGaussian;
            U = iso.U;
            repTypes = zeros(1, iso.nComponents);
            for r = 1:iso.nComponents
                d = iso.repDims(r);
                m = iso.repMuls(r);
                range = iso.compRange(r);
                if iso.repIsReal(r)
                    % use a second sample to put all irreducible components in the same basis
                    if iso.ordered
                        Urep = U(:, range);
                    else
                        Urep = iso.refinedBasis(r);
                    end
                    % only need the first row of blocks
                    repSample = Urep(:, 1:d)'*(sample + sample')*Urep;
                    % constructing the change of basis matrices
                    P = cell(1, m);
                    P{1} = eye(d);
                    for j = 2:m
                        P{j} = repSample(:, (j-1)*d+(1:d))';
                        P{j} = P{j} * sqrt(d/trace(P{j}*P{j}'));
                    end
                    Urep = Urep * blkdiag(P{:});
                    U(:, range) = Urep;
                    repTypes(r) = 1;
                else
                    Urep = IrrDec.orderedComplexBasis(iso, r);
                    if isequal(Urep, [])
                        % Quaternionic type
                        U(:, range) = NaN;
                        repTypes(r) = 3;
                    else
                        % Complex type
                        repSample = Urep(:, 1:d)'*sample*Urep;
                        P = cell(1, m);
                        P{1} = eye(d);
                        for j = 2:m
                            P{j} = repSample(:, (j-1)*d+(1:d))';
                            P{j} = P{j} * sqrt(d/trace(P{j}*P{j}'));
                        end
                        Urep = Urep * blkdiag(P{:});
                        U(:, range) = Urep;
                        repTypes(r) = 2;
                    end
                end
            end
            I = qdimsum.IrrDec(iso.group, iso.fromOrbit, U, iso.repDims, iso.repMuls, repTypes, iso.settings);
        end
        
        function I = forGroup(group, settings)
        % Constructs the irreducible decomposition of the given group
            import qdimsum.*
            if nargin < 2
                settings = NVSettings;
            end
            iso = IsoDec.forGroup(group);
            I = IrrDec.fromIsoDec(iso);
        end
                
    end
    
end
