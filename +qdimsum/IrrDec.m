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
        end
        
        function n = nComponents(self)
            n = length(self.compDims);
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
        
        function check(self)
        % Checks the validity of this irreducible decomposition
            import qdimsum.*
            % Perform isotypic checks
            self.toIsoDec.check;
            % TODO: irreducible checks
            
        end

    end

    methods (Static)
        
        function refinedBasis = orderedComplexBasis(iso, r)
        % Reorder the components of a complex representation so that the complex structure is visible
        % or returns [] if the representation is quaternionic
        % Optimized for precision
            import qdimsum.*
            tol = iso.settings.blockDiagEigTol;
            range = iso.repRange(r);
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
                sample = basis*Random.realGaussian(length(basisInd))*basis';
                sample = resGroup.phaseConfiguration.project(sample); % project in the invariant subspace
                [U T] = schur(sample);
                D = diag(T);
                nzInd = find(abs(D) > tol);
                zInd = find(abs(D) <= tol);
                clusters = zeros(1, n);
                clusters(zInd) = 1;
                nzD = D(nzInd(1:2:end));
                distD = abs(bsxfun(@minus, nzD, nzD'));
                maskD = distD <= tol;
                conD = findConnectedComponents(maskD);
                if length(conD) * 2 ~= realRank
                    % wrong rank found, it is a quaternionic representation
                    refinedBasis = [];
                    return
                end
                for i = 1:length(conD)
                    compD1 = nzInd(conD{i}*2-1);
                    compD2 = nzInd(conD{i}*2);
                    assert(length(compD1) == d/2);
                    clusters(compD1) = i + 1;
                    clusters(compD2) = i + 1;
                end
                [US TS] = ordschur(U, T, clusters);
                % force blocks corresponding to the same complex eigenvalue to express it
                % the same way
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
            repTypes = zeros(1, iso.nReps);
            for r = 1:iso.nReps
                d = iso.repDims(r);
                m = iso.repMuls(r);
                range = iso.repRange(r);
                if iso.repIsReal(r)
                    % use a second sample to put all irreducible components in the same basis
                    if iso.ordered
                        Urep = U(:, range);
                    else
                        Urep = iso.refinedRepresentationBasis(r);
                    end
                    % only need the first row of blocks
                    repSample = Urep(:, 1:d)'*(sample + sample')*Urep;
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
                        U(:, range) = NaN;
                        repTypes(r) = 3;
                    else
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
