classdef IsoDec < handle
% Isotypic decomposition of the natural action of a generalized permutation group
%
% The change of basis matrix is adapted to the orbits of the group
    
    properties (GetAccess = public, SetAccess = protected)
        group;        % Generalized permutation group of which we decompose the natural representation
        fromOrbit;    % fromOrbit(i) is the index of orbit in group.orbit from which
                      % the basis vector U(:,i) comes from
        U;            % Orthonormal change of basis matrix
        compDims;     % Isotypic component dimensions
        repDims;      % Representation dimensions
        repMuls;      % Representation multiplicities
        settings;
    end
    
    methods
        
        function n = nReps(self)
        % Number of isotypic components
            n = length(self.compDims);
        end
        
        function R = repRange(self, r)
        % Indices corresponding to the r-th representation
        % Correspond to columns of U, and to indices in fromOrbit
            from = sum(self.compDims(1:r-1)) + 1;
            to = sum(self.compDims(1:r));
            R = from:to;
        end
        
        function Urep = repBasis(self, r)
        % Returns the basis of the r-th isotypic components
            Urep = self.U(:, self.repRange(r));
        end
        
        function I = refine(self)
        % Refine representation bases
            import qdimsum.*
            U = self.U;
            for r = 1:self.nReps % refine each isotypic component
                range = self.repRange(r);
                orbits = self.fromOrbit(range); % orbits present in that component
                for o = 1:unique(orbits) % refine orbits individually
                    repOrbit = range(find(self.fromOrbit(range) == o)); % orbit indices for rep
                    n = length(repOrbit);
                    if n > 1
                        % the o-th orbit elements
                        fullOrbit = self.group.permOrbits.orbits{o};
                        % find restriction of group to the o-th orbit
                        resGroup = self.group.permOrbitRestriction(o);
                        % basis for the r-th representation in the o-th orbit
                        basis = U(fullOrbit, repOrbit);
                        % compute a sample
                        T = resGroup.phaseConfiguration.project(basis*Random.symmetricGaussian(n)*basis');
                        % force symmetry
                        T = T + T';
                        % compute eigenvalues, the n largest eigenvalues correspond to the representation basis
                        % TODO: would Gram-Schmidt be more precise?, or something else? EV decomposition seems
                        % wasteful
                        [refinedBasis, lambda] = eig(T);
                        lambda = diag(lambda);
                        [~, I] = sort(abs(lambda)); % sort eigenvalues (needed on Matlab 2015b)
                        I = flipud(I(:)); % largest magnitude first
                        lambda = lambda(I); % reorder accordingly
                        refinedBasis = refinedBasis(:, I);
                        refinedBasis = refinedBasis(:, 1:n); % cuts the possible additional eigenvectors
                        U(fullOrbit, repOrbit) = refinedBasis; % replace basis
                    end
                end
            end
            I = IsoDec(self.group, self.fromOrbit, U, self.repDims, self.repMuls, self.settings);
        end
        
        function o = smallestOrbitInRep(self, r)
            range = self.repRange(r);
            % need only to consider one orbit
            orbits = self.fromOrbit(range);
            O = unique(orbits(:));
            % Compute number of elements in the orbit
            O = [arrayfun(@(o) sum(orbits == o), O) O];
            O = sortrows(O);
            o = O(1, 2);
        end
        
        function r = repIsReal(self, r)
            import qdimsum.*
            % Eigenvalue tolerance
            tol = self.settings.blockDiagEigTol;
            % Find smallest group orbit to perform the test in
            o = smallestOrbitInRep(self, r);
            range = self.repRange(r);
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

        function self = IsoDec(group, fromOrbit, U, repDims, repMuls, settings)
        % Constructs an IsoDec from full data
            self.group = group;
            self.fromOrbit = fromOrbit;
            self.U = U;
            self.repDims = repDims(:)';
            self.repMuls = repMuls(:)';
            self.compDims = repDims(:)' .* repMuls(:)';
            self.settings = settings;
        end

    end
    
    methods (Static)
       
        function I = forGroup(group, settings)
        % Computes the isotypic decomposition of the given generalized permutation group
            [U reps fromOrbit] = qdimsum.Reps.isotypicComponents(group, settings);
            I = qdimsum.IsoDec(group, fromOrbit, U, reps(1, :), reps(2, :), settings);
        end
        
    end
    
end
