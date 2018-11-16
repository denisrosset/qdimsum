classdef IrrDec < handle
% Describes a decomposition of the group natural representation in irreducible representations
% over the reals.
%
% TODO: it identifies but does not handle quaternionic representations.
    
    properties (GetAccess = public, SetAccess = protected)
        group;        % Generalized permutation group of which we decompose the natural representation
        U;            % Orthonormal change of basis matrix
        compDims;     % Isotypic component dimensions
        repDims;      % Representation dimensions
        repMuls;      % Representation multiplicities
                      % all types are real
        repTypes;     % 1 - real, 2 - complex, 3 - quaternionic
    end
    
    methods
       
        function self = IrrDec(group, U, repDims, repMuls, repTypes)
            self.group = group;
            self.U = U;
            self.repDims = repDims(:)';
            self.repMuls = repMuls(:)';
            self.compDims = repDims(:)' .* repMuls(:)';
            self.repTypes = repTypes(:)';
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

    end
    
    methods (Static)
       
        function [V newFromOrbit] = decomposeRealComponent(self, sample1, sample2, dim, mult, fromOrbit)
        % Given two (symmetric matrix) samples from an isotypic component of real type for a representation
        % of the given dimension and multiplicity, computes the change of basis matrix that puts all copies
        % of the irreducible representation in the same basis.
            n = dim * mult;
            Q = [];
            newFromOrbit = [];
            for o = unique(fromOrbit)
                % for each orbit
                ind = find(fromOrbit == o);
                C = zeros(n, length(ind));
                % group eigenvalues
                [Qo, lambda] = eig(sample1(ind, ind));
                [~, I] = sort(abs(diag(lambda))); % sort eigenvalues (needed on Matlab 2015b)
                I = flipud(I(:)); % largest magnitude first
                Qo = Qo(:, I);
                C(ind, :) = Qo;
                Q = [Q C];
                newFromOrbit = [newFromOrbit o*ones(1, length(n))];
            end
            % use a second sample to put all irreducible components in the same basis
            T = Q(:,1:dim)'*sample2*Q;
            P = cell(1, mult);
            P{1} = eye(dim);
            for j = 2:mult
                P{j} = T(:, (j-1)*dim+(1:dim))';
                P{j} = P{j} * sqrt(dim/trace(P{j}*P{j}'));
            end
            V = Q * blkdiag(P{:});
        end
        
        function [t Urep] = splitComponent(self, r)
            import qdimsum.*
            tol = self.settings.blockDiagEigTol;
            Urep = self.repBasis(r);
            sampleGen = Urep'*self.group.phaseConfiguration.sampleRealGaussian*Urep;
            sampleSym = sampleGen + sampleGen';
            lambdaGen = eig(sampleGen);
            lambdaSym = eig(sampleSym);
            distGen = abs(bsxfun(@minus, lambdaGen, lambdaGen.'));
            distSym = abs(bsxfun(@minus, lambdaSym, lambdaSym.'));
            maskGen = distGen <= tol;
            maskSym = distSym <= tol;
            if ~isequal(self.settings.blockDiagEigHist, [])
                self.settings.blockDiagEigHist.register(distGen(maskGen));
                self.settings.blockDiagEigHist.register(distSym(maskSym));
            end
            conGen = Reps.findConnectedComponents(maskGen);
            conSym = Reps.findConnectedComponents(maskSym);
            lenGen = cellfun(@(x) length(x), conGen);
            lenSym = cellfun(@(x) length(x), conSym);
            assert(all(lenGen == lenGen(1)));
            assert(all(lenSym == lenSym(1)));
            if length(conGen) == length(conSym)
                t = 1;
                Urep = [];
            else
                t = 2;
                [U T] = schur(sampleGen);
                D = diag(T);
                D = D(1:2:end);
                distD = abs(bsxfun(@minus, D, D'));
                maskD = distD <= tol;
                conD = Reps.findConnectedComponents(maskD);
                n = self.compDims(r);
                d = self.repDims(r);
                m = self.repMuls(r); % multiplicity = number of distinct eigenvalue pairs
                assert(length(conD) == m);
                clusters = zeros(1, n);
                for i = 1:length(conD)
                    compD = conD{i};
                    assert(length(compD) == d/2);
                    clusters((compD-1)*2+1) = i;
                    clusters((compD-1)*2+2) = i;
                end
                [US TS] = ordschur(U, T, clusters);
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
                norm(US*TS*US' - sampleGen)
                sampleGen2 = Urep'*self.group.phaseConfiguration.sampleRealGaussian*Urep;
                TT = US'*sampleGen2*US;
                P = cell(1, m);
                P{1} = eye(d);
                for j = 2:m
                    P{j} = TT(1:d, (j-1)*d+(1:d))';
                    P{j} = P{j} * sqrt(d/trace(P{j}*P{j}'));
                end
                PP = blkdiag(P{:});
                Urep  = Urep * US * PP;
                %size(Urep)
                %Urep = Urep * V;
                %size(sampleGen)
                %size(Urep)
                %Urep*sampleGen*Urep'
                %norm(US*TS*US' - sampleGen)
            end
% $$$                 nCopies
% $$$                 distinct = cellfun(@(ind) mean(lambdaGen(ind)), conGen);
% $$$                 
% $$$                 ordered = [];
% $$$                 while length(distinct) > 0
% $$$                     l = distinct(1);
% $$$                     distinct = distinct(2:end);
% $$$                     lc = conj(l);
% $$$                     distinct = distinct(find(abs(distinct - lc) > tol));
% $$$                     ordered = [ordered l lc];
% $$$                 end
% $$$                 ordered
                
                %distReal = abs(bsxfun(@minus, real(lambdaGen), real(lambdaGen)'));
                %maskReal = distReal <= self.settings.blockDiagEigTol;
                %conReal = Reps.findConnectedComponents(maskReal);
                
                %list = cell2mat(conGen')
                %d = self.compDims(r); % component dimension
            %r = rand;
            %[~, I] = sort(check);
            %[~, Isym] = sort(checkSym);
            %runs = Reps.findRuns(check(I), settings);
            %runsSym = Reps.findRuns(checkSym(Isym), settings);
            %if length(runs) == length(runsSym)
            %    type = 'R';
            %else
            %    type = 'CH';
            %end
        end
            
        end
        
        function I = fromIsoDec(isoDec)
            
        end
        
        function I = forGroup(group, settings)
            [U reps] = qdimsum.Reps.irreducibleDecomposition(group, settings);
            types = ones(1, size(reps, 2));
            I = qdimsum.IrrDec(group, U, reps(1, :), reps(2, :), types);
        end
        
    end
    
end
