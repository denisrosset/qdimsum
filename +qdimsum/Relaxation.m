classdef Relaxation < handle
    
    properties (GetAccess = public, SetAccess = immutable)
        problem = [];                 % NVProblem instance to solve
        monomials = {};               % Set of monomials to use
        operatorsGroup = [];          % Group acting on operator variables
        settings = [];                % NVSettings
    end
    
    properties (Access = protected)
        monomialsGroup_ = [];         % Group acting on monomials
        monomialsIsoDec_ = [];
        monomialsIrrDec_ = [];
    end
    
    methods (Static)
        
        function R = npa(problem, level, settings)
            import qdimsum.*
            if nargin < 3
                settings = NVSettings;
            end
            monomials = Monomials.fromNPA(problem, level, settings);
            R = Relaxation(problem, monomials, settings);
        end
        
        function R = families(problem, families, settings)
            import qdimsum.*
            if nargin < 3
                settings = NVSettings;
            end
            monomials = Monomials.fromFamilies(problem, families, settings);
            R = Relaxation(problem, monomials, settings);
        end
        
    end
    
    methods
    
        function n = nMonomials(self)
        % Returns the number of elements in the monomials basis, i.e. the size of the moment matrix
            n = length(self.monomials);
        end
        
        function self = Relaxation(problem, monomials, settings, ...
                                   monomialsGroup, monomialsIsoDec, monomialsIrrDec)
        % Creates a relaxation for the given problem (NVProblem), with the given monomials and NVSettings
        %
        % moomials
        % TODO: bring help text from nvOptimize
            self.problem = problem;
            self.monomials = monomials;
            self.operatorsGroup = qdimsum.Group(problem.symmetryGroupGenerators);
            if nargin >= 6 && ~isequal(monomialsIrrDec, [])
                self.monomialsIrrDec_ = monomialsIrrDec;
            end
            if nargin >= 5 && ~isequal(monomialsIsoDec, [])
                self.monomialsIsoDec_ = monomialsIsoDec;
            end
            if nargin >= 4 && ~isequal(monomialsGroup, [])
                self.monomialsGroup_ = monomialsGroup;
            end            
            if nargin >= 3 && ~isequal(settings, [])
                self.settings = settings;
            else
                self.settings = NVSettings;
            end
        end
        
        function G = monomialsGroup(self)
        % Returns the group action on monomials
            import qdimsum.*
            if isequal(self.monomialsGroup_, [])
                self.monomialsGroup_ = self.operatorsGroup.monomorphism(@(g) self.monomials.action(g));
            end
            G = self.monomialsGroup_;
        end
        
        function I = monomialsIsoDec(self)
        % Returns the isotypic decomposition of the monomial action group
            import qdimsum.*
            if isequal(self.monomialsIsoDec_, [])
                if self.settings.blockDiagRefine
                    self.monomialsIsoDec_ = IsoDec.forGroup(self.monomialsGroup, self.settings).refine;
                else
                    self.monomialsIsoDec_ = IsoDec.forGroup(self.monomialsGroup, self.settings);
                end
            end
            I = self.monomialsIsoDec_;
        end
        
        function I = monomialsIrrDec(self)
        % Returns the irreducible decomposition of the monomial action group
            import qdimsum.*
            if isequal(self.monomialsIrrDec_, [])
                self.monomialsIrrDec_ = IrrDec.fromIsoDec(self.monomialsIsoDec);
            end
            I = self.monomialsIrrDec_;
        end

        function chi = momentMatrix(self, X, K)
        % Returns the moment matrix computed from operators X and state Kraus dec. K
            monos = self.monomials.computeKraus(X, K);
            stateDim = size(monos, 1);
            krausRank = size(monos, 2);
            monos = reshape(monos, [stateDim*krausRank self.nMonomials]);
            chi = monos'*monos; % compute moment matrix
            if self.problem.forceReal
                chi = real(chi);
            end
            chi = (chi + chi')/2; % to cater for possible rounding errors
        end
        
        function blocks = blockDiagMomentMatrix(self, method, X, K)
        % Returns a single copy for each block present in the symmetry adapted basis 
        % method can be 'isotypic' 'irreps' or 'blocks'/'fastest'
            blocks = {};
            switch method
              case {'blocks', 'fastest'}
                I = self.monomialsIrrDec;
                monos = self.monomials.computeKrausInBasis(X, K, I.U);
                monos = reshape(monos, [size(monos,1)*size(monos,2) size(monos,3)]);
                shift = 0;
                blocks = {};
                realDims = [1 2 4];
                for r = 1:self.monomialsIrrDec.nComponents
                    m = I.repMuls(r);
                    d = I.repDims(r);
                    f = realDims(I.repTypes(r));
                    block = zeros(m*f, m*f);
                    for i = 1:d
                        range = shift + (i:d/f:d*m);
                        block = block + monos(:,range)'*monos(:,range);
                    end
                    block = block/(d/f);
                    if self.problem.forceReal
                        block = real(block);
                    end
                    switch I.repTypes(r)
                      case 1
                        % do nothing
                      case 2
                        block = IrrDec.projectComplexBlocks(block, 1, false, false);
                      case 3
                        block = IrrDec.projectQuaternionicBlocks(block, 1, false, false);
                    end
                    blocks{r} = block;
                    shift = shift + d*m;
                end
              case 'isotypic'
                chi = self.symmetrizedMomentMatrix(X, K);
                blocks = self.monomialsIsoDec.projectInIsoBasis(chi);
              case 'irreps'
                chi = self.symmetrizedMomentMatrix(X, K);
                blocks = self.monomialsIrrDec.projectInIrrBasis(chi, false, false); % invariant matrices, do not preserve copies
              otherwise
                error(['Unsupported method ' method]);
            end
            for r = 1:length(blocks)
                blocks{r} = (blocks{r} + blocks{r}')/2;
            end
        end

        function chi = symmetrizedMomentMatrix(self, X, K)
        % Equivalent to the "reynolds" method
            import qdimsum.*
            chi = self.momentMatrix(X, K);
            chi = self.monomialsGroup.phaseConfiguration.project(chi);
            chi = (chi + chi')/2; % to cater for possible rounding errors
        end

        function [vec bs] = computeBasisElement(self, method, X, K)
        % Returns a realization of a moment matrix basis element, vectorized
        % 
        % method can be none, reynolds, isotypic, irreps or blocks/fastest
        % Optional: provide the sample X (operators) and K (state Kraus decomposition)
            import qdimsum.*
            switch method
              case 'none'
                block = self.momentMatrix(X, K);
                blocks = {block};
                blockSizes = size(block, 1);
              case 'reynolds'
                block = self.symmetrizedMomentMatrix(X, K); 
                blocks = {block};
                blockSizes = size(block, 1);
              otherwise
                blocks = self.blockDiagMomentMatrix(method, X, K);
                blockSizes = cellfun(@(x) size(x, 1), blocks);
            end
            bs = BlockStructure(blockSizes);
            vec = bs.packBlocks(blocks);
        end

        function R = computeRealization(self, method)
        % Computes an affine basis of samples
            import qdimsum.*
            X = self.problem.sampleOperators;
            K = self.problem.sampleStateKraus;
            [vec blockSizes] = self.computeBasisElement(method, X, K);
            objs = self.problem.computeObjective(X, K);
            samples = vec(:);
            sampleDim = length(samples);
            chunkSize = self.settings.sampleChunkSize;
            l = 1;
            while 1
                % preallocate chunk
                samples = [samples zeros(sampleDim, chunkSize)];
                objs = [objs zeros(1, chunkSize)];
                for i = 1:chunkSize
                    X = self.problem.sampleOperators;
                    K = self.problem.sampleStateKraus;
                    l = l + 1;
                    [vec, ~] = self.computeBasisElement(method, X, K);
                    samples(:, l) = vec;
                    objs(l) = self.problem.computeObjective(X, K);
                end
                r = rank(samples);
                if r < l
                    break
                end
            end
            samples = samples(:,1:r);
            objs = objs(1:r);
            R = Realization(self, method, samples, blockSizes, objs);
        end

    end
end
