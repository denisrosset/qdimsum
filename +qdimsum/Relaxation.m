classdef Relaxation < handle
    
    properties (GetAccess = public, SetAccess = protected)
        problem = [];                 % NVProblem instance to solve
        monomials = {};               % Set of monomials to use
        operatorsGroup = [];          % Group acting on operator variables
        settings = [];                % NVSettings
    end
    
    properties (Access = protected)
        monomialsGroup_ = [];         % Group acting on monomials
    end
    
    methods
    
        function n = nMonomials(self)
        % Returns the number of elements in the monomials basis, i.e. the size of the moment matrix
            n = length(self.monomials);
        end
        
        function self = Relaxation(problem, monomials, settings)
        % Creates a relaxation for the given problem (NVProblem), with the given monomials
        % and NVSettings
        %
        % TODO: bring help text from nvOptimize
            self.problem = problem;
            self.monomials = monomials;
            self.operatorsGroup = qdimsum.Group(problem.symmetryGroupGenerators);
            self.settings = settings;
        end

        function h = monomialAction(self, g)
        % For a generalized permutation "g" on the operator variables, returns
        % the corresponding generalized permutation "h" on the monomials
            import qdimsum.*
            h = Monomials.findMonomialAction(self.problem, self.monomials.indices, g, self.settings);
        end
        
        function G = monomialsGroup(self)
        % Returns the group action on monomials
            import qdimsum.*
            if isequal(self.monomialsGroup_, [])
                self.monomialsGroup_ = self.operatorsGroup.monomorphism(@(g) self.monomialAction(g));
            end
            G = self.monomialsGroup_;
        end

        function [chi obj] = sample(self)
        % Returns a sample
            X = self.problem.sampleOperators;
            K = self.problem.sampleStateKraus;
            monos = self.monomials.computeKraus(X, K);
            stateDim = size(monos, 1);
            krausRank = size(monos, 2);
            monos = reshape(monos, [stateDim*krausRank self.nMonomials]);
            chi = monos'*monos; % compute moment matrix
            if self.problem.forceReal
                chi = real(chi);
            end
            chi = (chi + chi')/2; % to cater for possible rounding errors
            obj = self.problem.computeObjective(X, K);
        end
        
        function [chi obj] = symmetrizedSample(self)
        % Returns a sample symmetrized under the symmetry group
            import qdimsum.*
            [chi obj] = self.sample;
            chi = self.monomialsGroup.phaseConfiguration.project(chi);
            chi = (chi + chi')/2; % to cater for possible rounding errors
        end
        
        function [samples objs] = computeBasis(self)
        % Computes an affine basis of samples
            import qdimsum.*
            blockStructure = BlockStructure({1:self.nMonomials});
            sampleDim = blockStructure.dimension;
            samples = zeros(sampleDim, 0);
            objs = zeros(1, 0);
            chunkSize = self.settings.sampleChunkSize;
            l = 0;
            while 1
                % preallocate chunk
                samples = [samples zeros(sampleDim, chunkSize)];
                objs = [objs zeros(1, chunkSize)];
                for i = 1:chunkSize
                    [chi obj] = self.sample;
                    sample = BlockStructure.matToVec(chi);
                    l = l + 1;
                    samples(:,l) = sample;
                    objs(l) = obj;
                end
                l = l + chunkSize;
                r = rank(samples);
                if r < l
                    break
                end
            end
            samples = samples(:,1:r);
            objs = objs(1:r);
        end

        function [samples objs] = computeSymmetrizedBasis(self)
            import qdimsum.*
            blockStructure = BlockStructure({1:self.nMonomials});
            sampleDim = blockStructure.dimension;
            samples = zeros(sampleDim, 0);
            objs = zeros(1, 0);
            chunkSize = self.settings.sampleChunkSize;
            l = 0;
            while 1
                % preallocate chunk
                samples = [samples zeros(sampleDim, chunkSize)];
                objs = [objs zeros(1, chunkSize)];
                for i = 1:chunkSize
                    [chi obj] = self.symmetrizedSample;
                    sample = BlockStructure.matToVec(chi);
                    l = l + 1;
                    samples(:,l) = sample;
                    objs(l) = obj;
                end
                l = l + chunkSize;
                r = rank(samples);
                if r < l
                    break
                end
            end
            samples = samples(:,1:r);
            objs = objs(1:r);
        end
   
    end
end
