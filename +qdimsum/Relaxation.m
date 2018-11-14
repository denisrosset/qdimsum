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
    
        function self = Relaxation(problem, userMonomials, settings)
            import qdimsum.*
            switch userMonomials{1}
              case 'npa'
                level = userMonomials{2};
                settings.log(sprintf('Computing monomials for NPA level %d', level));
                monomials = Monomials.npa(problem, level, settings);
              case 'families'
                families = userMonomials(2:end);
                familiesStr = cellfun(@(f) ['[' num2str(f) '],'], families, 'UniformOutput', false);
                familiesStr = strcat(familiesStr{:});
                familiesStr = familiesStr(1:end-1);
                settings.log(['Computing monomials for families ' familiesStr]);
                monomials = Monomials.families(problem, families, settings);
              otherwise
                monomials = userMonomials;
            end
            self.problem = problem;
            self.monomials = monomials;
            self.operatorsGroup = Group(problem.symmetryGroupGenerators);
            self.settings = settings;
        end

        function G = monomialsGroup(self)
            import qdimsum.*
            if isequal(self.monomialsGroup_, [])
                f = @(g) Monomials.findMonomialAction(self.problem, ...
                                                      self.monomials, ...
                                                      g, ...
                                                      self.settings);
                self.monomialsGroup_ = self.operatorsGroup.monomorphism(f);
            end
            G = self.monomialsGroup_;
        end

        function [chi obj] = sample(self)
            X = self.problem.sampleOperators;
            K = self.problem.sampleStateKraus;
            stateDim = size(K, 1);
            krausRank = size(K, 2);
            % compute monomials
            nMonomials = length(self.monomials);
            monos = zeros(stateDim, krausRank, nMonomials);
            for j = 1:nMonomials
                mono = K;
                indices = self.monomials{j};
                for k = length(indices):-1:1
                    mono = X{indices(k)} * mono;
                end
                monos(:,:,j) = mono;
            end
            monos = reshape(monos, [stateDim*krausRank nMonomials]);
            chi = monos'*monos; % compute moment matrix
            if self.problem.forceReal
                chi = real(chi);
            end
            chi = (chi + chi')/2; % to cater for possible rounding errors
            obj = self.problem.computeObjective(X, K);
        end
        
        function [chi obj] = symmetrizedSample(self)
            import qdimsum.*
            [chi obj] = self.sample;
            %chi1 = GenPerm.symmetrize(chi, self.monomialsGroup.decomposition);
            %chi1 = (chi1 + chi1')/2; % to cater for possible rounding errors
            chi = self.monomialsGroup.phaseConfiguration.project(chi);
            chi = (chi + chi')/2; % to cater for possible rounding errors
        end
   
    end
end
