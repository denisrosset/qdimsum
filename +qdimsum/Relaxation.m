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

         function chi = sample(self)
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
        end
        
        function chi = sampleAndSymmetrize(self)
            import qdimsum.*
            chi = self.sample;
            chi = GenPerm.symmetrize(chi, self.monomialsGroup.decomposition);
        end
% $$$         
% $$$         function [M P] = phasePartition(self)
% $$$             C = qdimsum.PhaseConfiguration.fromGenPerm(length(self.monomials), self.monoGenerators);
% $$$             [M P] = C.toPhasePartition;
% $$$         end
% $$$         
% $$$         function c = javaPC(self)
% $$$             G = self.monoGenerators;
% $$$             G(G > 0) = G(G > 0) - 1;
% $$$             c = com.faacets.jqdimsum.PhaseConfiguration.fromGenPerm(length(self.monomials), G);
% $$$         end
% $$$         
% $$$         function chi = symmetrize(self, chi)
% $$$             chi = qdimsum.GenPerm.symmetrize(chi, self.monoAction);
% $$$             chi = (chi + chi')/2; % to cater for possible rounding errors
% $$$         end
% $$$         
    end
    
end
