% Class that describes a family of optimization problems in noncommutative variables, where
% members of the family differ by rank constraints
classdef (Abstract) RankProblem
    
    properties (Abstract)
        operatorDimensions;   % 1xN vector providing the dimension of the operator variables
        includeDeterministic; % Whether to include the deterministic projectors in optimization
    end
    
    methods (Abstract)
        
        % Given a rank card, returns the corresponding NVProblem instance
        I = problemInstance(self, rankCard);

    end
    
    methods % Can/should be overriden
        
        function generators = symmetryGroupGenerators(self)
            generators = 1:self.numOperators;
        end
                
    end
    
    methods
       
        function n = numOperators(self)
            n = length(self.operatorDimensions);
        end
            
        function cards = allRankCards(self)
            import qdimsum.*
            cards = RankProblem.generateRankCards(self.operatorDimensions, self.includeDeterministic);
        end
        
        function img = rankCardAction(self, genPerm, rankCard)
            img = zeros(1, self.numOperators);
            for i = 1:self.numOperators
                ip = genPerm(i);
                if ip > 0
                    img(abs(ip)) = rankCard(i);
                else
                    img(abs(ip)) = self.operatorDimensions(i) - rankCard(i);
                end
            end
        end
        
        function orbits = rankCardOrbits(self, rankCards)
            remaining = rankCards;
            generators = self.symmetryGroupGenerators;
            orbits = {};
            while size(remaining, 1) > 0
                orbit = remaining(1, :);
                remaining = remaining(2:end, :);
                newElement = true;
                while newElement
                    newElement = false;
                    for i = 1:size(orbit, 1)
                        for j = 1:size(generators, 1)
                            el = self.rankCardAction(generators(j, :), orbit(i, :));
                            if length(findrows(orbit, el)) == 0
                                present = findrows(remaining, el);
                                if length(present) == 0
                                    warning(['Rank ' num2str(el) ' not present']);
                                end
                                remaining = remaining(setdiff(1:size(remaining, 1), present), :);
                                orbit = [orbit
                                         el];
                                newElement = true;
                            end
                        end
                    end
                end
                orbits = horzcat(orbits, {orbit});
                
            end
        end
        
        function [maxObj orbits objs errored] = optimize(self, monomials, method, settings)
            cards = self.allRankCards;
            orbits = self.rankCardOrbits(cards);
            nOps = length(self.operatorDimensions);
            nOrbits = length(orbits);
            objs = zeros(nOrbits, 1);
            errored = zeros(0, nOps);
            for i = 1:nOrbits
                orbit = orbits{i};
                card = orbit(1, :);
                try
                    problem = self.problemInstance(card);
                    objs(i) = nvOptimize(problem, monomials, method, settings);
                catch ME
                    ME
                    objs(i) = NaN;
                    errored = [errored
                               card];
                end
            end
            maxObj = max(objs);
        end
        
    end
    
    methods (Static)
       
        function cards = generateRankCards(dimensions, includeDeterministic)
            import qdimsum.*
            d = dimensions(1);
            if includeDeterministic
                ranks = (0:d)';
            else
                ranks = (1:d-1)';
            end
            if length(dimensions) == 1
                cards = ranks;
            else
                subCards = RankProblem.generateRankCards(dimensions(2:end), includeDeterministic);
                cards = zeros(0, length(dimensions));
                subN = size(subCards, 1);
                for i = 1:length(ranks)
                    cards = [cards
                             ranks(i)*ones(subN, 1) subCards];
                end
            end
        end
        
    end
    
end
