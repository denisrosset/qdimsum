classdef RankI3322c < RankProblem

    properties
        operatorDimensions = [];
        d = [];
        c = [];
        includeDeterministic = [];
    end
    
    methods
        
        function generators = symmetryGroupGenerators(self)
            id = 1:6;
            permPart = [4 5 6 1 2 3];
            flip1 = [2 1 3 4 5 -6];
            flip2 = [1 2 -3 5 4 6];
            generators = [permPart
                          flip1
                          flip2];
        end
        
        function self = RankI3322c(c, dimension, includeDeterministic)
            self.d = dimension;
            self.includeDeterministic = includeDeterministic;
            self.c = c;
            self.operatorDimensions = ones(1, 6) * dimension;
        end

        function problem = problemInstance(self, rankCard)
            problem = I3322c(rankCard, self.c, self.d);
        end
        
    end
        
end
