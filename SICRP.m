classdef SICRP < RankProblem

    properties
        d;XD
        G
        N;
        Nys;
        includeDeterministic = [];
        operatorDimensions;
    end
    
    methods
        
        function generators = symmetryGroupGenerators(self)
            generators = SIC.computeGenerators(self.N);
        end
        
        function self = SICRP(N, d, includeDeterministic)
            self.d = d;
            self.N = N;
            self.Nys = N*(N-1)/2;
            self.includeDeterministic = includeDeterministic;
            self.operatorDimensions = [2*ones(1, self.N) d*ones(1, self.Nys)];
        end

        function problem = problemInstance(self, rankCard)
            problem = SIC(self.N, self.d, rankCard);
        end
    end
end
