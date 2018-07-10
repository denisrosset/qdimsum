classdef RAC22 < NVProblem

    properties
       forceReal = true;
    end
    
    methods
        
        function C = operatorSDPConstraints(self, X)
            C = X; % all operators are SDP
        end
        
        function C = operatorEqualityConstraints(self, X)
            dim = 2;
            C = {eye(dim) - X{5} - X{6}
                 eye(dim) - X{7} - X{8}};
        end
        
        function C = scalarEqualityConstraints(self, X)
            dim = 2;
            mm = eye(dim)/dim;
            C = {mm - X{1}
                 mm - X{2}
                 mm - X{3}
                 mm - X{4}};
        end
        
        function X = sampleOperators(self, rank)
            dim = 2;
            rho = cell(2, 2);
            for x1 = 1:2
                for x2 = 1:2
                rho{x1,x2} = qdimsum.Random.pureNormalizedDensityMatrix(dim);
                end
            end
            M = cell(2, 2); % M(b,y)
            for y = 1:2
                U = qdimsum.Random.unitary(2);
                M{1,y} = U*[1 0; 0 0]*U';
                M{2,y} = U*[0 0; 0 1]*U';
            end
            X = self.pack(rho, M);
        end
        
        function types = operatorTypes(self)
            types = {1:4 5:8};
        end
        
        function X = pack(self, rho, M)
            X = {rho{1,1} rho{2,1} rho{1,2} rho{2,2} ...
                 M{1,1} M{2,1} M{1,2} M{2,2}};            
        end
        
        function [rho M] = unpack(self, X)
            rho = cell(2, 2);
            M = cell(2, 2);
            rho{1,1} = X{1};
            rho{2,1} = X{2};
            rho{1,2} = X{3};
            rho{2,2} = X{4};
            M{1,1} = X{5};
            M{2,1} = X{6};
            M{1,2} = X{7};
            M{2,2} = X{8};
        end
        
        function K = sampleStateKraus(self)
            dim = 2;
            K = eye(dim);
        end
        
        function obj = computeObjective(self, X, tau)
            [rho M] = self.unpack(X);
            obj = 0;
            for x1 = 1:2
                for x2 = 1:2
                    for y = 1:2
                        if y == 1
                            b = x1;
                        else
                            b = x2;
                        end
                        obj = obj + trace(tau * rho{x1,x2} * M{b,y});
                    end
                end
            end
            obj = real(obj / 8);
        end
        
        function generators = symmetryGroupGenerators(self)
            swapX1X2 = [1 3 2 4 7 8 5 6];
            flipX1 = [2 1 4 3 6 5 7 8];
            generators = [swapX1X2
                          flipX1];
        end
        
        function generators = ambientGroupGenerators(self)
            generators = [2 3 4 1 5 6 7 8
                          2 1 3 4 5 6 7 8
                          1 2 3 4 6 5 7 8
                          1 2 3 4 7 8 5 6];
        end
        
    end
end
