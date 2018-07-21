classdef RAC22c < NVProblem
    properties
       forceReal = true;
    end
    methods
        function X = sampleOperators(self)
        % the first 4 operators are the states, next 2 a projective measurement
        % the last 2 another projective measurement
            dim = 2;
            X = cell(1, 8);
            for i = 1:4
                X{i} = qdimsum.Random.pureNormalizedDensityMatrix(dim);
            end
            U = qdimsum.Random.unitary(2);
            X{5} = U*[1 0; 0 0]*U';
            X{6} = U*[0 0; 0 1]*U';
            U = qdimsum.Random.unitary(2);
            X{7} = U*[1 0; 0 0]*U';
            X{8} = U*[0 0; 0 1]*U';
        end
        function K = sampleStateKraus(self)
            K = eye(2); % dimension is 2
        end
        function obj = computeObjective(self, X, K)
            obj = 0;
            for x1 = 1:2
                for x2 = 1:2
                    for y = 1:2
                        if y == 1
                            b = x1;
                        else
                            b = x2;
                        end
                        rho = X{x1+(x2-1)*2};
                        M = X{4+b+(y-1)*2};
                        obj = obj + real(trace(M * rho)/8);
                    end
                end
            end
        end
        function generators = symmetryGroupGenerators(self)
            swapX1X2 = [1 3 2 4 7 8 5 6];
            flipX1 = [2 1 4 3 6 5 7 8];
            generators = [swapX1X2
                          flipX1];
        end        
        function types = operatorTypes(self)
            types = {1:4 5:8};
        end
    end
end
