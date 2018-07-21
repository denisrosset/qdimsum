classdef RAC22 < NVProblem

    properties
       forceReal = true;
    end
    
    methods
                
        function X = sampleOperators(self)
        % the first 4 operators are the states
        % the next 2 represent a projective measurement
        % the last 2 represent another projective measurement
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
        
        function types = operatorTypes(self)
            types = {1:4 5:8};
        end
        
        function K = sampleStateKraus(self)
            dim = 2;
            K = eye(dim);
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
                        obj = obj + trace(M * rho); % K is always identity
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
    
    methods
        
        function C = operatorSDPConstraints(self, X)
            C = X; % all operators are SDP
        end
        
        function C = operatorEqualityConstraints(self, X)
            dim = 2;
            C = {eye(dim) - X{5} - X{6}   % X{5}, X{6} form a projective measurement
                 eye(dim) - X{7} - X{8}}; % X{7}, X{8} form a projective measurement
        end
        
        function C = scalarEqualityConstraints(self, X)
            dim = 2;
            mm = eye(dim)/dim;
            C = {mm - X{1}
                 mm - X{2}
                 mm - X{3}
                 mm - X{4}};
        end

    end
    
end
