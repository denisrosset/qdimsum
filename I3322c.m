classdef I3322c < NVProblem

    properties
        rankVec = [];
        c = [];
        d = [];
        forceReal = true;
    end
    
    methods
        
        function self = I3322c(rankVec, c, d)
            self.rankVec = rankVec;
            self.c = c;
            self.d = d;
        end
        
        function types = operatorTypes(self)
            types = {1:3 4:6};
        end

        function X = pack(self, A, B)
            X = {A{1} A{2} A{3} B{1} B{2} B{3}};            
        end
        
        function [A B] = unpack(self, X)
            A = cell(1, 3);
            B = cell(1, 3);
            A{1} = X{1};
            A{2} = X{2};
            A{3} = X{3};
            B{1} = X{4};
            B{2} = X{5};
            B{3} = X{6};
        end
        
        function C = implicitOperatorSDPConstraints(self, X)
            [A B] = self.unpack(X);
            d = self.d;
            id = eye(d*d);
            C = {id + A{1}
                 id - A{1}
                 id + A{2}
                 id - A{2}
                 id + A{3}
                 id - A{3}
                 id + B{1}
                 id - B{1}
                 id + B{2}
                 id - B{2}
                 id + B{3}
                 id - B{3}};
        end
                        
        function X = sampleOperators(self)
            rankVec = self.rankVec;
            d = self.d;
            for i = 1:6
                r = rankVec(i);
                D = diag([ones(1, r), -ones(1, d - r)]);
                U = qdimsum.Random.unitary(d);
                X{i} = U * D * U';
            end
            for i = 1:3
                X{i} = kron(X{i}, eye(d));
            end
            for i = 4:6
                X{i} = kron(eye(d), X{i});
            end
        end
        
        function psi = sampleStateKraus(self)
            d = self.d;
            psi = qdimsum.Random.normalizedPureState(d*d);
        end
        
        function obj = computeObjective(self, X, K)
            c = self.c;
            d = self.d;
            id = eye(d);
            [A B] = self.unpack(X);
            AB = cell(3,3);
            for x = 1:3
                for y = 1:3
                    AB{x,y} = A{x}*B{y};
                end
            end
            bellOperator = (A{1} + A{2} + B{1} + B{2} ...
                            -(AB{1,1} + AB{1,2} + AB{2,1} + AB{2,2}) ...
                            +(AB{1,3} + AB{3,1} - AB{2,3} - AB{3,2})*c);
            obj = real(K' * bellOperator * K);
        end
        
        function generators = symmetryGroupGenerators(self)
            id = 1:6;
            permPart = [4 5 6 1 2 3];
            flip1 = [2 1 3 4 5 -6];
            flip2 = [1 2 -3 5 4 6];
            generators = [permPart
                          flip1
                          flip2];
        end
        
    end
end
