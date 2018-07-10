classdef CGLMP < NVProblem

    properties
        d; % number of outcomes
        forceReal = true;
    end
    
    methods
        
        function self = CGLMP(d)
            self.d = d;
        end
        
        function types = operatorTypes(self)
            d = self.d;
            types = {1:2*d 2*d+(1:2*d)}; % rho and M
        end
        
        function X = sampleOperators(self, rank)
            d=self.d;            
            rho = cell(d,2);
            for x = 1:2
                for  x0 = 1:d
                    r = (2*rand(d,1)-1) + 1i*(2*rand(d,1)-1);
                    rho{x0,x} = r*r';
                    rho{x0,x} = rho{x0,x}/norm(rho{x0,x});
                end
            end
            M = cell(d,2); % M(b,y)
            for y = 1:2
                U = qdimsum.Random.unitary(d);
                for b=1:d
                    M{b,y} = U*diag(full(sparse(1,b,1,1,d)))*U';
                end
            end
            X = self.pack(rho, M);
        end
        
        function K = sampleStateKraus(self)
            d = self.d;            
            K = eye(d);
        end
        
        function X = pack(self, rho, M)
            ind = 1;
            d = self.d;
            X = {};
            for x = 1:2
                for x0 = 1:d
                    X{ind} = rho{x0, x};
                    ind = ind + 1;
                end
            end
            for y = 1:2
                for b = 1:d
                    X{ind} = M{b,y};
                    ind = ind + 1;
                end
            end
        end
        
        function [rho M] = unpack(self, X)
            d = self.d;
            rho = cell(d, 2);
            M = cell(d, 2);
            ind = 1;
            for x = 1:2
                for x0 = 1:d
                    rho{x0, x} = X{ind};
                    ind = ind + 1;
                end
            end
            for y = 1:2
                for b = 1:d
                    M{b,y} = X{ind};
                    ind = ind + 1;
                end
            end
        end

        function obj = computeObjective(self, X, tau)
            d = self.d;
            [rho M] = self.unpack(X);
            obj = 0;
            for x = 0:1
                for x0 = 0:d-1
                    for y = 0:1
                        for K=0:floor(d/2)-1
                            b = mod(x0-x*y-(-1)^(x+y)*K,d);
                            obj = obj + (1-2*K/(d-1)) * trace(tau * rho{x0+1,x+1} * M{b+1,y+1});
                            b = mod(x0-x*y+(-1)^(x+y)*(K+1),d);
                            obj = obj - (1-2*K/(d-1))*trace(tau * rho{x0+1,x+1} * M{b+1,y+1});
                        end 
                    end
                end
            end
            obj = real(obj / (4*d));
        end
        
        function p = cyclicPermutation(self)
            d = self.d;
            cycle = [2:d 1];
            p = [cycle d+cycle 2*d+cycle 3*d+cycle];
        end
        
        function p = mirrorPermutation1(self)
            d = self.d;
            cycle = [d-1:-1:1 d];
            p = [2*d:-1:d+1 d:-1:1 3*d:-1:2*d+1 3*d+cycle];
        end
        
        function p = mirrorPermutation2(self)
            d = self.d;
            cycle = [1 d:-1:2];
            p = [d:-1:1 d+cycle 4*d:-1:2*d+1];
        end
        
        function generators = symmetryGroupGenerators(self)
            generators = [self.cyclicPermutation
                          self.mirrorPermutation1
                          self.mirrorPermutation2];
        end
    
    end
end
