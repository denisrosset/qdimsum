classdef RAC2d < NVProblem

    properties
        d;
        forceReal = true;
    end
    
    methods
        
        function self = RAC2d(d)
            self.d = d;
        end

        function types = operatorTypes(self)
            d = self.d;
            types = {1:d^2 d^2+(1:2*d)};
        end

        function X = sampleOperators(self, rank)
            d=self.d;            
            rho = cell(d, d);
            for x1 = 1:d
                for x2 = 1:d
                    r = (2*rand(d,1)-1) + 1i*(2*rand(d,1)-1);
                    rho{x1,x2} = r*r';
                    rho{x1,x2} = rho{x1,x2}/norm(rho{x1,x2});
                end
            end
            M = cell(d, 2); % M(b,y)
            for y = 1:2
                U = qdimsum.Random.unitary(d);
                for b=1:d
                    M{b,y} = U*diag(full(sparse(1,b,1,1,d)))*U';
                end
            end
            
            X={};
            for x2=1:d
                for x1=1:d
                    X = { X{:} rho{x1,x2} } ;
                end
            end
            for y=1:2
                for b=1:d
                    X={ X{:} M{b,y}};
                end
            end
            
        end
        
        function K = sampleStateKraus(self)
            d = self.d;            
            K = eye(d);
        end
        
        function X = pack(self, rho, M)
            ind = 1;
            d = self.d;
            X = cell(1, self.numOperatorVariables);
            for x2 = 1:d
                for x1 = 1:d
                    X{ind} = rho{x1, x2};
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
            rho = cell(d, d);
            M = cell(d, 2);
            ind = 1;
            for x2 = 1:d
                for x1 = 1:d
                    rho{x1, x2} = X{ind};
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
        
        function res = traceprod(A, B)
            A=A.';
            res = sum(A(:).*B(:));
        end

        function obj = computeObjective(self, X, tau)
            d = self.d;
            [rho M] = self.unpack(X);
            obj = 0;
            for x1 = 1:d
                for x2 = 1:d
                    for y = 1:2
                        if y == 1
                            b = x1;
                        else
                            b = x2;
                        end
                        A = rho{x1, x2};
                        B = M{b,y};
                        obj = obj + sum(conj(A(:)).*B(:)); % trace(rho{x1,x2} * M{b,y}); % don't need to include tau as it is identity
                    end
                end
            end
            obj = real(obj / (2*d^2));
        end
        
        function p = permutationX1(self, sigma)
            d = self.d;
            rho=reshape(1:d^2,[d d]);
            M=reshape(1:2*d,[d 2]);
            rho=reshape(rho(sigma,:),[1 d^2]);
            M(:,1)=M(sigma,1);
            M=reshape(M,[1 2*d])+d^2;
            p=[rho M];
        end
        
        function generators = symmetryGroupGenerators(self)
            d = self.d;
            id = 1:d^2+2*d;
            
            swapX1X2 = reshape(reshape(1:d^2, [d d])', d^2, 1)';
            Mv = reshape(1:2*d, [d 2]);
            M(:,1) = Mv(:,2);
            M(:,2) = Mv(:,1);
            
            swapX1X2 = [swapX1X2 reshape(M,[1 2*d])+d^2];         
            
            cyclicX1 = self.permutationX1([2:d 1]);
            swapX1 = self.permutationX1([2 1 3:d]);
            
            generators = [swapX1X2
                          cyclicX1
                          swapX1];
        end
        
        function generators = ambientGroupGenerators(self)
            d = self.d;
            cyclicX = [2:d^2 1 d^2+(1:2*d)];
            swapX = [2 1 3:d^2 d^2+(1:2*d)];
            swapY = [1:d^2 d^2+d+(1:d) d^2+(1:d)];
            cyclicB1 = [1:d^2 d^2+[2:d 1] d^2+d+(1:d)];
            swapB1 = [1:d^2 d^2+[2 1 3:d] d^2+d+(1:d)];
            generators = [cyclicX
                          swapX
                          swapY
                          cyclicB1
                          swapB1];
        end
        
    end
end
