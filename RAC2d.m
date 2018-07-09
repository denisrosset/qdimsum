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
                U = Random.unitary(d);
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
                        obj = obj + trace(tau * rho{x1,x2} * M{b,y});
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
        
        function p = permutationX2(self, sigma)
            d = self.d;
            rho=reshape(1:d^2,[d d]);
            M=reshape(1:2*d,[d 2]);
            rho=reshape(rho(:,sigma),[1 d^2]);
            M(:,2)=M(sigma,2);
            M=reshape(M,[1 2*d])+d^2;
            p=[rho M];    
        end
        
        function chain = groupDecomposition(self)
            d = self.d;
            id = 1:d^2+2*d;
            
            swapX1X2=reshape(reshape(1:d^2, [d d])', d^2, 1)';
            Mv=reshape(1:2*d, [d 2]);
            M(:,1)=Mv(:,2);
            M(:,2)=Mv(:,1);
            
            swapX1X2=[swapX1X2 reshape(M,[1 2*d])+d^2];         
            swapX1X2;
            
            flipX1={};        
            for PermOrder=2:d    
                dummy=[];        
                permlist=toeplitz([1 PermOrder:-1:2], 1:PermOrder);
                for PermutationIndex=1:length(permlist)
                    sigma=permlist(PermutationIndex,:);
                    sigma=[sigma [PermOrder+1:d]];
                    p=self.permutationX1(sigma);
                    dummy=[dummy; p];
                end
                flipX1 = {flipX1{:} dummy};
                dummy;
            end

            flipX2={};        
            for PermOrder=2:d    
                dummy=[];        
                permlist=toeplitz([1 PermOrder:-1:2], 1:PermOrder);
                for PermutationIndex=1:length(permlist)
                    sigma=permlist(PermutationIndex,:);
                    sigma=[sigma [PermOrder+1:d]];
                    p=self.permutationX2(sigma);
                    dummy=[dummy; p];
                end
                flipX2 = {flipX2{:} dummy};
                dummy;
            end
            
            chain = {[id; swapX1X2] flipX1{:} flipX2{:}};
        end
        
    end
end
