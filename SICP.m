classdef SICP < NVProblem

    properties
        d;
        N;
        forceReal = true;
        rankVec = [];
        Nys = [];
        ys = [];
    end
    
    methods
        
        function self = SICP(N, d, rankVec)
            self.d = d; % dimension
            self.N = N; % problem size
            self.rankVec = rankVec;
            % y1 < y2
            self.Nys = N*(N-1)/2;
            self.ys = SICP.computeY(N);
        end

        function [y flip] = findY(self, y1, y2)
            assert(y1 ~= y2);
            flip = y1 > y2;
            if flip
                y = findrows(self.ys, [y2 y1]);
            else
                y = findrows(self.ys, [y1 y2]);
            end
            assert(length(y) == 1);
        end
        
        function types = operatorTypes(self)
            d = self.d;
            N = self.N;
            types = {1:N N+(1:self.Nys*2)};
        end

        function X = sampleOperators(self)
            import qdimsum.*
            d = self.d;
            N = self.N;
            rho = cell(1, N);
            for x = 1:N
                r = (2*rand(d,1)-1) + 1i*(2*rand(d,1)-1);
                rho{x} = r*r';
                rho{x} = rho{x}/norm(rho{x});
            end
            M = cell(2, self.Nys);
            for i = 1:self.Nys
                r1 = self.rankVec(N + (i-1)*2 + 1);
                r2 = self.rankVec(N + (i-1)*2 + 2);
                P = Random.projectivePOVM([r1 r2]);
                M{1,i} = P{1};
                M{2,i} = P{2};
            end
            X = horzcat(rho, M(:)');
        end
        
        function K = sampleStateKraus(self)
            d = self.d;            
            K = eye(d);
        end
        
        function obj = computeObjective(self, X, K)
            d = self.d;
            N = self.N;
            rho = X(1:N);
            M0 = X(N+1:2:end);
            M1 = X(N+2:2:end);
            obj = 0;
            for i = 1:self.Nys
                y1 = self.ys(i, 1);
                y2 = self.ys(i, 2);
                obj = obj + trace(K * rho{y1} * M0{i});
                obj = obj + trace(K * rho{y2} * M1{i});
            end
            obj = real(obj);
        end
       
        function generators = symmetryGroupGenerators(self)
            generators = self.computeGenerators;
        end
        
    end
    
    methods (Static)
        
        function ys = computeY(N)
            Nys = N*(N-1)/2;
            ys = zeros(Nys, 2);
            i = 1;
            for y2 = 2:N
                for y1 = 1:y2-1
                    ys(i, 1) = y1;
                    ys(i, 2) = y2;
                    i = i + 1;
                end
            end
        end
        
    end
    
    methods
        
        function p = permutationN(self, omega)
            %% Permute (y1 y2) with [omega(y1) omega(y2)] while sending x --> omega(x) and compensate order flips in (omega(y1), omega(y2))  by letting b --> b+1 .
            d = self.d;
            N = self.N;
            Nys = self.Nys;
            p = zeros(1, N + 2*Nys);
            p(1:N) = omega;
            for y = 1:Nys
                y1 = self.ys(y, 1);
                y2 = self.ys(y, 2);
                [newy flip] = self.findY(omega(y1), omega(y2));
                if flip
                    p(N + (y-1)*2 + 1) = N + (newy-1)*2 + 2;
                    p(N + (y-1)*2 + 2) = N + (newy-1)*2 + 1;
                else
                    p(N + (y-1)*2 + 1) = N + (newy-1)*2 + 1;
                    p(N + (y-1)*2 + 2) = N + (newy-1)*2 + 2;
                end
            end
        end
        
        function generators = computeGenerators(self, N)
        % Note that cyclicN and swapN are identical when N=2
            N = self.N;
            cyclicN = self.permutationN([2:N 1]);
            swapN = self.permutationN([2 1 3:N]); 
            generators = [cyclicN
                          swapN];
        end
        
    end
end
