classdef SIC < NVProblem

    properties
        d;
        N;
        forceReal = true;
        rankVec = [];
        Nys = [];
        ys = [];
    end
    
    methods
        
        function self = SIC(N, d, rankVec)
            self.d = d; % dimension
            self.N = N; % problem size
            self.rankVec = rankVec;
            % y1 < y2
            self.Nys = N*(N-1)/2;
            self.ys = SIC.computeY(N);
        end

        function types = operatorTypes(self)
            d = self.d;
            N = self.N;
            types = {1:N N+(1:self.Nys)};
        end

        function X = sampleOperators(self)
            import qdimsum.*
            d = self.d;
            N = self.N;
            rho = cell(1, N);
            for x = 1:N
                rho{x} = pureNormalizedDensityMatrix(d);
            end
            M = cell(1, self.Nys);
            for i = 1:self.Nys
                r = self.rankVec(i + N);
                D = diag([ones(1, r) -ones(1,d-r)]);
                U = qdimsum.Random.unitary(d);
                M{i} = U*D*U';
            end
            X = horzcat(rho, M);
        end
        
        function K = sampleStateKraus(self)
            d = self.d;            
            K = eye(d);
        end
        
        function obj = computeObjective(self, X, K)
            d = self.d;
            N = self.N;
            rho = X(1:N);
            M = X(N+1:end);
            obj = 0;
            for i = 1:self.Nys
                y1 = self.ys(i, 1);
                y2 = self.ys(i, 2);
                obj = obj + trace(K * rho{y1} * M{i});
                obj = obj - trace(K * rho{y2} * M{i});
            end
            obj = real(obj);
        end
       
        function generators = symmetryGroupGenerators(self)
            generators = SIC.computeGenerators(self.N);
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
        
        function p = permutationN(omega)
        % Permute (y1 y2) with [omega(y1) omega(y2)] while sending x --> omega(x) and compensate order flips in (omega(y1), omega(y2))  by letting b --> b+1 .
            N = length(omega);
            Nys = N*(N-1)/2;
            ys = SIC.computeY(N);
            p = 1:(N+Nys);
            p(1:N) = omega;
            for i = 1:Nys
                y1 = ys(i, 1);
                y2 = ys(i, 2);
                y1p = omega(y1);
                y2p = omega(y2);
                if y2p > y1p
                    sign = 1;
                else
                    y1p = omega(y2);
                    y2p = omega(y1);
                    sign = -1;
                end
                ip = findrows(ys, [y1p y2p]);
                assert(length(ip) == 1);
                p(N+i) = (N+ip) * sign;
            end
        end
        
        function generators = computeGenerators(N)
        % Note that cyclicN and swapN are identical when N=2
            cyclicN = SIC.permutationN([2:N 1]);
            swapN = SIC.permutationN([2 1 3:N]); 
            generators=[cyclicN
                        swapN];
        end
        
    end
end
