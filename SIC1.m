classdef SIC1 < NVProblem

    properties
        d;
        N;
        mu;
        Nys;
        ys;
        component;
        forceReal = true;
    end
    
    methods
       
        function self = SIC1(N, d, mu, component)
            import qdimsum.*
            assert(N >= d);
            self.d = d;
            self.N = N;
            Nys = N*(N-1)/2;
            self.Nys = Nys;
            self.ys = zeros(Nys, 2);
            self.mu = mu;
            if nargin < 4
                d1 = floor(d/2);
                d2 = d - d1;
                component = [ones(1, Nys)*d1
                             ones(1, Nys)*d2];
                component = component(:)';
            end
            self.component = component;
            ind = 1;
            for y2 = 2:N
                for y1 = 1:y2-1
                    self.ys(ind, 1) = y1;
                    self.ys(ind, 2) = y2;
                    ind = ind + 1;
                end
            end
        end
        
        function types = operatorTypes(self)
            d = self.d;
            N = self.N;
            N1 = N;
            N2 = self.Nys*2;
            N3 = d;
            types = {1:N1 N1+(1:N2) (N1+N2)+(1:N3)};
        end
        
        function rho = sampleStates(self)
            import qdimsum.*
            N = self.N;
            rho = cell(1, N);
            for i = 1:N
                rho{i} = Random.pureNormalizedDensityMatrix(self.d);
            end
        end
        
        function B = sampleFirstMeasurement(self)
            import qdimsum.*
            d = self.d;
            Nys = self.Nys;
            component = reshape(self.component, [2 Nys]);
            B = cell(2, Nys);
            for i = 1:Nys
                B(:,i) = Random.projectivePOVM(component(:,i));
            end
            B = B(:)';
        end
        
        function C = sampleSecondMeasurement(self)
            import qdimsum.*
            d = self.d;
            N = self.N;
            C = Random.projectivePOVM(ones(1, d));
        end
                
        function X = sampleOperators(self)
            import qdimsum.*
            d = self.d;
            N = self.N;
            rho = self.sampleStates;
            B = self.sampleFirstMeasurement;
            C = self.sampleSecondMeasurement;
            X = horzcat(horzcat(rho, B), C);
        end

        function K = sampleStateKraus(self)
            K = eye(self.d);
        end
        
        function M = monomialsInObjective(self)
            M = {};
            d = self.d;
            N = self.N;
            Nys = self.Nys;
            ys = self.ys;
            mu = self.mu;
            for i = 1:Nys
                y1 = ys(i, 1);
                y2 = ys(i, 2);
                M = horzcat(M, {[y1 N+(i-1)*2+1] [y2 N+(i-1)*2+2]});
            end
        end
        
        function obj = computeObjective(self, X, K)
            d = self.d;
            N = self.N;
            Nys = self.Nys;
            ys = self.ys;
            mu = self.mu;
            N1 = N;
            N2 = 2*Nys;
            N3 = d;
            rho = X(1:N1);
            B = reshape(X(N1+(1:N2)), [2 Nys]);
            C = X(N1+N2+(1:N3));
            obj = 0;
            for i = 1:Nys
                y1 = ys(i, 1);
                y2 = ys(i, 2);
                obj = obj + trace(K' * rho{y1} * B{1,i} * K);
                obj = obj + trace(K' * rho{y2} * B{2,i} * K);
            end
            for i = 1:N3
                obj = obj + mu*trace(K' * rho{i} * C{i} * K);
            end
            obj = real(obj);
        end

        function p1 = firstMeasurementAction(self, omega1, omega2)
            import qdimsum.*
            Nys = self.Nys;
            ys = self.ys;
            g = zeros(1, Nys);
            omega = GenPerm.directSum(omega1, omega2);
            for i = 1:Nys
                y1 = self.ys(i, 1);
                y2 = self.ys(i, 2);
                y1p = omega(y1);
                y2p = omega(y2);
                if y2p > y1p
                    sign = 1;
                    ip = findrows(ys, [y1p y2p]);
                else
                    sign = -1;
                    ip = findrows(ys, [y2p y1p]);
                end
                g(i) = sign*ip;
            end
            p1 = GenPerm.toUnsigned(g);
        end
        
        function g = computeGenerator(self, omega1, omega2)
            import qdimsum.*
            O1 = GenPerm.directSum(omega1, omega2);
            O2 = self.firstMeasurementAction(omega1, omega2);
            O3 = omega1;
            g = GenPerm.directSum(O1, O2, O3);
        end
        
        function G = symmetryGroupGenerators(self)
            import qdimsum.*
            d = self.d;
            N = self.N;
            e = N - d;
            G  = [self.computeGenerator(1:d, [(2:e) 1])
                  self.computeGenerator([(2:d) 1], 1:e)];
            if d > 1
                G = [G
                     self.computeGenerator([2 1 (3:d)], 1:e)];
            end
            if e > 1
                G = [G
                     self.computeGenerator(1:d, [2 1 (3:e)])];
            end
        end
        
    end
        
end
