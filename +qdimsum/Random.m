classdef Random
    
    methods (Static)
                
        function U = unitary(dim)
        % Returns a random unitary matrix of dimension dim
        % construct the Ginibre ensemble
            gin = randn(dim);
            gin = gin + 1i*randn(dim);
            % QR decomposition of the Ginibre ensemble
            [Q,R] = qr(gin);
            % compute U from the QR decomposition
            R = sign(diag(R));
            R(R==0) = 1; % protect against potentially zero diagonal entries
            U = bsxfun(@times,Q,R.'); % much faster than the naive U = Q*diag(R) 
        end
        
        function M = realGaussian(d)
            M = randn(d);
        end
        
        function M = complexGaussian(d)
            M = (randn(d) + randn(d) * 1i)/sqrt(2);
        end
        
        function M = hermitianGaussian(d)
        % Generates a Hermitian matrix with measure invariant under unitary transformations,
        % sampled from the Gaussian Unitary Ensemble
        % see http://staff.math.su.se/shapiro/UIUC/random_matrices.pdf
        %
        % d is the dimension            
            M = zeros(d, d);
            for r = 1:d
                M(r, r) = randn;
                M(r, r+1:end) = (randn(1, d-r) + randn(1, d-r)*1i)/sqrt(2);
                M(r+1:end, r) = conj(M(r, r+1:end));
            end
        end
        
        function M = symmetricGaussian(d)
        % Generates a symmetric matrix with measure invariant under orthogonal transformations,
        % sampled from the Gaussian Orthogonal Ensemble
        % see http://staff.math.su.se/shapiro/UIUC/random_matrices.pdf
        %
        % d is the dimension
            M = zeros(d, d);
            for r = 1:d
                % diagonal elements are scaled up by sqrt(2)
                M(r, r) = randn * sqrt(2);
                % while other elements are standard normals
                M(r, r+1:end) = randn(1, d-r);
                M(r+1:end, r) = M(r, r+1:end);
            end
        end
        
        function M = symmetricTraceOneGem(d)
            M = 2*gemRand(d, d) - 1;
            M = M + M';
            M = M/trace(M);
        end
        
        function M = realTraceOneGem(d)
            M = 2*gemRand(d, d) - 1;
            M = M/trace(M);
        end
        
        function M = symmetricTraceOne(d)
            M = 2*rand(d, d) - 1;
            M = M + M';
            M = M/trace(M);
        end
        
        function M = realTraceOne(d)
            M = 2*rand(d, d) - 1;
            M = M/trace(M);
        end
        
        function M = complexTraceOne(d)
            M = (2*rand(d, d) - 1) + 1i*(2*rand(d, d) - 1);
            M = M/trace(M);            
        end
        
        function M = hermitianTraceOne(d)
            M = (2*rand(d, d) - 1) + 1i*(2*rand(d, d) - 1);
            M = M + M';
            M = M/trace(M);
        end
        
        function M = psdTraceOne(d)
            M = (2*rand(d, d) - 1) + 1i*(2*rand(d, d) - 1);
            M = M*M';
            M = M/trace(M);
        end
        
        function ket = unnormalizedPureState(d)
            ket = (2*rand(d,1)-1) + 1i*(2*rand(d,1)-1);
        end
        
        function ket = normalizedPureState(d)
            ket = qdimsum.Random.unnormalizedPureState(d);
            ket = ket / norm(ket);
        end
        
        function rho = pureNormalizedDensityMatrix(d)
            ket = qdimsum.Random.unnormalizedPureState(d);
            rho = ket*ket';
            rho = rho / trace(rho);
        end
        
        function P = projectivePOVM(ranks)
            import qdimsum.*
            d = sum(ranks);
            U = Random.unitary(d);
            n = length(ranks);
            P = cell(1, n);
            shift = 0;
            for i = 1:n
                r = ranks(i);
                D = [zeros(1, shift) ones(1, r) zeros(1, d-shift-r)];
                P{i} = U*diag(D)*U';
                shift = shift + r;
            end
        end
        
    end
    
end
