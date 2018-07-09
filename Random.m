classdef Random
    
    methods (Static)
        
        function g = groupElement(chain)
            g = chain{1};
            g = g(1,:);
            C = length(chain);
            for i = C:-1:1
                c = chain{i};
                r = randi(size(c, 1));
                g = GenPerm.compose(c(r,:), g);
            end
        end
        
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
        
        function M = realGaussian(d, sigma)
            if nargin < 2
                sigma = 1;
            end
            M = randn(d) * sigma;
        end
        
        function M = complexGaussian(d, sigma)
            if nargin < 2
                sigma = 1;
            end
            M = (randn(d) + randn(d) * 1i)*sigma/sqrt(2); % TODO: is that correct?
        end
        
        function M = hermitianGaussian(d, sigma)
            if nargin < 2
                sigma = 1;
            end
            M = zeros(d, d);
            for r = 1:d
                M(r, r:end) = randn(1, d-r+1) * (sigma / sqrt(2));
                M(r, r+1:end) = randn(1, d-r) * (1i*sigma / sqrt(2));
                M(r, r) = M(r, r) * sqrt(2);
                M(r+1:end, r) = conj(M(r, r+1:end));
            end
        end
        
        function M = symmetricGaussian(d, sigma)
        % Generates a symmetric matrix with measure invariant under orthogonal transformes
        % see http://reference.wolfram.com/language/ref/GaussianOrthogonalMatrixDistribution.html
        % however, that reference is ambiguous on how to sample the matrix elements
        % (should the diagonal elements be scaled by sqrt(2) or not? we choose the affirmative)
        %
        % d is the dimension and sigma the standard deviation
            if nargin < 2
                sigma = 1;
            end
            M = BlockStructure.vecToMat(randn(1, d*(d+1)) * (sigma / sqrt(2)), d);
            for i = 1:d
                M(i, i) = M(i, i) * sqrt(2);
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
            ket = Random.unnormalizedPureState(d);
            ket = ket / norm(ket);
        end
        
        function rho = pureNormalizedDensityMatrix(d)
            ket = Random.unnormalizedPureState(d);
            rho = ket*ket';
            rho = rho / trace(rho);
        end
        
    end
    
end
