% Class describing a Hermitian/symmetric matrix with a block structure, 
% and enables packing the matrix elements in a vector removing the degeneracies
classdef BlockStructure
    
    properties
        ranges = {}; % nonoverlapping range of indices corresponding to blocks
                     % for example, a matrix blkdiag([1 1 1; 1 1 1; 1 1 1], [1 1; 1 1], [1])
                     % has ranges {1:3 4:5 6}
    end
    
    methods
        
        % Constructor using the given ranges (see property ranges above)
        function self = BlockStructure(ranges)
            self.ranges = ranges;
        end
        
        % Total complex vector space dimension of the block matrix
        function d = dimension(self)
            d = sum(arrayfun(@(x) x*(x+1)/2, self.blockSizes));
        end
        
        % Packs a block diagional matrix in a complex vector
        function vec = pack(self, mat)
            dims = self.blockDimensions;
            vec = zeros(1, sum(dims));
            for i = 1:self.numBlocks
                range = self.ranges{i};
                block = mat(range, range);
                vec(self.blockRange(i)) = BlockStructure.matToVec(block);
            end
        end
        
        % Returns the (width/height) size of the i-th block
        function n = blockSize(self, i)
            range = self.ranges{i};
            n = length(range);
        end
        
        % Returns the complex vector space dimension of the i-th block
        function d = blockDimension(self, i)
            n = self.blockSize(i);
            d = n*(n+1)/2;
        end
        
        % Returns the number of blocks
        function n = numBlocks(self)
            n = length(self.ranges);
        end
        
        % Returns the block sizes
        function sizes = blockSizes(self)
            sizes = arrayfun(@(i) self.blockSize(i), 1:self.numBlocks);
        end
        
        % Returns the complex vector dimension of all blocks
        function dims = blockDimensions(self)
            dims = arrayfun(@(i) self.blockDimension(i), 1:self.numBlocks);
        end
        
        % Returns the range of indices in the packed vector corresponding to the i-th block
        function vecRange = blockRange(self, i)
            cumdims = [0 cumsum(self.blockDimensions)];
            vecRange = cumdims(i)+(1:self.blockDimension(i));
        end
        
        % Extracts the i-th block from the packed vector
        function mat = extractBlock(self, vec, i)
            data = vec(self.blockRange(i));
            mat = BlockStructure.vecToMat(data, self.blockSize(i));
        end
        
        % Extracts a block diagonal matrix from a packed vector
        function mat = unpack(self, vec)
            blocks = arrayfun(@(i) self.extractBlock(vec, i), 1:self.numBlocks, 'UniformOutput', false);
            mat = blkdiag(blocks{:});
        end
    end
    
    methods (Static)
        % Unpacks a lower triangle vector into a Hermitian matrix
        function mat = vecToMat(vec, n)
            mat = zeros(n, n);
            shift = 0;
            for i = 1:n
                nel = n - i + 1;
                mat(i,i:n) = vec(shift+(1:nel));
                mat(i+1:n, i) = conj(vec(shift+(2:nel)));
                shift = shift + nel;
            end
        end
        % Packs the lower triangle in a vector, using col-major ordering
        function vec = matToVec(mat)
            n = size(mat, 1);
            assert(size(mat, 2) == n);
            vec = zeros(1, n*(n+1)/2);
            shift = 0;
            for i = 1:n
                nel = n - i + 1;
                vec(shift+(1:nel)) = mat(i,i:n);
                shift = shift + nel;
            end
        end
    end
    
end
