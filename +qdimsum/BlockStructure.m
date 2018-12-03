% Class describing a Hermitian/symmetric matrix with a block structure, 
% and enables packing the matrix elements in a vector removing the degeneracies
classdef BlockStructure
    
    properties
        sizes = [];
        colRanges = {}; % nonoverlapping range of indices corresponding to blocks
                        % for example, a matrix blkdiag([1 1 1; 1 1 1; 1 1 1], [1 1; 1 1], [1])
                        % has ranges {1:3 4:5 6}
    end
    
    methods
        
        function self = BlockStructure(sizes)
        % Constructor using the given ranges (see property ranges above)
            n = length(sizes);
            colRanges = cell(1, n);
            shift = 0;
            for i = 1:n
                colRanges{i} = shift + (1:sizes(i));
                shift = shift + sizes(i);
            end
            self.sizes = sizes;
            self.colRanges = colRanges;
        end
        
        function d = dimension(self)
        % Total complex vector space dimension of the block matrix
            d = sum(arrayfun(@(x) x*(x+1)/2, self.blockSizes));
        end
        
        function vec = packBlocks(self, blocks)
        % Packs a blocks in a complex vector
            import qdimsum.*
            dims = self.blockDimensions;
            vec = zeros(1, sum(dims));
            for i = 1:self.nBlocks
                block = blocks{i};
                vec(self.blockRange(i)) = BlockStructure.matToVec(block);
            end
        end
        
        function vec = packMatrix(self, mat)
        % Packs a block diagional matrix in a complex vector
            import qdimsum.*
            dims = self.blockDimensions;
            vec = zeros(1, sum(dims));
            for i = 1:self.nBlocks
                range = self.colRanges{i};
                block = mat(range, range);
                vec(self.blockRange(i)) = BlockStructure.matToVec(block);
            end
        end
        
        function n = blockSize(self, i)
        % Returns the (width/height) size of the i-th block
            range = self.colRanges{i};
            n = length(range);
        end
        
        function d = blockDimension(self, i)
        % Returns the complex vector space dimension of the i-th block
            n = self.blockSize(i);
            d = n*(n+1)/2;
        end
        
        function n = nBlocks(self)
        % Returns the number of blocks
            n = length(self.sizes);
        end
        
        function sizes = blockSizes(self)
        % Returns the block sizes
            sizes = arrayfun(@(i) self.blockSize(i), 1:self.nBlocks);
        end
        
        function dims = blockDimensions(self)
        % Returns the complex vector dimension of all blocks
            dims = arrayfun(@(i) self.blockDimension(i), 1:self.nBlocks);
        end
        
        function vecRange = blockRange(self, i)
        % Returns the range of indices in the packed vector corresponding to the i-th block
            cumdims = [0 cumsum(self.blockDimensions)];
            vecRange = cumdims(i)+(1:self.blockDimension(i));
        end
        
        function mat = extractBlock(self, vec, i)
        % Extracts the i-th block from the packed vector
            import qdimsum.*
            data = vec(self.blockRange(i));
            mat = BlockStructure.vecToMat(data, self.blockSize(i));
        end
        
        function mat = unpackMatrix(self, vec)
        % Extracts a block diagonal matrix from a packed vector
            blocks = self.unpackBlocks(vec);
            mat = blkdiag(blocks{:});
        end
        
        function blocks = unpackBlocks(self, vec)
        % Extracts the blocks of a block diagonal matrix from a packed vector
            blocks = arrayfun(@(i) self.extractBlock(vec, i), 1:self.nBlocks, 'UniformOutput', false);
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
