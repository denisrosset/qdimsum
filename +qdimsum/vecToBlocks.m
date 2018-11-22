function blocks = vecToBlocks(vec, blockSizes, isHermitian)
% blocksSizes is n x 1 for Hermitian, and n x 2 for non Hermitian
    nBlocks = size(blockSizes, 1);
    blocks = cell(1, nBlocks);
    shift = 0;
    for i = 1:nBlocks
        if isHermitian
            n = blockSizes(i);
            block = zeros(n, n);
            for j = 1:n
                nel = n - j + 1;
                block(j,j:n) = vec(shift+(1:nel));
                block(j+1:n,j) = conj(vec(shift+(2:nel)));
                shift = shift + nel;
            end
        else
            nR = blockSizes(i, 1);
            nC = blockSizes(i, 2);
            block = reshape(vec(shift + (1:(nR*nC))), [nR nC]);
            shift = shift + nR*nC;
        end
        blocks{i} = block;
    end
end
