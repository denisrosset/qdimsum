function [vec blockSizes] = blocksToVec(blocks, isHermitian)
    nBlocks = length(blocks);
    d = 0;
    blockSizes = zeros(nBlocks, 2 - isHermitian);
    for i = 1:nBlocks
        if isHermitian
            n = size(blocks{i}, 1);
            d = d + n*(n+1)/2;
            blockSizes(i) = n;
        else
            d = d + prod(size(blocks{i}));
            blockSizes(i, :) = size(blocks{i});
        end
    end
    vec = zeros(d, 1);
    shift = 0;
    for i = 1:nBlocks
        block = blocks{i};
        if isHermitian
            n = size(block, 1);
            for j = 1:n
                nel = n - j + 1;
                vec(shift+(1:nel)) = block(j, j:n);
                shift = shift + nel;
            end
        else
            n = prod(size(blocks{i}));
            vec(shift + (1:n)) = block(:);
            shift = shift + n;
        end
    end
end
