settings = NVSettings('sampleChunkSize', 200);

rank = [2 2 1 1 2 2];
prob = I3322c(rank, 2, 4);

level = {'npa' 4};
bound1 = nvOptimize(prob, level, 'reynolds', settings);
bound2 = nvOptimize(prob, level, 'isotypic', settings);
bound3 = nvOptimize(prob, level, 'irreps', settings);
bound4 = nvOptimize(prob, level, 'blocks', settings);
