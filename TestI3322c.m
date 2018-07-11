settings = NVSettings('sampleChunkSize', 200, 'yalmipSettings', NVSettings.yalmipMOSEK(1e-6), 'verbosityLevel', 2);

rank = [2 2 1 1 2 2];
prob = I3322c(rank, 2, 4);

tic
level = {'families' [] [1] [1 1] [1 1 1] [1 1 1 1] [1 1 1 1 1] [2] ...
         [2 2] [2 2 2] [2 2 2 2] [2 2 2 2 2] [1 2] [1 1 2 2]};
%level = {'npa' 4};
%bound1 = nvOptimize(prob, level, 'reynolds', settings);
time1 = toc;
tic
%bound2 = nvOptimize(prob, level, 'isotypic', settings);
time2 = toc;
tic
%bound3 = nvOptimize(prob, level, 'irreps', settings);
time3 = toc;
tic
bound4 = nvOptimize(prob, level, 'blocks', settings);
time4 = toc;
