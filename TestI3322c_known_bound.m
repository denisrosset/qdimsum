settings = NVSettings('sampleChunkSize', 200, 'yalmipSettings', NVSettings.yalmipMOSEK(0), 'verbosityLevel', 2);

rank = [1 1 1 1 1 1];
prob = I3322c(rank, 1, 2);
tolerance = 1e-6;

level = {'families' [] [1] [2] [1 2] [1 1] [2 2] [1 1 1] [2 2 2]};
bound1 = nvOptimize(prob, level, 'none', settings);
bound2 = nvOptimize(prob, level, 'reynolds', settings);
bound3 = nvOptimize(prob, level, 'isotypic', settings);
bound4 = nvOptimize(prob, level, 'irreps', settings);
bound5 = nvOptimize(prob, level, 'blocks', settings);
sol = 5;
assert(abs(bound1 - sol) < tolerance, 'Result obtained is outside tolerance');
assert(abs(bound2 - sol) < tolerance, 'Result obtained is outside tolerance');
assert(abs(bound3 - sol) < tolerance, 'Result obtained is outside tolerance');
assert(abs(bound4 - sol) < tolerance, 'Result obtained is outside tolerance');
assert(abs(bound5 - sol) < tolerance, 'Result obtained is outside tolerance');
