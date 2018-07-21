settings = NVSettings('yalmipSettings', NVSettings.yalmipMOSEK(0, 'verbose', 0), 'verbosityLevel', 1);
problem = I3322c([1 1 1 1 1 1], 1, 2);
sol = 5;
N = 20;
objs = zeros(N, 4);
timings = cell(N, 4);
monomials = {'npa' 4};
for i = 1:N
    yalmip clear;
    [bound2,~,timings2] = nvOptimize(problem, monomials, 'reynolds', settings);
    yalmip clear;
    [bound3,~,timings3] = nvOptimize(problem, monomials, 'isotypic', settings);
    yalmip clear;
    [bound4,~,timings4] = nvOptimize(problem, monomials, 'irreps', settings);
    yalmip clear;
    [bound5,~,timings5] = nvOptimize(problem, monomials, 'blocks', settings);
    objs(i, 1) = bound2;
    objs(i, 2) = bound3;
    objs(i, 3) = bound4;
    objs(i, 4) = bound5;
    timings{i, 1} = timings2;
    timings{i, 2} = timings3;
    timings{i, 3} = timings4;
    timings{i, 4} = timings5;
    save BenchmarkI3322big.mat objs timings
end
