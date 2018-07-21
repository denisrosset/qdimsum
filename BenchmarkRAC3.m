settings = NVSettings('yalmipSettings', NVSettings.yalmipMOSEK(0, 'verbose', 0), 'verbosityLevel', 1);
d = 3;
[monomials U reps] = RAC2d_rep_basis(d);
rac = RAC2d(d);
sol = 1/2*(1+1/sqrt(d));
N = 20;
objs = zeros(N, 7);
timings = cell(N, 7);
for i = 1:N
    yalmip clear;
    [bound1,~,timings1] = nvOptimize(rac, monomials, 'none', settings);
    yalmip clear;
    [bound2,~,timings2] = nvOptimize(rac, monomials, 'reynolds', settings);
    yalmip clear;
    [bound3,~,timings3] = nvOptimize(rac, monomials, 'isotypic', settings);
    yalmip clear;
    [bound4,~,timings4] = nvOptimize(rac, monomials, 'irreps', settings);
    yalmip clear;
    [bound5,~,timings5] = nvOptimize(rac, monomials, 'blocks', settings);
    yalmip clear;
    [bound6,~,timings6] = nvOptimize(rac, monomials, 'irreps', settings, U, reps);
    yalmip clear;
    [bound7,~,timings7] = nvOptimize(rac, monomials, 'blocks', settings, U, reps);
    objs(i, 1) = bound1;
    objs(i, 2) = bound2;
    objs(i, 3) = bound3;
    objs(i, 4) = bound4;
    objs(i, 5) = bound5;
    objs(i, 6) = bound6;
    objs(i, 7) = bound7;
    timings{i, 1} = timings1;
    timings{i, 2} = timings2;
    timings{i, 3} = timings3;
    timings{i, 4} = timings4;
    timings{i, 5} = timings5;
    timings{i, 6} = timings6;
    timings{i, 7} = timings7;
end
save BenchmarkRAC3.mat objs timings