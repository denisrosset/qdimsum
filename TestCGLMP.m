settings = NVSettings('yalmipSettings', sdpsettings('solver', 'mosek', 'verbose', 0)); 

dims = [3 4 5 6 7];
res  = [0.7287 0.7432 0.7569 0.8000 0.8333];
tolerance = 1e-3;

% Test non symmetrized variant
prob = CGLMP(dims(1));
bound = nvOptimize(prob, {'families' [] [1] [2] [1 2] [2 2]}, 'none', settings);
assert(abs(res(1) - bound) < tolerance, 'Result obtained is outside tolerance');

% Test symmetrized variants
for i = 1:length(dims)
    prob = CGLMP(dims(i));
    bound1 = nvOptimize(prob, {'families' [] [1] [2] [1 2] [2 2]}, 'reynolds', settings);
    bound2 = nvOptimize(prob, {'families' [] [1] [2] [1 2] [2 2]}, 'isotypic', settings);
    bound3 = nvOptimize(prob, {'families' [] [1] [2] [1 2] [2 2]}, 'irreps', settings);
    bound4 = nvOptimize(prob, {'families' [] [1] [2] [1 2] [2 2]}, 'blocks', settings);
    assert(abs(res(i) - bound1) < tolerance, 'Result obtained is outside tolerance');
    assert(abs(res(i) - bound2) < tolerance, 'Result obtained is outside tolerance');
    assert(abs(res(i) - bound3) < tolerance, 'Result obtained is outside tolerance');
    assert(abs(res(i) - bound4) < tolerance, 'Result obtained is outside tolerance');
end
