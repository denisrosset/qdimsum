settings = NVSettings('yalmipSettings', NVSettings.yalmipMOSEK(1e-16));
problem = RAC22;
monomials = {'families' [] [1] [2] [1 2]};
objMax1 = nvOptimize(problem, monomials, 'none', settings);
objMax2 = nvOptimize(problem, monomials, 'reynolds', settings);
objMax3 = nvOptimize(problem, monomials, 'isotypic', settings);
objMax4 = nvOptimize(problem, monomials, 'irreps', settings);
objMax5 = nvOptimize(problem, monomials, 'blocks', settings);
tolerance = 1e-6;
sol = 1/2*(1+1/sqrt(2));
assert(abs(objMax1 - sol) < tolerance, 'Result obtained is outside tolerance');
assert(abs(objMax2 - sol) < tolerance, 'Result obtained is outside tolerance');
assert(abs(objMax3 - sol) < tolerance, 'Result obtained is outside tolerance');
assert(abs(objMax4 - sol) < tolerance, 'Result obtained is outside tolerance');
assert(abs(objMax5 - sol) < tolerance, 'Result obtained is outside tolerance');
abs(objMax1 - sol)
abs(objMax2 - sol)
abs(objMax3 - sol)
abs(objMax4 - sol)
abs(objMax5 - sol)