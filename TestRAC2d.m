settings = NVSettings('yalmipSettings', sdpsettings('solver', 'mosek', 'verbose', 0), ...
                      'verbosityLevel', 0);
tolerance = 1e-6;
monomials = {'families' [] [1] [2] [1 2]};
for d = 2:5
    rac = RAC2d(d);
    sol = 1/2*(1+1/sqrt(d));
    bound1 = nvOptimize(rac, monomials, 'reynolds', settings);
    bound2 = nvOptimize(rac, monomials, 'isotypic', settings);
    bound3 = nvOptimize(rac, monomials, 'irreps', settings);
    bound4 = nvOptimize(rac, monomials, 'blocks', settings);
    assert(abs(bound1 - sol) < tolerance, 'Result obtained is outside tolerance');
    assert(abs(bound2 - sol) < tolerance, 'Result obtained is outside tolerance');
    assert(abs(bound3 - sol) < tolerance, 'Result obtained is outside tolerance');
    assert(abs(bound4 - sol) < tolerance, 'Result obtained is outside tolerance');
end
