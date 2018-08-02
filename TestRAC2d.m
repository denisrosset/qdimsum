settings = NVSettings;
tolerance = 1e-6;
monomials = {'families' [] [1] [2] [1 2]};
for d = 3:6
    [monomialsF U reps] = RAC2d_rep_basis(d);
    rac = RAC2d(d);
    sol = 1/2*(1+1/sqrt(d));
    bound1 = nvOptimize(rac, monomials, 'reynolds', settings);
    bound2 = nvOptimize(rac, monomials, 'isotypic', settings);
    bound3 = nvOptimize(rac, monomials, 'irreps', settings);
    bound4 = nvOptimize(rac, monomials, 'blocks', settings);
    bound5 = nvOptimize(rac, monomialsF, 'irreps', settings, U, reps);
    bound6 = nvOptimize(rac, monomialsF, 'blocks', settings, U, reps);
    assert(abs(bound1 - sol) < tolerance, 'Result obtained is outside tolerance');
    assert(abs(bound2 - sol) < tolerance, 'Result obtained is outside tolerance');
    assert(abs(bound3 - sol) < tolerance, 'Result obtained is outside tolerance');
    assert(abs(bound4 - sol) < tolerance, 'Result obtained is outside tolerance');
    assert(abs(bound5 - sol) < tolerance, 'Result obtained is outside tolerance');
    assert(abs(bound6 - sol) < tolerance, 'Result obtained is outside tolerance');
end
