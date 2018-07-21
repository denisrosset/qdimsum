settings = NVSettings('yalmipSettings', NVSettings.yalmipMOSEK(0, 'verbose', 1), ...
                      'verbosityLevel', 1);
d = 6;
[monomials U reps] = RAC2d_rep_basis(d);
rac = RAC2d(d);
sol = 1/2*(1+1/sqrt(d));
diff = zeros(2, 10);
for i = 1:10
    [bound5,~,timings5] = nvOptimize(rac, monomials, 'blocks', settings);
    [bound6,~,timings6] = nvOptimize(rac, monomials, 'blocks', settings, U, reps);
    diff(1,i) = abs(bound5-sol);
    diff(2,i) = abs(bound6-sol);
end
