settings = NVSettings;
problem = RAC22a;
monomials = {'npa' 2};
disp('The objective value is:')
upperBoundSDP = nvOptimize(problem, monomials, 'none', settings)
disp('While the analytical value is:')
sol = 1/2*(1+1/sqrt(2))
