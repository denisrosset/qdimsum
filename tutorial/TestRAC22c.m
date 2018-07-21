settings = NVSettings;
problem = RAC22c;
monomials = {'families' [] [1] [2] [1 2]};
nvOptimize(problem, monomials, 'blocks', settings)
