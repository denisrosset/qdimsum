settings = NVSettings;
problem = RAC22b;
monomials = {'npa' 2};
nvOptimize(problem, monomials, 'none', settings)
nvOptimize(problem, monomials, 'reynolds', settings)
nvOptimize(problem, monomials, 'isotypic', settings)
nvOptimize(problem, monomials, 'irreps', settings)
nvOptimize(problem, monomials, 'blocks', settings)
