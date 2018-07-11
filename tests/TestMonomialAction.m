import qdimsum.*                        
settings = NVSettings;
d = 5;
% Test non symmetrized variant
problem = CGLMP(d);
families = {[] [1] [2] [1 2] [2 2]};
monos = Monomials.families(problem, families, settings);
chain = Chain.fromGenerators(problem.symmetryGroupGenerators);
odec = chain.groupDecomposition;
mdec = Monomials.actionDecomposition(problem, odec, monos, settings);
nMonos = length(monos);
x = chain.random;
y = chain.random;
z = GenPerm.compose(x, y);

xm = Monomials.findMonomialAction(problem, monos, x, settings);
ym = Monomials.findMonomialAction(problem, monos, y, settings);
zm = Monomials.findMonomialAction(problem, monos, z, settings);
zm1 = GenPerm.compose(xm, ym);
assert(isequal(zm, zm1));
