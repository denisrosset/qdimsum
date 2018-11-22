import qdimsum.*                        
settings = NVSettings;
d = 5;
% Test non symmetrized variant
problem = CGLMP(d);
families = {[] [1] [2] [1 2] [2 2]};
monomials = Monomials.fromFamilies(problem, families);
relaxation = Relaxation(problem, monomials, settings);
x = relaxation.operatorsGroup.randomElement;
y = relaxation.operatorsGroup.randomElement;
z = GenPerm.compose(x, y);

xm = relaxation.monomials.action(x);
ym = relaxation.monomials.action(y);
zm = relaxation.monomials.action(z);
zm1 = GenPerm.compose(xm, ym);
assert(isequal(zm, zm1));
