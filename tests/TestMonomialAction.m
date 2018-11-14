import qdimsum.*                        
settings = NVSettings;
d = 5;
% Test non symmetrized variant
problem = CGLMP(d);
families = {[] [1] [2] [1 2] [2 2]};
relaxation = Relaxation(problem, {'families' families{:}}, settings);
x = relaxation.operatorsGroup.randomElement;
y = relaxation.operatorsGroup.randomElement;
z = GenPerm.compose(x, y);

xm = relaxation.monomialAction(x);
ym = relaxation.monomialAction(y);
zm = relaxation.monomialAction(z);
zm1 = GenPerm.compose(xm, ym);
assert(isequal(zm, zm1));
