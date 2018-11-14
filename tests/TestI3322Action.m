import qdimsum.*                        
settings = NVSettings;
% Test non symmetrized variant
rankCard = [2 2 2 2 2 2];
problem = I3322c(rankCard, 1, 4);
families = {[] [1] [2] [1 2] [1 1] [2 2]};
relaxation = Relaxation(problem, {'families' families{:}}, settings);
x = relaxation.operatorsGroup.randomElement;
y = relaxation.operatorsGroup.randomElement;
z = GenPerm.compose(x, y);

xm = relaxation.monomialAction(x);
ym = relaxation.monomialAction(y);
zm = relaxation.monomialAction(z);
zm1 = GenPerm.compose(xm, ym);
assert(isequal(zm, zm1));

X = GenPerm.orthogonalMatrix(xm);
X1 = GenPerm.slowOrthogonalMatrix(xm);
Y = GenPerm.orthogonalMatrix(ym);
Z = GenPerm.orthogonalMatrix(zm);
assert(isequal(X*Y, Z));
assert(isequal(X, X1));
