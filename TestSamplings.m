import qdimsum.*
S = NVSettings('blockDiagRefine', false, 'checkLevel', 0);
P = RAC2d(7); R = Relaxation(P, Monomials.fromFamilies(P, {[] [1] [2] [1 2]}), S);
X = P.sampleOperators; K = P.sampleStateKraus;
A = R.blockDiagMomentMatrix('blocks', X, K);
B = R.blockDiagMomentMatrix('irreps', X, K);
C = R.blockDiagMomentMatrix('isotypic', X, K);
D = R.symmetrizedMomentMatrix(X, K);
