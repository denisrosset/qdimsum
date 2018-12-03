import qdimsum.*
N = 4;
d = 2;
Nys = N*(N-1)/2;
rankVec = ones(1, 4+12);
S = NVSettings;
P = SICP(N, d, rankVec);
M = Monomials.fromFamilies(P, {[] [1] [2] [1 2] [1 1]}, S);
R = Relaxation(P, M, S);
Z = R.computeRealization('fastest');
O = Z.solveWithYalmip;
assert(abs(O.obj - 10.8989) < 1e-3);
