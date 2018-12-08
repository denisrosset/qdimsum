import qdimsum.*
clear all
N = 4;
d = 2;
mu = 1;
P = SIC1(N, d, mu);
NS = NVSettings('yalmipSettings', NVSettings.yalmipMOSEK(0));
F1 = Monomials.families(P, {[] [1] [2] [3]}, NS);
M = Monomials(P, F1, NS, 'Custom');
R = Relaxation(P, M, NS);

Z1 = R.computeRealization('none');
Z2 = R.computeRealization('reynolds');
Z3 = R.computeRealization('fastest');

[I1 time1] = Z1.solveWithYalmip;
[I2 time2] = Z2.solveWithYalmip;
[I3 time3] = Z3.solveWithYalmip;
