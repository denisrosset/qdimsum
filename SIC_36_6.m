import qdimsum.*
clear all
N = 36;
d = 6;
mu = 1;
P = SIC1(N, d, mu);
NS = NVSettings('yalmipSettings', NVSettings.yalmipMOSEK(0));
F1 = Monomials.families(P, {[] [1] [2] [3] [1 1]}, NS);
F2 = P.monomialsInObjective;
M = Monomials(P, horzcat(F1, F2), NS, 'Custom');
R = Relaxation(P, M, NS);

Z3 = R.computeRealization('fastest');
[I3 time3] = Z3.solveWithYalmip;

save res_36_6.mat
