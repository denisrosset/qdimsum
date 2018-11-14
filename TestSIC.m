settings = NVSettings('yalmipSettings', sdpsettings('solver', 'mosek', 'verbose', 0), ...
                      'verbosityLevel', 1);
tolerance = 1e-10;
monomials = {'families' [] [1] [2]};

d = 3;
N = 4;
Nys = N*(N-1)/2;
rp = SICRP(N, d, 0);
arc = rp.allRankCards;
nsSamples = [];
sSamples = [];
for i = 1:size(arc, 1)
    rankCard = arc(i, :);
    problem = rp.problemInstance(rankCard);
    relaxation = qdimsum.Relaxation(problem, monomials, settings);
    ns = relaxation.computeBasis;
    s = relaxation.computeSymmetrizedBasis;
    nsSamples = [nsSamples ns];
    sSamples = [sSamples s];
end
%[maxObj orbits objs errored] = SIC.optimize(monomials, 'none', settings)
% randomized ranks 13.8564 no symmetrization
