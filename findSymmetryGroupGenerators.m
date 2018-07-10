function generators = findSymmetryGroupGenerators(problem, settings)
% findSymmetryGroupGenerators Find the generators of a problem from its ambient group
%
% Inputs
%
%    problem: A NV optimization problem with a defined set of generators for its ambient group
%   settings: An instance of NVSettings (the parameter findSymTol is used to check for invariance)
%
% Output
%
% generators: A matrix of generators for the symmetry group of the problem, one generator per row
    import qdimsum.*
    X = problem.sampleOperators;
    n = length(X);
    K = problem.sampleStateKraus;
    obj = problem.computeObjective(X, K);
    symStart = Start.emptyChain(n);
    ambGens = problem.ambientGroupGenerators;
    if length(findrows(ambGens, 1:n)) == size(ambGens, 1)
        settings.log('Warning: ambient group of the given problem is trivial');
        generators = 1:n;
        return
    end
    ambChain = Chain.schreierSims(ambGens);
    
    walkChain(ambChain, 1:n);
    generators = symStart.next.strongGeneratingSet;
    function walkChain(chain, g)
        if chain.isTerm
            X1 = GenPerm.operatorsImage(g, X);
            obj1 = problem.computeObjective(X1, K);
            diff = abs(obj1 - obj);
            if diff < settings.findSymTol
                if ~isequal(settings.findSymHist, [])
                    settings.findSymHist.register(diff);
                end
                Chain.siftAndAddStrongGenerator(symStart, g);
            end
        else
            for i = 1:chain.orbitSize
                walkChain(chain.next, GenPerm.compose(g, chain.u(i,:)));
            end
        end
    end
end
