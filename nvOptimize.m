function [objMax data timings] = nvOptimize(problem, monomials, method, settings, basis, repStructure)
% NVOPTIMIZE Implementation of the Navascues-Vertesi hierarchy
%
%
% problem              - Problem description of class NVProblem, with relevant method implemented
%                        as described below
%
% monomials            - List of monomial basis indices, or {'npa', level}, or {'families' family1 family2 ...}
%                        The list must not include duplicates, even under dimension/commutativity constraints
%
%                        The first level of the hierarchy for 4 operators would be given by
%                        {[] [1] [2] [3] [4]}, where the indices are the indices of variables.
%
%                        For {'npa', level}, constructs all monomials of bounded degree "level", removing
%                        potential duplicates.
%
%                        For {'families' family1 family2 ...}, returns the union of all monomials generated
%                        by all families, removing potential duplicates.
%
%                        A family is given by a list of operator types [t1 t2 ...], where the values of t1, t2, ...
%                        index the ordered partitions returned by problem.operatorTypes.
%
%                        Let Pi = problem.operatorTypes(ti). Then we generate, for each family, all monomials
%                        p1 * p2 * ..., where p1 \in P1, p2 \in P2, etc...
%
%                        To use that method, the problem must override the method "operatorTypes" with a proper
%                        definition; see RAC22.m for an example.
%    
%                        In a bipartite Bell scenario, the two types would be the operators for Alice and Bob
%                        respectively, and the families {'families' [] [1] [2] [1 2]} would correspond to the
%                        definition of the almost quantum set (1 + A + B + AB).
%                        It is a good idea to include the family of the identity "[]" in the list.
%
% method               - The symmetry method used (or lack thereof). Can take one of the values:
%
%                        'none':        Do not symmetrize
%                        'reynolds':    Average the samples without block-diagonalization.
%                        'isotypic':    Average the samples, then block-diagonalize into isotypic components.
%                        'irreps':      Average the samples, then fully block-diagonalize.
%                        'blocks':      Sample in the block-diagonal basis explicitly
%                                       If any of the 'changeOfBasis' or 'representations' parameters is missing,
%                                       the change of basis is computed numerically.
%                        'fastest':     Equivalent to 'blocks'
%
%                        These functions require certain methods to be available on the given problem.
%                        
%                        'reynolds', 'isotypic', 'full' require 'problem.groupDecomposition'
%                        'blocks' needs 'problem.groupDecomposition' if changeOfBasis+representations are not provided
%
% settings             - An instance of NVSettings
% 
% basis                - Optional explicit change of basis matrix provided by the user
%
% repStructure         - Let nR be the number of inequivalent irreducible representations present in the
%                        representation induced by the action of the group on the provided monomial basis.
%
%                        Let permutationMatrix(g) be the (generalized) permutation matrix corresponding to the action
%                        of a group element g.
%
%                        repStructure is either:
%
%                        - a 1xnR vector providing the size of the isotypic components, such that
%                          basis' * permutationMatrix(g) * basis is block-diagonal with blocks of size
%                          given by "repStructure"
%
%                        - a 2xnR integer matrix, where 
%                             - repStructure(1,i) is the dimension of the representation and
%                             - repStructure(2,i) its multiplicity.
%
%                          Then "basis" has to be the change of basis such that basis' * permutationMatrix(g) * basis 
%                          is block-diagonal, with blocks of type kron(eye(d), actionInRep_i(g))
%
%                          Note that the moment matrix is NOT block diagonal in that basis, rather has blocks
%                          corresponding to kron(B, eye(d)), where B is a m x m matrix (m is the multiplicity) and
%                          d is the representation dimension.
    import qdimsum.*
    import qdimsum.group.*
    
    timings = struct;
    data = [];
    
    if settings.checkLevel > 0
        Check.problem(problem, settings);
    end
      
    if isequal(method, 'fastest')
        method = 'blocks';
    end
    
    start = tic;
    if isequal(monomials{1}, 'npa')
        level = monomials{2};
        settings.log(sprintf('Computing monomials for NPA level %d', level));
        monomials = Monomials.npa(problem, level, settings);
    elseif isequal(monomials{1}, 'families')
        families = monomials(2:end);
        familiesStr = cellfun(@(f) ['[' num2str(f) '],'], families, 'UniformOutput', false);
        familiesStr = strcat(familiesStr{:});
        familiesStr = familiesStr(1:end-1);
        settings.log(['Computing monomials for families ' familiesStr]);
        monomials = Monomials.families(problem, families, settings);
    end
    
    nMonomials = length(monomials);
    settings.log(sprintf('Working with %d monomials', nMonomials));
    timings.monomials = toc(start);
    
    UB_NONE = 0;
    UB_ISOTYPIC = 1;
    UB_IRREDUCIBLE = 2;
    userBasis = UB_NONE; % 0: none, 1: isotypic, 2: irreducible
    if nargin == 6
        switch size(repStructure, 1)
          case 1
            userBasis = UB_ISOTYPIC;
          case 2
            userBasis = UB_IRREDUCIBLE;
          otherwise
            error('Wrong format for repStructure');
        end
    end
    
    % compute prerequisities
    switch method
      case 'none'
        needsGroupDecomposition = false;
        needsBasis = UB_NONE;
      case 'reynolds'
        needsGroupDecomposition = true;
        needsBasis = UB_NONE;
      case 'isotypic'
        needsGroupDecomposition = true;
        needsBasis = UB_ISOTYPIC;
      case 'irreps'
        needsGroupDecomposition = true;
        needsBasis = UB_IRREDUCIBLE;
      case 'blocks'
        needsBasis = UB_IRREDUCIBLE;
        needsGroupDecomposition = (userBasis ~= UB_IRREDUCIBLE);
      otherwise
        error('Method should be either: none, reynolds, isotypic, irreps, blocks, fastest')
    end

    start = tic;
    if needsGroupDecomposition || settings.checkLevel > 0
        group = Group(problem.symmetryGroupGenerators);
        f = @(g) Monomials.findMonomialAction(problem, monomials, g, settings);
        monoGroup = group.monomorphism(f);
        gd = Chain.fromGenerators(problem.symmetryGroupGenerators).groupDecomposition;
        monoAction = Monomials.actionDecomposition(problem, gd, monomials, settings);
    end
    timings.groupDecomposition = toc(start);

    start = tic;

    if userBasis < needsBasis
        switch needsBasis
          case UB_ISOTYPIC
            settings.log('Computing isotypic decomposition...');
            [basis reps] = Reps.isotypicComponents(monoGroup, settings);
            repStructure = reps(1,:).*reps(2,:);
          case UB_IRREDUCIBLE
            settings.log('Computing irreducible decomposition...');
            [basis repStructure] = Reps.irreducibleDecomposition(monoGroup, settings);
        end
    end

    if needsBasis > UB_NONE
        if settings.checkLevel > 0
            nRepresentations = size(repStructure, 2);
            if size(repStructure, 1) == 2
                isotypicSizes = repStructure(1,:) .* repStructure(2,:);
            else
                isotypicSizes = repStructure;
            end
            for i = 1:settings.numSamplesChecks
                g = GenPerm.randomFromChain(monoAction);
                G = GenPerm.orthogonalMatrix(g);
                G = basis'*G*basis;
                shift = 0;
                for i = 1:nRepresentations
                    bs = isotypicSizes(i);
                    G(shift+(1:bs), shift+(1:bs)) = 0;
                    shift = shift + bs;
                end
                assert(norm(G) < settings.blockDiagMatTol, 'Error in isotypic components');
            end
            if size(repStructure, 1) == 2
                for i = 1:settings.numSamplesChecks
                    g = GenPerm.randomFromChain(monoAction);
                    G = GenPerm.orthogonalMatrix(g);
                    G = basis'*G*basis;
                    shift = 0;
                    for i = 1:nRepresentations
                        for j = 1:repStructure(2, i)
                            bs = repStructure(1, i);
                            G(shift+(1:bs), shift+(1:bs)) = 0;
                            shift = shift + bs;
                        end
                    end
                    assert(norm(G) < settings.blockDiagMatTol, 'Error in irreducible components');
                end
            end
        end
        switch size(repStructure, 1)
          case 1
            settings.log(['Found isotypic components of size ' num2str(repStructure)]);
            blockSizes = repStructure;
            nRepresentations = length(blockSizes);
            sampleDim = sum(arrayfun(@(x) x*(x+1)/2, blockSizes));
          case 2
            dim = repStructure(1, :);
            mult = repStructure(2, :);
            settings.log('Found representations of dimensions')
            settings.log(num2str(dim))
            settings.log('and multiplicities')
            settings.log(num2str(mult));
            blockSizes = mult;
            nRepresentations = length(blockSizes);
            sampleDim = sum(arrayfun(@(x) x*(x+1)/2, blockSizes));
            % reorder elements so that the moment matrix is fully block-diagonal
            shift = 0;
            for i = 1:nRepresentations
                repSize = dim(i) * mult(i);
                trs = reshape(1:repSize, [dim(i) mult(i)])';
                trs = trs(:);
                range = shift + (1:repSize);
                basis(:, range) = basis(:, range(trs));
                shift = shift + repSize;
            end
        end
    else
        blockSizes = nMonomials;
        sampleDim = nMonomials*(nMonomials+1)/2;
    end
    blockStructure = BlockStructure(blockSizes);
    timings.blockDiagonalization = toc(start); 
    
    timings.sampling = 0;
    timings.rank = 0;
    settings.log('Computing samples');
    samples = zeros(sampleDim, 0);
    objContribs = zeros(1, 0);
    l = 0;
    while 1
        samples = [samples zeros(sampleDim, settings.sampleChunkSize)];
        start = tic;
        for i = 1:settings.sampleChunkSize
            X = problem.sampleOperators;
            K = problem.sampleStateKraus;
            stateDim = size(K, 1);
            krausRank = size(K, 2);
            % compute monomials
            monos = zeros(stateDim, krausRank, nMonomials);
            for j = 1:nMonomials
                mono = K;
                indices = monomials{j};
                for k = length(indices):-1:1
                    mono = X{indices(k)} * mono;
                end
                monos(:,:,j) = mono;
            end
            blocks = {};
            if isequal(method, 'blocks')
                % performs monomial change of basis
                monos = reshape(monos, [stateDim*krausRank nMonomials]) * basis;
                shift = 0;
                for i = 1:nRepresentations
                    block = zeros(mult(i), mult(i));
                    for j = 1:dim(i)
                        range = shift + (1:mult(i));
                        block = block + monos(:,range)'*monos(:,range);
                        shift = shift + mult(i);
                    end
                    block = (block + block')/(2*dim(i));
                    if problem.forceReal
                        block = real(block);
                    end

                    blocks = horzcat(blocks, {block});
                end
            else
                monos = reshape(monos, [stateDim*krausRank nMonomials]);
                chi = monos'*monos; % compute moment matrix
                if problem.forceReal
                    chi = real(chi);
                end
                if ~isequal(method, 'none')
                    chi = GenPerm.symmetrize(chi, monoAction);
                end
                switch method
                  case {'none', 'reynolds'}
                    chi = (chi + chi')/2; % to cater for possible rounding errors
                    blocks = {chi};
                  case 'isotypic'
                    chi = basis'*chi*basis;
                    chi = (chi + chi')/2;
                    shift = 0;
                    for i = 1:nRepresentations
                        range = shift + (1:repStructure(i));
                        blocks = horzcat(blocks, {chi(range, range)});
                        shift = shift + repStructure(i);
                    end
                  case 'irreps'
                    chi = basis'*chi*basis;
                    shift = 0;
                    for i = 1:nRepresentations
                        block = zeros(mult(i), mult(i));
                        for j = 1:dim(i)
                            range = shift + (1:mult(i));
                            block = block + chi(range, range);
                            shift = shift + mult(i);
                        end
                        block = (block + block')/(2*dim(i));
                        blocks = horzcat(blocks, {block});
                    end
                end
            end
            sample = cellfun(@(b) BlockStructure.matToVec(b), blocks, 'UniformOutput', false);
            l = l + 1;
            samples(:,l) = horzcat(sample{:});
            objContribs(l) = problem.computeObjective(X, K);
        end
        timings.sampling = timings.sampling + toc(start);
        start = tic;
        r = rank(samples);
        timings.rank = timings.rank + toc(start);
        if r < l
            break
        end
    end
    settings.log(['Total of ' num2str(r) ' samples']);
    timings.sampling = toc(start);
    
    settings.log('Formulating optimization problem');
    if r == 1
        settings.log('Only one sample, not running semidefinite solver');
        objMax = objContribs(1);
        if settings.checkLevel > 0
            assert(abs(objContribs(2) - objContribs(1)) < settings.checksObjTol, ...
                   ['Failed test of linear dependence of the objective on the moment matrix. ' ...
                    'Increase the degree of monomials in the generating set.']);
        end
        data = 'Did not run optimization';
    else
        start = tic;
        Cons = [];
        x = sdpvar(r-1, 1);
        objMax = objContribs(1);
        objMax = objMax + (objContribs(2:r) - objContribs(1)) * x;
        if settings.checkLevel > 0
            sampleDep = samples(:,1:r) \ samples(:, r+1);
            recObj = dot(objContribs(1:r), sampleDep);
            newObj = objContribs(r+1);
            assert(abs(recObj - newObj) < settings.checksObjTol, ...                   
                   ['Failed test of linear dependence of the objective on the moment matrix. ' ...
                    'Increase the degree of monomials in the generating set.']);
        end
        for i = 1:blockStructure.nBlocks
            vecRange = blockStructure.blockRange(i);
            n = blockStructure.blockSize(i);
            C = BlockStructure.vecToMat(samples(vecRange,1), n);
            A = zeros(n, n, r - 1);
            for j = 1:r-1
                sample = samples(vecRange, j+1) - samples(vecRange, 1);
                A(:,:,j) = BlockStructure.vecToMat(sample, n);
            end
            block = C + reshape(reshape(A, n*n, r - 1) * x, n, n);
            Cons = [Cons
                    block >= 0];
        end
        timings.formulation = toc(start);
        settings.log('Starting optimization');
        data = solvesdp(Cons, -objMax, settings.yalmipSettings);
        objMax = double(objMax);
        timings.yalmip = data.yalmiptime;
        timings.solver = data.solvertime;
    end
    
end
