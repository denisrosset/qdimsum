% Methods to handle generalized permutations, i.e. permutations
% with an eventual sign flip (i.e. monomial actions with scalars
% equal to +/-1).
classdef GenPerm

    methods (Static)

        % Compute the closure of a set of generators
        function G = closure(generators, n)
            import qdimsum.*
            % TODO: replace by Dimino's algorithm
            id = 1:n;
            R = rand(n, 1);
            G = unique(generators, 'rows');
            nG = size(G, 1);
            G = [G
                 id];
            hash = G * R;
            tol = 1e-12;
            assert(length(find(abs(hash - id * R) < tol)) == 1, 'Identity cannot be part of generators');
            newOne = true;
            while newOne
                size(G, 1)
                newOne = false;
                for i = 1:nG
                    for j = 1:size(G, 1)
                        test = GenPerm.compose(G(i,:), G(j,:));
                        candidates = abs(hash - test*R) < tol;
                        isNew = true;
                        for c = find(candidates)
                            if isequal(test, G(c, :))
                                isNew = false;
                            end
                        end
                        if isNew
                            G = [G
                                 test];
                            hash = [hash
                                    test * R];
                            newOne = true;
                        end
                    end
                end
            end
        end
        
        function g = randomFromChain(genPermChain)
        % Given a group decomposition (cell array), returns a random element
            import qdimsum.*
            C = length(genPermChain);
            dummy = genPermChain{1};
            N = size(dummy, 2);
            g = 1:N;
            for i = 1:C
                c = genPermChain{i};
                ui = c(randi(size(c, 1)), :);
                g = GenPerm.compose(g, ui);
            end
        end

        function matSym = symmetrize(mat, genPermChain)
        % Symmetrizes a matrix under the action of a symmetry group (jointly on the rows and columns)
        % The group is provided as a cell array decomposition
            import qdimsum.*
            C = length(genPermChain);
            N = size(mat, 1);
            matSym = mat;
            for i = C:-1:1
                av = zeros(N, N);
                u = genPermChain{i};
                for j = 1:size(u, 1)
                    c = u(j, :);
                    absC = abs(c);
                    mS = find(c < 0);
                    if isempty(mS)
                        invc = 1:N;
                        invc(c) = 1:N;
                        av = av + matSym(invc, invc);
                    else
                        pS = find(c > 0);
                        pT = absC(pS);
                        mT = absC(mS);
                        av(pT, pT) = av(pT, pT) + matSym(pS, pS);
                        av(mT, mT) = av(mT, mT) + matSym(mS, mS);
                        av(pT, mT) = av(pT, mT) - matSym(pS, mS);
                        av(mT, pT) = av(mT, pT) - matSym(mS, pS);
                    end
                end
                av = av / size(u, 1);
                matSym = av;
            end
        end

        function matSym = slowSymmetrize(mat, genPermChain)
        % Symmetrizes a matrix under the action of a symmetry group (jointly on the rows and columns)
        % The group is provided as a cell array decomposition
            import qdimsum.*
            C = length(genPermChain);
            N = size(mat, 1);
            matSym = mat;
            for i = C:-1:1
                av = zeros(N, N);
                c = genPermChain{i};
                for j = 1:size(c, 1)
                    av = av + GenPerm.matrixImage(c(j,:), matSym);
                end
                av = av / size(c, 1);
                matSym = av;
            end
        end


        function v = vectorImage(gp, v)
        % Image of a column vector under a group element
            v = v .* sign(gp(:));
            v(abs(gp)) = v;
        end
        
        function M = orthogonalMatrix(gp)
        % Permutation/orthogonal matrix corresponding to a group element
            n = length(gp);
            M = sparse(abs(gp), 1:n, sign(gp), n, n);
            %M = sparse([], [], [], n, n);
            %M(abs(gp), :) = D;
        end
        
        function M = slowOrthogonalMatrix(gp)
        % Corresponds to the definition in the paper
            n = length(gp);
            M = zeros(n, n);
            for r = 1:n
                for c = 1:n
                    M(r,c) = (r == abs(gp(c)))*sign(gp(c));
                end
            end
        end
        
        function zeta = matrixImage(gp, chi)
            d = size(chi, 1);
            zeta = chi;
            minusSign = find(gp < 0);
            if length(minusSign) > 0
                zeta(minusSign, :) = -zeta(minusSign, :);
                zeta(:, minusSign) = -zeta(:, minusSign);
            end
            zeta(abs(gp), abs(gp)) = zeta;
        end
                
        function Y = operatorsImage(gp, X)
            import qdimsum.*
            assert(length(gp) == length(X));
            n = length(X);
            Y = cell(1, n);
            for i = 1:n
                j = GenPerm.image(gp, i);
                Y{abs(j)} = sign(j) * X{i};
            end
        end
        
        function j = image(gp, i)
            s = sign(i);
            absI = abs(i);
            j = gp(absI) * s;
        end
        
        function j = findMovedPoint(gp)
        % Returns a point ( > 0 ) moved by the 
            j = 0;
            for i = 1:length(gp)
                if i ~= gp(i)
                    j = i;
                    return
                end
            end
        end
        
        function p = fromCycles(n, varargin)
            import qdimsum.*
            p = 1:n;
            for i = length(varargin):-1:1
                cycle = varargin{i};
                % cycle 2 3 1 means that 2 -> 3, 3 -> 1, 1 -> 2
                cycleImage = [cycle(2:end) cycle(1)];
                newEl = 1:n;
                newEl(cycle) = cycleImage;
                p = GenPerm.compose(newEl, p);
            end
        end
        
        function z = compose(x, y)
            import qdimsum.*
            assert(length(x) == length(y));
            %z = x(y);
            n = length(x);
            z = zeros(1, n);
            for i = 1:n
                z(i) = GenPerm.image(x, GenPerm.image(y, i));
            end
        end
    
        function y = inverse(x)
            n = length(x);
            y = zeros(1, n);
            xAbs = abs(x);
            y(xAbs) = 1:n;
            flip = find(x < 0);
            invFlip = xAbs(flip);
            y(invFlip) = -y(invFlip);
        end
        
        function x = random(n)
            x = randperm(n);
            flip = randi(2, 1, n) == 1;
            x(flip) = -x(flip);
        end
        
        function p1 = toUnsigned(p)
            n = length(p);
            p1 = zeros(1, 2*n);
            for k = 1:n
                k1 = (k-1)*2 + 1;
                i1 = (abs(p(k))-1)*2 + 1;
                if p(k) > 0
                    p1(k1) = i1;
                    p1(k1+1) = i1+1;
                else
                    p1(k1) = i1+1;
                    p1(k1+1) = i1;
                end
            end
        end
        
        function c = directSum(a, b, varargin)
        % Computes the direct sum of permutations
            a = a(:)';
            b = b(:)';
            nA = length(a);
            nB = length(b);
            c = zeros(1, nA + nB);
            c(nA+(1:nB)) = (nA + abs(b)).*sign(b);
            c(1:nA) = a;
            if length(varargin) > 0
                c = qdimsum.GenPerm.directSum(c, varargin{1}, ...
                                              varargin{2:end});
            end
        end
        
    end
    
end
