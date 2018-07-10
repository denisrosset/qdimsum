% Methods to handle generalized permutations, i.e. permutations
% with an eventual sign flip (i.e. monomial actions with scalars
% equal to +/-1).
classdef GenPerm

    methods (Static)

        % Compute the closure of a set of generators
        function G = closure(generators, n)
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
        
        function matSym = symmetrize(mat, genPermChain)
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
        
        function x = randomInDecomposition(chain)
            n = size(chain{1}, 2);
            x = 1:n;
            for i = 1:length(chain)
                c = chain{i};
                s = randi(size(c, 1));
                x = GenPerm.compose(x, c(s, :));
            end
        end
    end
    
end
