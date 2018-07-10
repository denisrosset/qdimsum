classdef Monomials
% MONOMIALS - functions generating monomial basis indices
%
% The list of monomials is represented by 'indices', which is a 1 x N cell array containing 
% the N monomial basis elements.
% Each element is a 1 x D vector of operator variable indices, where D is the degree of monomial.

   methods (Static)
       
       function indices = npa(problem, level, settings)
       % NPA - a generator for a monomial basis that enumerates all products of bounded degree
       %
       % The maximal degree enumerated is given by "level", and the function removes duplicates.
           X = problem.sampleOperators;
           nOp = length(X);
           if level == 0
               indices = {[]};
           elseif level == 1
               indices = Monomials.removeDuplicates(problem, horzcat({[]}, num2cell(1:nOp)), settings);
           else
               prev = Monomials.npa(problem, level - 1, settings);
               nPrev = length(prev);
               newIndices = cell(1, nPrev * (nOp + 1));
               ind = 1;
               for i = 1:nPrev
                   newIndices{ind} = prev{i};
                   ind = ind + 1;
                   for j = 1:nOp
                       newIndices{ind} = [prev{i} j];
                       ind = ind + 1;
                   end
               end
               indices = Monomials.removeDuplicates(problem, newIndices, settings);
           end
       end
       
       function indices = families(problem, families, settings)
       % FAMILIES - a generator for a monomial basis from families of monomials
       %
       % A family is a string of operator types, and the operator types are problem dependent
       %
       % For example, in a Bell scenarios, the two types would be the operators for Alice and Bob,
       % and the families would be A, B, AB, AAB, ABB, etc...
       %
       % problem   - Problem considered. Needs to implement the 'operatorTypes' method.
       %
       % families  - Families of operators involved
       %             For example, the families AA AB BB would be written {[1 1] [1 2] [2 2]}
           indices = {};
           for i = 1:length(families)
               newIndices = Monomials.indicesFromFamily(problem, families{i}, settings);
               indices = horzcat(indices, newIndices);
           end
           indices = Monomials.removeDuplicates(problem, indices, settings);
       end
       
   end
   
   methods (Static) % Helpers
       
       % Computes X{indices(1)} * X{indices(2)} ...
       function op = computeFromIndices(X, indices)
           X1 = X{1};
           d = size(X1, 1);
           if length(indices) == 0
               op = eye(d);
           else
               op = X{indices(1)};
               for i = 2:length(indices)
                   op = op * X{indices(i)};
               end
           end
       end
       
       function monos = computeMonomials(X, monomials)
           monos = cellfun(@(I) Monomials.computeFromIndices(X, I), monomials, 'UniformOutput', false);
       end
       
       function monoChain = actionDecomposition(problem, groupDecomposition, monomials, settings)
           C = length(groupDecomposition);
           N = length(monomials);
           monoChain = cell(1, C);
           X = problem.sampleOperators; % verify that the computation goes in the right order
           monoX = Monomials.computeMonomials(X, monomials);
           for i = 1:C
               c = groupDecomposition{i};
               E = size(c, 1);
               mc = zeros(E, N);
               for j = 1:E
                   Y = GenPerm.operatorsImage(c(j, :), X);
                   monoY = Monomials.computeMonomials(Y, monomials);
                   mc(j, :) = Monomials.findOperatorsImage(monoY, monoX, settings);
               end
               monoChain{i} = mc;
           end
       end
       
       function monoGenPerm = findMonomialAction(problem, monomials, opGenPerm, settings)
           X = problem.sampleOperators; % verify that the computation goes in the right order
           monoX = Monomials.computeMonomials(X, monomials);
           Y = GenPerm.operatorsImage(opGenPerm, X);
           monoY = Monomials.computeMonomials(Y, monomials);
           monoGenPerm = Monomials.findOperatorsImage(monoY, monoX, settings);
       end

       function mon = findOperatorsImage(Y, X, settings)
           assert(length(Y) == length(X));
           n = length(X);
           d = size(X{1}, 1);
           H = Random.complexGaussian(d);
           xh = zeros(1, n);
           yh = zeros(1, n);
           for i = 1:n
               xh(i) = sum(sum(H.*X{i}));
               yh(i) = sum(sum(H.*Y{i}));
           end
           mon = zeros(1, n);
           for i = 1:n
               j = find(abs(yh - xh(i)) < settings.monomialTol);
               jInv = find(abs(yh + xh(i)) < settings.monomialTol);
               if length(j) == 1 && length(jInv) == 1
                   assert(j == jInv, 'Cannot find unique image');
               else
                   assert(length(j) + length(jInv) == 1, 'Cannot find unique image');
               end
               if length(j) == 1
                   mon(i) = j;
                   if ~isequal(settings.monomialHist, [])
                       settings.monomialHist.register(yh(j) - xh(i));
                   end
               else
                   mon(i) = -jInv;
                   if ~isequal(settings.monomialHist, [])
                       settings.monomialHist.register(yh(jInv) + xh(i));
                   end
               end
           end
       end
        
        function indices = removeDuplicates(problem, indicesWithDuplicates, settings)
        % TODO use unique rows
           X = problem.sampleOperators;
           X1 = X{1};
           d = size(X1, 1);
           tau = Random.complexGaussian(d);
           N = length(indicesWithDuplicates);
           hash = zeros(1, N);
           for i = 1:N               
               current = indicesWithDuplicates{i};
               h = trace(tau * Monomials.computeFromIndices(X, current));
               if real(h) < 0
                   h = -h;
               end
               hash(i) = h;
           end
           vals = sort(hash);
           tolerance = settings.monomialTol;
           blocks = [1 find(abs(vals(1:end-1) - vals(2:end)) > tolerance)+1 length(vals)+1];
           uvals = [];
           for i = 1:length(blocks)-1
               uvals = [uvals mean(vals(blocks(i):blocks(i+1)-1))];
           end
           uniques = [];
           for i = 1:length(uvals)
               if abs(uvals(i)) > tolerance
                   pick = find(abs(hash - uvals(i)) < tolerance);
                   uniques = [uniques pick(1)];
               end
           end
           uniques = sort(uniques);
           indices = indicesWithDuplicates(uniques);
       end

       function indices = indicesFromFamily(problem, family, settings)
           types = problem.operatorTypes;
           if length(family) == 0
               indices = {[]};
           elseif length(family) == 1
               indices = Monomials.removeDuplicates(problem, num2cell(types{family(1)}), settings);
           else
               prev = Monomials.indicesFromFamily(problem, family(1:end-1), settings);
               nPrev = length(prev);
               type = types{family(end)};
               newIndices = cell(1, nPrev * length(type));
               ind = 1;
               for i = 1:nPrev
                   for j = type
                       newIndices{ind} = [prev{i} j];
                       ind = ind + 1;
                   end
               end
               indices = Monomials.removeDuplicates(problem, newIndices, settings);
           end
       end

   end
   
end
