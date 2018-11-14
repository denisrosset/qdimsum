classdef TrieNode < handle
% Represents the monomial operator variable indices as a trie, used to reuse
% common subexpressions at the tail to avoid duplicate computations
%
% Leaves (i.e. nChildren = 0) correspond always to a monomial position
    properties
        position;      % monomial position, either 1..nMonomials or 0
                       % this position can be "borrowed", in the sense
                       % that it contains a subexpression that is later
                       % overwritten
        childNodes;    % child nodes of type Trie
        childIndices;  % indices of operator variables for each child
    end
    
    methods
        
        function n = nChildren(self)
            n = length(self.childNodes);
        end
        
        function self = TrieNode(positions, indices)
            n = length(indices);
            % Computes all suffixes, special 0 = no suffix as element is empty
            last = zeros(1, n);
            for i = 1:n
                ind = indices{i};
                if length(ind) == 0
                    last(i) = 0;
                else
                    last(i) = ind(end);
                end
            end
            % Do we have an empty element somewhere?
            thisOne = find(last == 0);
            switch length(thisOne)
              case 0
                % no: the current node will need to borrow its position
                self.position = 0;
              case 1
                % yes: the current node will be present in the final monomial list
                self.position = positions(thisOne);
              otherwise
                % by construction in the case of monomials, duplicates were removed
                error('Duplicate elements are not allowed');
            end
            % finds unique suffixes
            childIndices = unique(last);
            childIndices = childIndices(childIndices > 0);
            nChildren = length(childIndices);
            % compute the child nodes
            childNodes = cell(1, nChildren);
            for i = 1:nChildren
                I = find(last == childIndices(i));
                newIndices = cell(1, length(I));
                for j = 1:length(I)
                    x = indices{I(j)};
                    newIndices{j} = x(1:end-1);
                end
                childNodes{i} = qdimsum.TrieNode(positions(I), newIndices);
            end
            self.childIndices = childIndices;
            self.childNodes = childNodes;
            % if the current node does not appear in the final monomial list, borrow its position
            % as the first monomial in the list is computed last: that position can be reused until
            % the last minute
            if self.position == 0
                self.position = self.borrowPosition;
            end
        end
        
        function p = borrowPosition(self)
        % Finds a children position, always uses the left-most branch as
        % we always have a position on a leaf node
            if self.position == 0
                p = self.childNodes{1}.borrowPosition;
            else
                p = self.position;
            end
        end
        
        function P = program(self, P, operatorVariable, parentPosition)
        % Constructs the program that computes monomials, see "Monomials.program" for a description
        %
        % The first call will be program([], 0, 0)
            P = [P
                 self.position operatorVariable parentPosition];
            for i = self.nChildren:-1:1
                P = self.childNodes{i}.program(P, self.childIndices(i), self.position);
            end
        end
        
    end
end
