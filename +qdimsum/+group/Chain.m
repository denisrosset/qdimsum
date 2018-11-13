classdef Chain < handle
% Describes a (generalized) permutation group using a BSGS chain
% constructed using randomized algorithms
    properties
        n; % domain size (not including negative elements)
    end
    
    methods (Abstract)
        b = isTerm(self) % whether it is a terminal node
    end

    methods
        
        % Note: most of those abstract methods are tail-recursive
        % and could be rewritten as loops
        
        function l = length(self)
            l = 0;
            it = self;
            while ~it.isTerm
                l = l + 1;
                it = it.next;
            end
        end
        
        function r = random(self)
            import qdimsum.GenPerm;
            it = self;
            r = 1:self.n;
            while ~it.isTerm
                r = GenPerm.compose(r, it.randomU);
                it = it.next;
            end
        end

        function o = order(self)
            it = self;
            o = java.math.BigInteger(1);
            while ~it.isTerm
                o = o.multiply(java.math.BigInteger(it.orbitSize));
                it = it.next;
            end
        end
        
        function remaining = sift(self, el)
            import qdimsum.GenPerm;
            it = self;
            remaining = el;
            while ~it.isTerm
                b = GenPerm.image(remaining, it.beta);
                i = it.orbitIndex(b);
                if i == 0
                    s = [];
                    return
                else
                    remaining = GenPerm.compose(it.uInv(i, :), remaining);
                end
                it = it.next;
            end
        end
        
        function s = strongGeneratingSet(self)
            s = zeros(0, self.n);
            it = self;
            while ~it.isTerm
                s = vertcat(s, it.ownStrongGenerators);
                it = it.next;
            end
        end
                        
        function gd = groupDecomposition(self)
            gd = {};
            it = self;
            while ~it.isTerm
                gd = horzcat(gd, {it.u});
                it = it.next;
            end
        end
        
    end
    
    methods (Static)
      
        function res = siftAndUpdateBaseFrom(at, g)
        % Returns an element obtained by sifting g through the chain starting at "at"
        % inserting new base points as required, returns either [],
        % when g can be completely sifted, or
        %
        % {remaining node} where "remaining ~= identity" is the element
        % resulting from an incomplete sift, and node is where
        % this element should be inserted as a strong generator
        %
        % Based on Holt (2005) RANDOMSCHREIER procedure, page 98.
            n = length(g);
            node = at.next; % we have the chain at -> node
            if node.isTerm
                beta = qdimsum.GenPerm.findMovedPoint(g);
                if beta > 0
                    newNode = qdimsum.group.Node(beta, n);
                    % insert node so that: at -> newNode -> node
                    at.next = newNode;
                    newNode.next = node;
                    res = {g newNode};
                else
                    res = [];
                end
            else
                b = qdimsum.GenPerm.image(g, node.beta);
                i = node.orbitIndex(b);
                if i == 0
                    res = {g node};
                else
                    h = qdimsum.GenPerm.compose(node.uInv(i, :), g);
                    % recursive call
                    res = qdimsum.group.Chain.siftAndUpdateBaseFrom(at.next, h);
                end
            end
        end
        
        function addStrongGenerator(start, node, g)
        % Adds a strong generator to a node, and updates the orbits
        % of all its parents
            node.addOwnStrongGenerator(g);
            c = start.next;
            while ~isequal(c, node.next)
                c.update;
                c = c.next;
            end
        end
        
        function b = siftAndAddStrongGenerator(start, g)
        % Sifts the given element through the chain and returns whether
        % a new strong generator has been discovered
            res = qdimsum.group.Chain.siftAndUpdateBaseFrom(start, g);
            if isequal(res, [])
                b = false;
            else
                g = res{1};
                node = res{2};
                qdimsum.group.Chain.addStrongGenerator(start, node, g);
                b = true;
            end
        end
        
        function chain = randomConstruction(randomElement, order, numTests)
        % Constructs a BSGS chain from an oracle that returns random group elements.
        % 
        % If the order of the group is provided, the algorithm is randomized but
        % cannot fail. If the order of the group is not provided (or = -1), then
        % a number of tests numTests is performed so that the probability of failure
        % is 2^-numTests (provided that randomElement returns group elements sampled
        % uniformly at random).
            if nargin < 3
                numTests = 128;
            end
            if nargin < 2 || isequal(order, [])
                hasOrder = false;
            else
                hasOrder = true;
            end
            n = length(randomElement);
            start = qdimsum.group.Start.emptyChain(n);
            numSifted = 0;
            while (~hasOrder && numSifted <= numTests) || (hasOrder && start.next.order.compareTo(order) < 0)
                s = randomElement();
                b = qdimsum.group.Chain.siftAndAddStrongGenerator(start, s);
                if b
                    numSifted = 0;
                else
                    numSifted = numSifted + 1;
                end
            end
            chain = start.next;
        end
        
        function chain = fromGenerators(generators, order, numTests)
        % Constructs a BSGS chain from a list of generators, given a row vectors
        % in a matrix.
        %
        % If the order of the group is provided, the algorithm is randomized but
        % cannot fail.
        %
        % Random elements are produced using the product replacement algorithm;
        % that algorithm can be unsatisfactory when the group is a direct product
        % of a large number of copies of the same finite simple group 
        % (see Holt et al. Handbook of Computational Group Theory (2005), p. 69).
            if nargin < 3
                numTests = 128;
            end
            if nargin < 2 || isequal(order, [])
                order = [];
            end
            bag = qdimsum.group.RandomBag(generators);
            chain = qdimsum.group.Chain.randomConstruction(@() bag.sample, order, numTests);
        end
        
    end
    
end
