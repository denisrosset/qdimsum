classdef Chain < handle
% Describes a (generalized) permutation group using a BSGS chain
% constructed using randomized algorithms
    properties
        n; % domain size (not including negative elements)
    end
    
    methods (Abstract)
        
        % Note: most of those abstract methods are tail-recursive
        % and could be rewritten as loops
        
        l = length_(self, acc)
        
        g = random_(self, acc)
        
        o = order_(self, acc)
        
        s = sift(self, remaining)
        
        s = strongGeneratingSet_(self, acc)
        
        b = isTerm(self) % whether it is a terminal node
        
        gd = groupDecomposition_(self, acc)
        
    end
    
    methods
       
        function l = length(self)
            l = self.length_(0);
        end
        
        function r = random(self)
            r = self.random_(1:self.n);
        end
        
        function o = order(self)
            o = self.order_(1);
        end
        
        function s = strongGeneratingSet(self)
            s = self.strongGeneratingSet_(zeros(0, self.n));
        end
        
        function gd = groupDecomposition(self)
            gd = self.groupDecomposition_({});
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
            import qdimsum.*
            n = length(g);
            node = at.next; % we have the chain at -> node
            if node.isTerm
                beta = GenPerm.findMovedPoint(g);
                if beta > 0
                    newNode = Node(beta, n);
                    % insert node so that: at -> newNode -> node
                    at.next = newNode;
                    newNode.next = node;
                    res = {g newNode};
                else
                    res = [];
                end
            else
                b = GenPerm.image(g, node.beta);
                i = node.orbitIndex(b);
                if i == 0
                    res = {g node};
                else
                    h = GenPerm.compose(node.uInv(i, :), g);
                    % recursive call
                    res = Chain.siftAndUpdateBaseFrom(at.next, h);
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
            import qdimsum.*
            res = Chain.siftAndUpdateBaseFrom(start, g);
            if isequal(res, [])
                b = false;
            else
                g = res{1};
                node = res{2};
                Chain.addStrongGenerator(start, node, g);
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
            import qdimsum.*
            if nargin < 3
                numTests = 128;
            end
            if nargin < 2 || isequal(order, [])
                order = -1;
            end
            n = length(randomElement);
            start = Start.emptyChain(n);
            numSifted = 0;
            while (order == -1 && numSifted <= numTests) || (order >= 1 && start.next.order < order)
                s = randomElement();
                b = Chain.siftAndAddStrongGenerator(start, s);
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
            import qdimsum.*
            if nargin < 3
                numTests = 128;
            end
            if nargin < 2 || isequal(order, [])
                order = -1;
            end
            bag = RandomBag(generators);
            chain = Chain.randomConstruction(@() bag.sample, order, numTests);
        end
        
    end
    
end
