classdef Chain < handle

    properties
        n; % domain size (not including negative elements)
    end
    
    methods (Abstract)
        
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
        
        function chain = schreierSims(generators, order, numTests)
            if nargin < 3
                numTests = 128;
            end
            if nargin < 2 || isequal(order, [])
                order = -1;
            end
            n = size(generators, 2);
            start = Start.emptyChain(n);
            bag = RandomBag(generators);
            numSifted = 0;
            while (order == -1 && numSifted <= numTests) || (order >= 1 && start.next.order < order)
                s = bag.sample;
                b = Chain.siftAndAddStrongGenerator(start, s);
                if b
                    numSifted = 0;
                else
                    numSifted = numSifted + 1;
                end
            end
            chain = start.next;
        end
        
    end
    
end
