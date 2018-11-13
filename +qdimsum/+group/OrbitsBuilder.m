classdef OrbitsBuilder < handle
% Implementation of a disjoint-set data structure, using forests.
% 
% See http://en.wikipedia.org/wiki/Disjoint-set_data_structure
% We implement path compression, but not union by rank.

    properties
        n;
        parent;
        size;
        index;
        nOrbits;
        orbits;
    end
    
    methods (Static)
        function B = fromPerm(generators)
            n = size(generators, 2);
            B = qdimsum.group.OrbitsBuilder(n);
            for i = 1:size(generators, 1)
                for x = 1:n
                    xi = generators(i, x);
                    B.union(x, xi);
                end
            end
            B.computeOrbits;
        end
    end
    
    methods
        
        function self = OrbitsBuilder(n)
            self.n = n;
            self.parent = 1:n;
            self.nOrbits = n;
            self.size = ones(1, n);
            self.index = zeros(1, n);
        end
        
        function r = find(self, x)
        % Finds the representative of the set in which x is contained.
            while self.parent(x) ~= x
                % uses path splitting
                nxt = self.parent(x);
                self.parent(x) = self.parent(nxt);
                x = nxt;
            end
            r = x;
        end
        
        function union(self, x, y)
        % Combines the trees containing x and y. The representative of the union is chosen 
        % as the minimal representative of x and y.
            xRoot = self.find(x);
            yRoot = self.find(y);
            if (xRoot < yRoot)
                self.parent(yRoot) = xRoot;
                self.nOrbits = self.nOrbits - 1;
                self.size(xRoot) = self.size(xRoot) + self.size(yRoot);
                self.size(yRoot) = 0;
            elseif (xRoot > yRoot)
                self.parent(xRoot) = yRoot;
                self.nOrbits = self.nOrbits - 1;
                self.size(yRoot) = self.size(yRoot) + self.size(xRoot);
                self.size(xRoot) = 0;
            end
        end
        
        function computeOrbits(self)
        % Returns the partition represented by this disjoint set forest
        % A partition is a cell(1, nBlocks) where each block is a double(1, nElements)
        % where the elements are sorted
            orbit = 1;
            orbitIndex = ones(1, self.nOrbits);
            self.orbits = cell(1, self.nOrbits);
            for x = 1:self.n
                x0 = self.find(x);
                if self.index(x0) == 0
                    % new orbit
                    self.index(x0) = orbit;
                    self.orbits{orbit} = zeros(self.size(x), 1);
                    orbit = orbit + 1;
                end
                o = self.index(x0);
                self.index(x) = o;
                orb = self.orbits{o};
                orb(orbitIndex(o)) = x;
                self.orbits{o} = orb;
                orbitIndex(o) = orbitIndex(o) + 1;
            end
        end
        
    end
    
end
