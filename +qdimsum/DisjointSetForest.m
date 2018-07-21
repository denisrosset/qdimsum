classdef DisjointSetForest < handle
% Implementation of a disjoint-set data structure, using forests.
% 
% See http://en.wikipedia.org/wiki/Disjoint-set_data_structure
% We implement path compression, but not union by rank.

    properties
        parent;
    end
    
    methods
        
        function self = DisjointSetForest(size)
            self.parent = 1:size;
        end
        
        function r = find(self, x)
        % Finds the representative of the set in which x is contained.
            while self.parent(x) ~= x
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
            elseif (xRoot > yRoot)
                self.parent(xRoot) = yRoot;
            end
        end
        
        function p = toPartition(self)
        % Returns the partition represented by this disjoint set forest
        % A partition is a cell(1, nBlocks) where each block is a double(1, nElements)
        % where the elements are sorted
            for x = 1:length(self.parent)
                self.find(x);
            end
            values = unique(self.parent);
            nBlocks = length(values);
            p = cell(1, nBlocks);
            for i = 1:nBlocks
                p{i} = find(self.parent == values(i));
            end
        end
        
    end
    
end
