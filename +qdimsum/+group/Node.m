classdef Node < qdimsum.group.Chain
% A node with non-trivial orbit in the BSGS chain
    
    properties
        beta; % base point
        orbit; % orbit of beta, must be sorted, a row vector
        u; % transversal elements, stored as rows in the order given by orbit
        uInv; % inverse of transversal elements
        ownStrongGenerators; % strong generators found at this node
        next; % next node in chain, or [] for the terminal node
    end
        
    methods
        
        function self = Node(beta, n)
            assert(beta > 0);
            assert(n > 0);
            assert(beta <= n);
            self.beta = beta;
            self.n = n;
            self.orbit = beta;
            self.u = 1:n;
            self.uInv = 1:n;
            self.ownStrongGenerators = zeros(0, n);
        end
       
        function l = orbitSize(self)
            l = length(self.orbit);
        end
        
        function r = randomU(self)
            r = self.u(randi(self.orbitSize), :);
        end
                
        function i = orbitIndex(self, b)
        % Returns index of orbit element, or 0 if not found
            i = ismembc2(b, self.orbit);
        end
       
        function s = sift(self, remaining)
            import qdimsum.*
            b = GenPerm.image(remaining, self.beta);
            i = self.orbitIndex(b);
            if i == 0
                s = [];
            else
                self.next.sift(GenPerm.compose(self.uInv(i,:), remaining));
            end            
        end
        
        function addOwnStrongGenerator(self, g)
            self.ownStrongGenerators = [self.ownStrongGenerators; g];
        end
        
        function update(self)
            import qdimsum.*
            S = self.strongGeneratingSet;
            nS = size(S, 1);
            newOrbit = [];
            newU = zeros(0, self.n);
            newUInv = zeros(0, self.n);
            for i = 1:self.orbitSize
                b = self.orbit(i);
                for j = 1:nS
                    g = S(j,:);
                    newB = GenPerm.image(g, b);
                    if self.orbitIndex(newB) == 0 && ~ismember(newB, newOrbit)
                        newOrbit = [newOrbit newB];
                        newEl = GenPerm.compose(g, self.u(i,:));
                        newU = [newU; newEl];
                        newUInv = [newUInv; GenPerm.inverse(newEl)];
                    end
                end
            end
            if length(newOrbit > 0)
                self.orbit = [self.orbit newOrbit];
                self.u = [self.u; newU];
                self.uInv = [self.uInv; newUInv];
                [~, I] = sort(self.orbit);
                self.orbit = self.orbit(I);
                self.u = self.u(I,:);
                self.uInv = self.uInv(I,:);
                self.update
            end
        end
        
        function b = isTerm(self)
            b = false;
        end
        
    end
    
end
