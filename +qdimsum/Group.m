% Group of generalized permutations
classdef Group < handle
    
    properties (GetAccess = public, SetAccess = protected)
        n;
        generators;
    end
    
    properties (Access = ?qdimsum.Group)
        chain_ = [];
        order_ = []; % java.math.BigInteger
        decomposition_ = [];
        permOrbits_ = [];
        permOrbitRestrictions_ = [];
        phaseConfiguration_ = [];
        randomBag_ = [];
    end
        
    methods
        
        function G1 = monomorphism(self, f)
        % Given a monomorphism f: GenPerm -> GenPerm, returns the
        % image f(G). The monomorphism must be able to operate on
        % elements in parallel, given as matrix rows.
            generators1 = f(self.generators);
            G1 = qdimsum.Group(generators1);
            if ~isequal(self.decomposition_, [])
                G1.decomposition_ = cellfun(f, self.decomposition_, 'UniformOutput', false);
            end
            if ~isequal(self.order_, [])
                G1.order_ = self.order_;
            end
        end
                
        function m = nGenerators(self)
        % Returns the number of elements generating this group
            m = size(self.generators, 1);
        end
        
        function self = Group(generators, order)
        % Constructs a generalized permutation group from a list of generators
        %
        % The generators are provided, one per row, in a matrix.
        %
        % For example, the permutation group S_4 is generated by
        % Group([2 3 4 1; 2 1 3 4])
        %
        % Optionally, the order of the group can be provided if known.
            self.n = size(generators, 2);
            assert(self.n >= 1, 'The column dimension of generators defines the domain size, cannot be 0');
            self.generators = generators;
            if size(generators, 1) == 0
                self.order_ = java.math.BigInteger.ONE;
                self.chain_ = qdimsum.group.Term(self.n);
            end
            if nargin > 1
                if ~isa(order, 'java.math.BigInteger')
                    order = java.math.BigInteger.valueOf(order);
                end
                self.order_ = order;
            end
        end
        
        function d = decomposition(self)
        % Computes a decomposition of the group
        % as a cell array of size m x 1, corresponding
        % to sets U_1, ..., U_m, so that every group element
        % is uniquely written as u_1 * u_2 * ... u_m and
        % u_i \in U_i
        %
        % The sets U_i are provided as matrices, where
        % each row is a group element
            if isequal(self.decomposition_, [])
                self.decomposition_ = self.chain.groupDecomposition;
            end
            d = self.decomposition_;
        end
        
        function c = chain(self)
        % Returns the group BSGS chain
            if isequal(self.chain_, [])
                if isequal(self.order_, [])
                    self.chain_ = qdimsum.group.Chain.fromGenerators(self.generators);
                else
                    self.chain_ = qdimsum.group.Chain.fromGenerators(self.generators, self.order_);
                end
            end
            c = self.chain_;
        end
        
        function r = randomElement(self)
        % Returns a random group element with "pretty good" distribution;
        % faster than "uniformRandomElement" but does not provide the same
        % guarantees
            if isequal(self.randomBag_, [])
                self.randomBag_ = qdimsum.group.RandomBag(self.generators);
            end
            r = self.randomBag_.sample;
        end
        
        function r = uniformRandomElement(self)
        % Returns a random group element guaranteed to be chosen uniformly at random
            r = self.chain.random;
        end
        
        function o = order(self)
        % Returns the order of the group as a java.math.BigInteger
        %
        % use G.order.doubleValue to get a Matlab double
            if isequal(self.order_, [])
                self.order_ = self.chain.order;
            end
            o = self.order_;
        end
        
        function O = permOrbits(self)
        % Returns the orbits of the domain 1:self.n under the group
            if isequal(self.permOrbits_, [])
                self.permOrbits_ = qdimsum.group.Orbits.fromGenPerm(self.generators);
            end
            O = self.permOrbits_;
        end
        
        function G = permOrbitRestriction(self, o)
        % Returns the permutation group that acts of the o-th permutation orbit
            if isequal(self.permOrbitRestrictions_, [])
                m = length(self.permOrbits.orbits);
                self.permOrbitRestrictions_ = cell(1, m);
                for oo = 1:m
                    orbit = self.permOrbits.orbits{oo};
                    m = length(orbit);
                    n = self.n;
                    backIndex = zeros(1, n);
                    backIndex(orbit) = 1:m;
                    computeBackIndex = @(el) sign(el).*backIndex(abs(el));
                    generators1 = zeros(0, m);
                    for i = 1:self.nGenerators
                        g = computeBackIndex(self.generators(i, orbit));
                        if ~isequal(g, 1:m)
                            generators1 = [generators1
                                           g];
                        end
                    end
                    self.permOrbitRestrictions_{oo} = qdimsum.Group(generators1);
                end
            end
            G = self.permOrbitRestrictions_{o};
        end

        
        function P = phaseConfiguration(self)
        % Returns the phase configuration generated by the group
            if isequal(self.phaseConfiguration_, [])
                self.phaseConfiguration_ = qdimsum.group.PhaseConfiguration.fromGenPerm(self.generators);
            end
            P = self.phaseConfiguration_;
        end
        
        function [U reps] = irreducibleDecomposition(self, settings)
        % TODO remove
            import qdimsum.*
            [U reps] = Reps.irreducibleDecomposition(self, settings);
        end
        
        function [U reps fromOrbit] = isotypicComponents(self, settings)
function P = phaseConfiguration(self)
        % Returns the phase configuration generated by the group
            if isequal(self.phaseConfiguration_, [])
                self.phaseConfiguration_ = qdimsum.group.PhaseConfiguration.fromGenPerm(self.generators);
            end
            P = self.phaseConfiguration_;
        end
        
        function [U reps] = irreducibleDecomposition(self, settings)
        % TODO remove
            import qdimsum.*
            [U reps] = Reps.irreducibleDecomposition(self, settings);
        end
        
        function [U reps fromOrbit] = isotypicComponents(self, settings)
        % TODO remove
            import qdimsum.*
            [U reps fromOrbit] = Reps.isotypicComponents(self, settings);
        end
        
    end
   
    methods (Static)
        
        function G = trivial(n)
            G = qdimsum.Group(zeros(0, n));
        end
        
        function G = cyclic(d)
        % The cyclic group of degree (and order) d
            G = qdimsum.Group([2:d 1], d);
        end
        
        function G = symmetric(d)
        % The symmetric group of degree d
            if d == 1
                G = qdimsum.Group.trivial(d);
            else
                o = java.math.BigInteger.ONE;
                for i = 2:d
                    o = o.multiply(java.math.BigInteger.valueOf(i));
                end
                G = qdimsum.Group([2:d 1; 2 1 3:d], o);
            end
        end
        
        function G = quaternion
        % The quaternion group of order 8
            import qdimsum.*
            G = Group([GenPerm.fromCycles(8, [1 2 4 7], [3 6 8 5])
                       GenPerm.fromCycles(8, [1 3 4 8], [2 5 7 6])], 8);
        end
        
        function G = binaryTetrahedralGroup
            import qdimsum.*
            g1 = GenPerm.fromCycles(24, [1 2 6], [3 8 20], [4 16 13], [5 9 15], [7 14 10], [11 18 24], [12 23 21], [17 22 19]);
            g2 = GenPerm.fromCycles(24, [1 11 5 3], [2 17 9 7], [4 10 12 19], [6 21 15 13], [8 16 18 23], [14 20 22 24]);
            G = Group([g1;g2], 24);
        end
        
    end
    
end
