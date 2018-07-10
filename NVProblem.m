% Class that describes an optimization problem in noncommutative variables to be solved by the Navascues-Vertesi hierarchy based
% on sampling feasible moment matrices
classdef (Abstract) NVProblem < handle

    properties
        cacheNumOperators = [];  % Number of operators, computed by taking a single sample
    end
    
    properties(Abstract) % Configuration
        
        forceReal;  % Whether to use only the real part of the sampled moment matrix
                    % If the objective is real, this can safely be set to true
        
    end
    
    methods (Abstract) % Problem definition: methods MUST BE implemented
        
        % Returns a 1xI cell array of operators representing a typical sample of
        % the feasible matrix realizations of those operators
        X = sampleOperators(self)
        
        % Returns the state such that elements of the moment matrix are computed as
        % chi(j,k) = K'*mono(j)'*mono(k)*K
        %
        % i.e. tau = K*K' is the density matrix such that 
        % chi(j,k) = trace(tau * mono(j)'*mono(k))
        K = sampleStateKraus(self)

        % Computes the objective value for the given operators and state;
        % please ensure that the returned value is real
        obj = computeObjective(self, X, K)

    end
    
    methods % For family monomial construction: override to use 'families' as a monomial parameter
        
        function types = operatorTypes(self)
        % This method should return a partition of the operator indices
        %
        % For example, for a Bell scenario, return {[indicesAlice] [indicesBob]}
        %
        % See RAC22.m for an example
            types = {1:self.numOperators};
        end

    end
    
    methods % Generators of the symmetry/ambient group (required only for symmetrized optimization)
        
        % Returns the generators of the symmetry group of the problem
        %
        % All elements {g} in the symmetry group must satisfy two conditions
        % 
        % - if X is a valid sample, then g(X) must be a valid sample as well
        % - the objective must be invariant under g: objective(X, K) == objective(g(x), K)
        %
        % Generators are given in a nG x nOp matrix, where each row is a
        % generalized permutation acting on the operator variables.
        %
        % nG is the number of generators, and nOp the number of operators
        % as returns by "sampleOperators"
        %
        % If the symmetry group is trivial, set generators to the identity element
        function generators = symmetryGroupGenerators(self)
            generators = 1:self.numOperators;
        end
        
        % Returns the generators of the ambient group of the problem
        %
        % All elements {g} of the ambient group must leave the feasible set
        % invariant (i.e. if X is a valid sample, then g(X) is a valid sample).
        %
        % The format of the matrix providing those generators is the same as
        % for "symmetryGroupGenerators" above.
        function generators = ambientGroupGenerators(self)
            generators = 1:self.numOperators;
        end
        
    end
    
    methods % Constraint definitions, OPTIONAL
        
        % These methods do not need to be implemented, they are only
        % useful to check the well-behavedness of the user implementation
        
        % From the given feasible realization of operator variables X,
        % computes the terms corresponding to semidefinite positive
        % constraints and returns them in a 1xm cell array of Hermitian
        % operators, each element being a product of the original operator
        % variables.
        %
        % The constraint is that C{i} >=sdp 0, for all i
        %
        % By default, no constraints are defined
        function C = operatorSDPConstraints(self, X)
            C = {};
        end

        % From the given feasible realization of operator variables X,
        % computes the terms corresponding to operator equality constraints
        %
        % The constraint is that C{i} == 0, for all i
        %
        % By default, no constraints are defined
        function C = operatorEqualityConstraints(self, X)
            C = {};
        end
        
        % From the given feasible realization of operator variables,
        % computes the terms corresponding to scalar equality constraints
        %
        % The constraint is that trace(C{i}*rho) == 0, for all i and all rho
        % given by 'sampleState'.
        function C = scalarEqualityConstraints(self, X)
            C = {};
        end
        
    end
    
    
    methods
                        
        function n = numOperators(self)
            if isequal(self.cacheNumOperators, [])
                self.cacheNumOperators = length(self.sampleOperators);
            end
            n = self.cacheNumOperators;
        end
        
    end
   
end
