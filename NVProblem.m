% Class that describes an optimization problem in noncommutative variables to be solved by the Navascues-Vertesi hierarchy based
% on sampling feasible moment matrices
classdef (Abstract) NVProblem < handle

    properties(Abstract) % Configuration
        forceReal; % Whether to use only the real part of the sampled moment matrix
                   % If the objective is real, this can safely be set to true
    end
    
    methods (Abstract) % Problem definition
        
        % All these methods need to be implemented

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
    
    methods
                
        % Returns the chain of cosets representatives able to describe
        % any element in the group[
        %
        % chain is a 1xC cell array, where the i-th cell element is
        % a E(i) x n matrix, whose first row is the identity Mon
        % and E(i) is the number of representatives in the i-th cell
        %
        % n is the size of the mon, i.e. the number of operator variables
        %
        % All group elements are written as the composition of a single row
        % from each cell, i.e.
        % 
        % g = chain{1}(i1,:) o chain{2}(i2,:) o chain{3}(i3,:) ...
        % for all combinations of indices (i1,i2,i3, ...) 
        %
        % Note: this function will be called several times by the
        % optimization framework, so if the result is expensive,
        % it should be cached for reuse.
        function chain = groupDecomposition(self)
            chain = 'Not implemented';
        end
        
    end
    
    methods
        
        function tau = sampleState(self)
            K = self.sampleStateKraus;
            tau = K*K';
            tau = (tau + tau')/2;
        end
        
        function types = operatorTypes(self)
            types = 'Not implemented';
        end
                
        function chi = computeMomentMatrixForIndices(self, X, tau, indices)
            monos = cellfun(@(I) Monomials.computeFromIndices(X, I), indices, 'UniformOutput', false);
            N = length(monos);
            d = size(tau, 1);
            % special case: tau is the identity
            % TODO: accept scalar multiples of identity as well
            if isequal(tau, eye(d))
                % conj(S) = (mono').'
                
                colMonos = zeros(d*d, N);
                for i = 1:N
                    mono = monos{i};
                    colMonos(:, i) = mono(:);
                end
                chi = zeros(N, N);
                for j = 1:N
                    chi(j,j:N) = colMonos(:,j)'*colMonos(:,j:N);
                    chi(j+1:N,j) = conj(chi(j,j+1:N));
                end
            else
                for j = 1:N
                    for k = 1:N
                        chi(j,k) = trace(monos{j}'*monos{k}*tau);
                    end
                end
            end
            if self.forceReal
                chi = real(chi);
            end
            % force Hermitian/symmetric if rounding errors etc...
            chi = (chi + chi')/2;
        end
        
        % Compute a single sample for the given problem
        function [chi objContrib] = sample(self)
            N = self.numMonomials;
            X = self.sampleOperators;
            tau = self.sampleState;
            chi = self.computeMomentMatrix(X, tau);
            objContrib = self.computeObjective(X, tau);            
        end

        % Compute a single sample for the given problem
        function [chi objContrib] = sampleForIndices(self, indices)
            X = self.sampleOperators;
            tau = self.sampleState;
            chi = self.computeMomentMatrixForIndices(X, tau, indices);
            objContrib = self.computeObjective(X, tau);            
        end

    end

    methods % Constraint definitions
        
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
   
end
