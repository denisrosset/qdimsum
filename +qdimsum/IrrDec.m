classdef IrrDec < handle
% TODO: handle complex or quaternionic representations
    
    properties (GetAccess = public, SetAccess = protected)
        group;        % Generalized permutation group of which we decompose the natural representation
        U;            % Orthonormal change of basis matrix
        compDims;     % Isotypic component dimensions
        repDims;      % Representation dimensions
        repMuls;      % Representation multiplicities
                      % all types are real
        repTypes;     % 1 - real, 2 - complex, 3 - quaternionic
    end
    
    methods
       
        function self = IrrDec(group, U, repDims, repMuls, repTypes)
            self.group = group;
            self.U = U;
            self.repDims = repDims(:)';
            self.repMuls = repMuls(:)';
            self.compDims = repDims(:)' .* repMuls(:)';
            self.repTypes = repTypes(:)';
        end
        
        function n = nComponents(self)
            n = length(self.compDims);
        end
        
        function perm = swap(self)
            shift = 0;
            perm = [];
            for r = 1:self.nComponents
                d = self.repDims(r);
                m = self.repMuls(r);
                s = reshape(reshape(1:d*m, [d m])', [1 d*m]);
                perm = [perm, shift + s];
                shift = shift + d*m;
            end
        end
        
    end
    
    methods (Static)
       
        function I = forGroup(group, settings)
            [U reps] = qdimsum.Reps.irreducibleDecomposition(group, settings);
            types = ones(1, size(reps, 2));
            I = qdimsum.IrrDec(group, U, reps(1, :), reps(2, :), types);
        end
        
    end
    
end
