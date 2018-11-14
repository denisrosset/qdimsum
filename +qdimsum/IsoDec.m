classdef IsoDec < handle
    
    properties (GetAccess = public, SetAccess = protected)
        group;        % Generalized permutation group of which we decompose the natural representation
        U;            % Orthonormal change of basis matrix
        compDims;     % Isotypic component dimensions
        repDims;      % Representation dimensions
        repMuls;      % Representation multiplicities
        repTypes;        % 1 - real, 2 - complex or quaternionic
    end
    
    methods
       
        function self = IsoDec(group, U, repDims, repMuls, repTypes)
            self.group = group;
            self.U = U;
            self.repDims = repDims(:)';
            self.repMuls = repMuls(:)';
            self.compDims = repDims(:)' .* repMuls(:)';
            self.repTypes = repTypes(:)';
        end
                
    end
    
    methods (Static)
       
        function I = forGroup(group, settings)
            [U reps] = qdimsum.Reps.isotypicComponents(group, settings);
            shift = 0;
            sample = qdimsum.GenPerm.symmetrize(randn(group.n, group.n), group.decomposition);
            types = zeros(1, size(reps, 2));
            shift = 0;
            for i = 1:size(reps, 2)
                compDim = reps(1, i) * reps(2, i);
                Urep = U(:, shift+(1:compDim));
                t = qdimsum.Reps.identifyRepresentationType(Urep'*sample*Urep, settings);
                switch t
                  case 'R'
                    types(i) = 1;
                  case 'CQ'
                    types(i) = 2;
                end
                shift = shift + compDim;
            end
            I = qdimsum.IsoDec(group, U, reps(1, :), reps(2, :), types);
        end
        
    end
    
end
