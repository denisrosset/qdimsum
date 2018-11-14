classdef IsoDec < handle
    
    properties (GetAccess = public, SetAccess = protected)
        group;        % Generalized permutation group of which we decompose the natural representation
        U;            % Orthonormal change of basis matrix
        compDims;     % Isotypic component dimensions
        repDims;      % Representation dimensions
        repMuls;      % Representation multiplicities
        types;        % 0 - real, 1 - complex or quaternionic
    end
    
    methods
       
        function self = IsoDec(group, U, repDims, repMuls, types)
            self.group = group;
            self.U = U;
            self.repDims = repDims(:)';
            self.repMuls = repMuls(:)';
            self.compDims = repDims(:)' .* repMuls(:)';
            self.types = types(:)';
        end
        
        function P = componentProjector(self)
            
        end
        
    end
    
end
