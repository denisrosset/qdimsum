classdef IsoDec < handle
    
    properties (GetAccess = public, SetAccess = protected)
        group;    % Generalized permutation group of which we decompose the natural representation
        U;        % Orthonormal change of basis matrix
    end
    
end
