classdef Dec
    
    properties (Abstract, SetAccess = immutable)
        group;
        fromOrbit;
        U;
        compDims;
        repDims;
        repMuls;
        settings;
    end
    
    methods (Abstract)
        
        n = nComponents(self);
        
    end
    
end
