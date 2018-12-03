classdef Solution
    
    properties
        realization;
        obj;
        x;
    end
    
    methods
        
        function self = Solution(realization, obj, x)
            self.realization = realization;
            self.obj = obj;
            self.x = x;
        end
        
    end
    
end
