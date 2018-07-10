classdef Start < handle
    
    properties
        next;
    end
   
    methods (Static)
        
        function start = emptyChain(n)
            start = Start;
            start.next = Term(n);
        end
        
    end
end
