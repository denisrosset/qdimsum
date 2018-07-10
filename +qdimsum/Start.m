classdef Start < handle
    
    properties
        next;
    end
   
    methods (Static)
        
        function start = emptyChain(n)
            import qdimsum.*
            start = Start;
            start.next = Term(n);
        end
        
    end
end
