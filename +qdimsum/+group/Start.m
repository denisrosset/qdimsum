classdef Start < handle
    
    properties
        next;
    end
   
    methods (Static)
        
        function start = emptyChain(n)
            start = qdimsum.group.Start;
            start.next = qdimsum.group.Term(n);
        end
        
    end
end
