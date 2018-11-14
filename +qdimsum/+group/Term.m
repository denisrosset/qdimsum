classdef Term < qdimsum.group.Chain
% The terminal node in a BSGS chain
    
   methods
       
       function self = Term(n)
           self.n = n;          
       end
       
       function s = sift(self, remaining)
           s = remaining;
       end
       
       function b = isTerm(self)
           b = true;
       end
       
   end
   
end
