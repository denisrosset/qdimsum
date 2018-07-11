classdef Term < qdimsum.Chain
% The terminal node in a BSGS chain
    
   methods
       
       function self = Term(n)
           self.n = n;          
       end
       
       function l = length_(self, acc)
           l = acc;
       end
       
       function r = random_(self, acc)
           r = acc;
       end
       
       function o = order_(self, acc)
           o = acc;
       end
       
       function s = sift(self, remaining)
           s = remaining;
       end
       
       function s = strongGeneratingSet_(self, acc)
           s = acc;
       end
       
       function b = isTerm(self)
           b = true;
       end
       
       function gd = groupDecomposition_(self, acc)
           gd = acc;
       end
       
   end
   
end
