classdef Realization

    properties
        relaxation;
        method;
        samples;
        blockStructure;
        objs;
    end
    
    methods
        
        function r = rank(self)
            r = size(self.samples, 2);
        end
        
        function self = Realization(relaxation, method, samples, blockStructure, objs)
            self.relaxation = relaxation;
            self.method = method;
            self.samples = samples;
            self.blockStructure = blockStructure;
            self.objs = objs;
        end
        
        function [S time] = solveWithYalmip(self)
            import qdimsum.*
            r = self.rank;
            x = sdpvar(r - 1, 1);
            nb = self.blockStructure.nBlocks;
            obj = self.objs(1) + (self.objs(2:r) - self.objs(1)) * x;
            CONS = [];
            for b = 1:nb
                range = self.blockStructure.blockRange(b);
                bs = self.blockStructure.blockSize(b);
                C = self.blockStructure.extractBlock(self.samples(:, 1), b);
                d = size(C, 1);
                A = zeros(d, d, r - 1);
                for i = 2:r
                    v = self.samples(range, i) - self.samples(range, 1);
                    A(:,:,i-1) = BlockStructure.vecToMat(v, bs);
                end
                A = C + reshape(reshape(A, [d*d r-1]) * x, [d d]);
                CONS = [CONS; A >= 0];
            end
            tic;
            solvesdp(CONS, -obj);
            time = toc;
            xx = zeros(r, 1);
            xx(2:end) = double(x);
            xx(1) = 1 - sum(xx(2:end));
            S = Solution(self, double(obj), xx);
        end
        
    end
    
end
