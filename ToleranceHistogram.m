classdef ToleranceHistogram < handle
    
    properties
        title = '';
        logMin = [];
        logMax = [];
        nBins = [];
        counts = [];
        rejected = 0;
    end
    
    methods
        
        function self = ToleranceHistogram(title, logMin, logMax, nBins)
            if nargin < 4
                logMin = -30;
                logMax = 10;
                nBins = 80;
            end
            self.title = title;
            self.logMin = logMin;
            self.logMax = logMax;
            self.nBins = nBins;
            self.counts = zeros(1, nBins);
        end
                
        function register(self, values)
            values = values(:);
            lmin = self.logMin;
            lmax = self.logMax;
            ldiff = lmax - lmin;
            for i = 1:length(values)
                bin = (log10(abs(values(i))) - lmin)/ldiff*self.nBins;
                bin = floor(bin + 1);
                if bin >= 1 && bin <= self.nBins
                    self.counts(bin) = self.counts(bin) + 1;
                else
                    self.rejected = self.rejected + 1;
                end
            end
        end
        
        function edges = binEdges(self)
           edges = self.logMin + (0:self.nBins)/self.nBins * (self.logMax - self.logMin); 
       end
       
        function bm = middlePoints(self)
            bm = self.logMin + ((0:self.nBins-1) + 0.5)/self.nBins*(self.logMax - self.logMin);
        end

        function newHist = coarseGrain(self, div)
            assert(mod(self.nBins, div) == 0, 'The given divisor must divide the bin count');
            newCounts = sum(reshape(self.counts, [div self.nBins/div]), 1);
            newHist = ToleranceHistogram(self.logMin, self.logMax, self.nBins/div);
            newHist.counts = newCounts;
            newHist.rejected = self.rejected;
        end
        
        function h = plot(self, axes, range)
            nonZeroCounts = find(self.counts > 0);
            if ~isequal(nonZeroCounts, [])
                if nargin < 3
                    range = min(nonZeroCounts):max(nonZeroCounts);
                end
                x = self.middlePoints;
                y = self.counts;
                if nargin > 1
                    range = intersect(find(y > 0), range);
                    h = bar(axes, x(range), log10(y(range)));
                    title(axes, self.title);
                    xlabel(axes, 'log10(abs(epsilon))');
                    ylabel(axes, 'log10(Counts)');
                else
                    range = intersect(find(y > 0), range);
                    h = bar(axes, x(range), log10(y(range)));
                    title(self.title);
                    xlabel('log10(abs(epsilon))');
                    ylabel('log10(Counts)');
                end
            end
        end
        
    end
end
