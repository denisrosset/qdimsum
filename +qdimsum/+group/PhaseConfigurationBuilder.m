% Constructs a phase configuration using a variant of a disjoint
% set forest that implements path splitting, but not union by rank/size
%
% A phase configuration represent a family (algebra?) of matrices of size n x n
% obeying relations of the form
% M(xr, xc) = M(yr, yc) * phase
%
% The phase is expressed using roots of unity exp(2*pi*i*k/maxPhase),
% "maxPhase". For example, the roots {-1, 1} can be expressed with
% maxPhase = 2.
%
% The data structure is a disjoint set forest that stores the parent
% of each element in the matrix.
%
% The relation above
% M(xr, xc) = M(yr, yc) * phase
% is encoded using yr = parentRow(xr, xc) and yc = parentCol(xr, xc)
% and phase = exp(2*pi*i*phase(xr, xc)/maxPhase)
%
% An element (r, c) can either 
% - be a root, in which case parentRow(r, c) = r and parentCol(r, c) = c
%   and size(r, c) gives the size of the cell of which (r, c) is a
%   representative
% - not be a root, in which case size(r, c) = 0
%
% After all relations have been established, one calls "computeOrbits"
% which provides a numbering of all cells in the matrix "index",
% and a list of orbits in the cell array "orbits"
%
% The cell corresponding to zero elements is encoded as other cells,
% but its root is stored in (zr, zc). Additionally, for all zero elements
% we have phase(r, c) = 0
classdef PhaseConfigurationBuilder < handle
    
    properties
        n;         % matrix size
        maxPhase;  % maximal phase
        parentRow; % parent row index of size n x n
        parentCol; % parent col index of size n x n
        phase;     % phase index of size n x n, 0 ... maxPhase-1
        size;      % size of cell, n x n
        nOrbits;   % number of orbits
        zr;        % (zr, zc) is the root corresponding to
        zc;        % a zero element
        
        % Properties computed by "computeOrbits"
        index;     % cell index,   n x n
        orbits;    % cell array of orbits of size orbitSize x 2
    end
    
    methods (Static)
        
        function C = fromGenPerm(generators)
            n = size(generators, 2);
            C = qdimsum.group.PhaseConfigurationBuilder(n, 2);
            for i = 1:size(generators, 1)
                g = generators(i, :);
                for xr = 1:n
                    for xc = 1:n
                        yr = g(xr);
                        yc = g(xc);
                        yp = double((yr > 0 && yc < 0) || (yr < 0 && yc > 0));
                        C.union(xr, xc, abs(yr), abs(yc), yp); % beware of inverse when generalizing
                    end
                end
            end
            C.computeOrbits;
        end
        
    end
    
    methods
        
        function self = PhaseConfigurationBuilder(n, maxPhase)
            self.n = n;
            self.maxPhase = maxPhase;
            self.parentRow = zeros(n, n);
            self.parentCol = zeros(n, n);
            self.phase = zeros(n, n);
            self.size = ones(n, n);
            self.index = zeros(n, n);
            self.nOrbits = n * n;
            self.zr = 0;
            self.zc = 0;
            for i = 1:n
                self.parentRow(i, :) = i;
                self.parentCol(:, i) = i;
            end
        end
        
        function computeOrbits(self)
            self.orbits = cell(1, self.nOrbits);
            orbitIndex = ones(1, self.nOrbits);
            orbit = 1;
            for xr = 1:self.n
                for xc = 1:self.n
                    [xr0, xc0, ~] = self.find(xr, xc);
                    if xr0 ~= self.zr || xc0 ~= self.zc
                        % nonzero element
                        if self.index(xr0, xc0) == 0
                            % new orbit
                            self.index(xr0, xc0) = orbit;
                            self.orbits{orbit} = zeros(self.size(xr0, xc0), 2);
                            orbit = orbit + 1;
                        end
                        o = self.index(xr0, xc0);
                        self.index(xr, xc) = o;
                        orb = self.orbits{o};
                        orb(orbitIndex(o), :) = [xr xc];
                        self.orbits{o} = orb;
                        orbitIndex(o) = orbitIndex(o) + 1;
                    end
                end
            end
        end
                  
        % Finds the representative of the cell to which (xr, xc) belongs
        % and returns the representative (yr, yc) and the multiplicative sign
        % such that M(xr, xc) = M(yr, yc) * rootOfUnity(yp, maxPhase)
        function [yr yc p] = find(self, xr, xc)
            p = 0;
            while 1
                % path splitting variant
                yr = self.parentRow(xr, xc);
                yc = self.parentCol(xr, xc);
                yp = self.phase(xr, xc);
                if (xr == yr) && (xc == yc)
                    % x.parent == x
                    return
                end
                self.parentRow(xr, xc) = self.parentRow(yr, yc);
                self.parentCol(xr, xc) = self.parentCol(yr, yc);
                np = yp + self.phase(yr, yc);
                if np >= self.maxPhase
                    np = np - self.maxPhase;
                end
                self.phase(xr, xc) = np;
                xr = yr;
                xc = yc;
                p = p + yp;
                if p >= self.maxPhase
                    p = p - self.maxPhase;
                end
            end
        end
        
        function nz = nZeros(self)
            if self.zr == 0 || self.zc == 0
                nz = 0;
            else
                nz = self.size(self.zr, self.zc);
            end
        end
        
        function mergeTo(self, xr0, xc0, yr0, yc0, phase)
        % Merge the root (xr0, xc0) to (yr0, yc0)
        % so that M(xr0, xc0) = M(yr0, yc0) * phase
        % By convention, we require (xr0, xc0) > (yr0, yc0)
            assert(xr0 ~= yr0 || xc0 ~= yc0);
            self.parentRow(xr0, xc0) = yr0;
            self.parentCol(xr0, xc0) = yc0;
            self.phase(xr0, xc0) = phase;
            self.size(yr0, yc0) = self.size(yr0, yc0) + self.size(xr0, xc0);
            self.size(xr0, xc0) = 0;
            self.nOrbits = self.nOrbits - 1;
        end
        
        function phase = normalizePhase(self, phase)
            while phase < 0
                phase = phase + self.maxPhase;
            end
            while phase >= self.maxPhase
                phase = phase - self.maxPhase;
            end
        end
        
        function c = compareIndices(self, xr, xc, yr, yc)
            if xr < yr || (xr == yr && xc < yc)
                c = -1;
            elseif xr > yr || (xr == yr && xc > yc)
                c = 1;
            else
                c = 0;
            end                
        end
        
        % Merges the sets from which (xr, xc) and (yr, yc) are members
        % with a phase difference such that M(xr, xc) = M(yr, yc) * rootOfUnity(p, maxPhase)
        function union(self, xr, xc, yr, yc, p)
            [xr0 xc0 xp0] = self.find(xr, xc);
            [yr0 yc0 yp0] = self.find(yr, yc);
            % phase difference between roots, i.e.
            % M(xr0, xc0) = M(yr0, yc0) * rootOfUnity(p0, maxPhase)
            p0 = self.normalizePhase(yp0 + p - xp0);
            if (xr0 == yr0 && xc0 == yc0 && p0 ~= 0)
                % new zero element discovered
                if self.zr == 0 || self.zc == 0
                    % no zero element known until now
                    self.zr = xr0;
                    self.zc = xc0;
                    self.phase(self.zr, self.zc) = 0;
                    % the zero orbit does not count
                    self.nOrbits = self.nOrbits - 1;
                    return
                else
                    % zero element present, so merge the current cell
                    % with the zero cell
                    salf.phase(xr0, xc0) = 0;
                    switch self.compareIndices(self.zr, self.zc, xr0, xc0)
                      case -1
                        % the zero representative stays the same
                        self.mergeTo(xr0, xc0, self.zr, self.zc, 0);
                      case 1
                        % the new zero cell has the zero representative
                        self.mergeTo(self.zr, self.zc, xr0, xc0, 0);
                        self.zr = xr0;
                        self.zc = xc0;
                    end
                end
            else
                switch self.compareIndices(xr0, xc0, yr0, yc0)
                  case -1
                    self.mergeTo(yr0, yc0, xr0, xc0, self.normalizePhase(-p0));
                  case 1
                    self.mergeTo(xr0, xc0, yr0, yc0, p0);
                end
            end
        end
        
    end
        
end
