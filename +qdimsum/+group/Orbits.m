classdef Orbits
   
    properties
        n;      % domain size
        index;  % index of the cell for each element
        orbits; % cell array of orbits
    end
    
    methods
        
        function self = Orbits(n, index, orbits)
            self.n = n;
            self.index = index;
            self.orbits = orbits;
        end
        
    end

    methods (Static)
        
        function O = fromGenPermMatlab(generators)
            B = qdimsum.group.OrbitsBuilder.fromPerm(abs(generators));
            O = qdimsum.group.Orbits(B.n, B.index, B.orbits);
        end
        
        function O = fromGenPermJava(generators)
            n = size(generators, 2);
            generators = abs(generators) - 1;
            B = com.faacets.qdimsum.OrbitsBuilder.fromPerm(n, generators);
            orbits = cell(1, B.nOrbits);
            for i = 1:B.nOrbits
                orbit = double(B.getOrbit(i - 1)) + 1;
                orbits{i} = orbit(:)';
            end
            index = double(B.index) + 1;
            O = qdimsum.group.Orbits(B.n, index, orbits);
        end
        
        function O = fromGenPerm(generators)
            if exist('com.faacets.qdimsum.OrbitsBuilder')
                O = qdimsum.group.Orbits.fromGenPermJava(generators);
            else
                warning('Warning: qdimsum.jar not in the Java path, using Matlab code as fallback.');
                O = qdimsum.group.Orbits.fromGenPermMatlab(generators);
            end
            
        end
        
    end

end
