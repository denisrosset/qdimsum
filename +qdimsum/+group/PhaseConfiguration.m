classdef PhaseConfiguration
    properties
        n;        % matrix size
        phase;    % phase array, containing roots of unity
        orbits;   % cell array of orbits, each of size 2 x orbitSize
                  % note: orbits are transposed compared to the Builders
        index;    % index of the cell for each element, or 0 if element is zero
    end
    
    methods
        
        function self = PhaseConfiguration(n, index, shift, order, orbits)
            import qdimsum.group.PhaseConfiguration
            self.n = n;
            self.index = index;
            self.phase = zeros(n, n);
            for r = 1:n
                for c = 1:n
                    self.phase(r, c) = PhaseConfiguration.rootOfUnity(shift(r, c), order);
                end
            end
            self.orbits = orbits;
        end
        
        function m = nOrbits(self)
        % Returns the number of orbits
            m = length(self.orbits);
        end
                
        function M1 = project(self, M)
            n = self.n;
            M1 = zeros(n, n);
            for i = 1:length(self.orbits)
                o = self.orbits{i};
                lini = o(1,:) + n*(o(2,:)-1);
                f = length(o);
                % undo the phases and average
                s = dot(M(lini), conj(self.phase(lini)))/f;
                % put back the phases
                M(lini) = s * self.phase(lini);
            end
        end
        
        function m = orbitSize(self, o)
        % Returns the size of the o-th orbit
            m = length(self.orbits{o});
        end
        
        function t = orbitTranspose(self, o)
        % Returns the index of the transpose of this orbit, i.e. the orbit
        % corresponding to {(j,i) for all (i,j) in the o-th orbit}
            orbit = self.orbits{o};
            t = self.index(orbit(1, 1), orbit(2, 1));
        end
        
        function t = orbitTransposeScalar(self, o, s)
        % Returns the scalar t corresponding to the transpose of the o-th orbit
        % when the o-th orbit is assigned the scalar s
        %
        % Enables 
        % setOrbit(M, o, s)
        % setOrbit(M, orbitTranspose(o), orbitTransposeScalar(o, s))
            orbit = self.orbits{o};
            r = orbit(1, 1);
            c = orbit(2, 1);
            t = s / self.phase(r, c) * self.phase(c, r);
        end
        
        function M = setOrbit(self, M, o, s)
        % Sets all elements of the o-th orbit in the matrix M to the value s
            orbit = self.orbits{o};
            linIndices = orbit(1,:) + self.n*(orbit(2,:)-1);
            phases = 1 - 2*self.phase(linIndices);
            M(linIndices) = phases * s;
        end
        
        function b = isOrbitDiagonal(self, o)
            orbit = self.orbits{o};
            b = (orbit(1, 1) == orbit(2, 1)); % if the first element is diagonal, all elements are diagonal 
        end
        
        function b = isOrbitSelfTranspose(self, o)
        % Returns whether the o-th orbit contains its transpose, i.e.
        % for every (i,j) in the orbit, (j,i) is also contained
            b = self.orbitTranspose(o) == o;
        end
        
        function M = sampleRealGaussian(self)
        % Samples from a matrix with real entries distributed
        % according to the standard normal distribution, and
        % then symmetrized under the current phase configuration
        %
        % Similar to Random.realGaussian followed by self.project
        %
        % Standard deviation of average of n random normal
        % variables is origStdDev/sqrt(n)
            n = self.n;
            M = zeros(n, n);
            for o = 1:self.nOrbits
                s = randn / sqrt(self.orbitSize(o));
                M = self.setOrbit(M, o, s);
            end
        end
        
        function M = sampleComplexGaussian(self)
        % Similar to sampleRealGaussian except it samples
        % from complex normal entries
        %
        % Standard deviation of real/imag part is 1/sqrt(2)
            n = self.n;
            M = zeros(n, n);
            for o = 1:self.nOrbits
                s = (randn + randn*1i) / sqrt(self.orbitSize(o) * 2);
                M = self.setOrbit(M, o, s);
            end
        end
        
        function M = sampleSymmetricGaussian(self)
        % Works as if it samples from a matrix from the 
        % Gaussian Orthogonal Ensemble
        % and then symmetrizes it under the current phase configuration
        %
        % Similar to Random.symmetricGaussian followed by self.project
            n = self.n;
            M = zeros(n, n);
            for o = 1:self.nOrbits
                if self.isOrbitDiagonal(o)
                    % diagonal elements are scaled up by sqrt(2)
                    s = randn * sqrt(2/self.orbitSize(o));
                    M = self.setOrbit(M, o, s);
                elseif self.isOrbitSelfTranspose(o)
                    % elements are averaged over standard normals
                    s = randn / sqrt(self.orbitSize(o));
                    M = self.setOrbit(M, o, s);
                else
                    ot = self.orbitTranspose(o);
                    if o < ot
                        % the orbit index o is minimal under transpose, so do it
                        % factor 2 due to the transpose part
                        s = randn / sqrt(self.orbitSize(o)*2);
                        st = self.orbitTransposeScalar(o, s);
                        M = self.setOrbit(M, o, s);
                        M = self.setOrbit(M, ot, st);
                    end % else do nothing, will be handled by the orbit transpose
                end
            end
        end
        
        function M = sampleHermitianGaussian(self)
        % Similar to sampleSymmetricGaussian for matrices
        % from the Gaussian Unitary Ensemble
        %
        % Similar to Random.hermitianGaussian followed by self.project
            n = self.n;
            M = zeros(n, n);
            for o = 1:self.nOrbits
                if self.isOrbitDiagonal(o)
                    % diagonal elements are real
                    s = randn / sqrt(self.orbitSize(o));
                    M = self.setOrbit(M, o, s);
                elseif self.isOrbitSelfTranspose(o)
                    % self transpose so real
                    % off-diagonal elements are scaled down by sqrt(2)
                    s = randn / sqrt(2*self.orbitSize(o));
                    M = self.setOrbit(M, o, s);
                else
                    ot = self.orbitTranspose(o);
                    if o < ot
                        % factor 2 for off-diagonal, and factor 2 for
                        % averaging over the orbit transpose
                        s = (randn + randn*1i) / sqrt(self.orbitSize(o)*4);
                        st = self.orbitTransposeScalar(o, s);
                        M = self.setOrbit(M, o, s);
                        M = self.setOrbit(M, ot, conj(st));
                    end
                end
            end
       end
        
        
    end
    
    methods (Static)
        
        function C = fromGenPermMatlab(generators)
            n = size(generators, 2);
            B = qdimsum.group.PhaseConfigurationBuilder.fromGenPerm(generators);
            orbits = cellfun(@(x) x', B.orbits, 'UniformOutput', false);
            C = qdimsum.group.PhaseConfiguration(n, B.index, B.phase, 2, orbits);
        end
        
        function C = fromGenPermJava(generators)
            n = size(generators, 2);
            % convert 1-based indices to 0-based indices
            generators(generators > 0) = generators(generators > 0) - 1;
            B = com.faacets.qdimsum.PhaseConfigurationBuilder.fromGenPerm(n, generators);
            orbits = cell(1, B.nOrbits);
            for i = 1:B.nOrbits
                orbits{i} = double(B.getOrbit(i - 1))' + 1;
            end
            index = double(B.index) + 1;
            phase = double(B.phase);
            C = qdimsum.group.PhaseConfiguration(n, index, phase, 2, orbits);
        end
        
        function C = fromGenPerm(generators)
            if exist('com.faacets.qdimsum.PhaseConfigurationBuilder')
                C = qdimsum.group.PhaseConfiguration.fromGenPermJava(generators);
            else
                warning('Warning: qdimsum.jar not in the Java path, using Matlab code as fallback.');
                C = qdimsum.group.PhaseConfiguration.fromGenPermMatlab(generators);
            end
        end
        
        function u = rootOfUnity(k, n)
            switch n
              case 1
                u = 1;
              case 2
                switch k
                  case 0
                    u = 1;
                  case 1
                    u = -1;
                end
              case 4
                switch k
                  case 0
                    u = 1;
                  case 1
                    u = 1i;
                  case 2
                    u = -1;
                  case 3
                    u = -1i;
                end
              otherwise
                u = exp(2*pi*1i*k/n);
            end
        end
    end


end
