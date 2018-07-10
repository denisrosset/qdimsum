classdef NVSettings < handle & matlab.mixin.SetGet
    
    properties
                                % General settings
                                % ================
                                %
        yalmipSettings = [];    % Settings to pass to YALMIP
        sampleChunkSize = 20;   % How many samples to generate at once
        parallel = false;       % Whether to use parallelism
        checkLevel = 1;         % Level of sanity checks
                                % 0 - disable
                                % 1 - enable, low overhead
                                % 2 - enable, high overhead
        numSamplesChecks = 10;  % number of samples to use in random checks
        verbosityLevel = 1;     % Level of display output
                                % 0 - disable
                                % 1 - warnings and high level information only
                                % 2 - detailed output
                                %
                                % Settings used in randomized checks
                                % ==================================
                                %
        checksEigTol = 1e-9;    % Tolerance for eigenvalues of SDP constraints
        checksTol = 1e-12;      % Tolerance for linear equality constraints
                                %
                                % Settings for monomial samples comparison
                                % ========================================
                                %
                                % Used to remove duplicates when constructing the monomial basis, 
                                % and to discover the action of the symmetry group on monomials
                                %
        monomialTol = 1e-12;    % tolerance used to check for duplicates among monomials
        monomialHist = [];      % histogram of differences below tolerance
                                %
                                % Settings for symmetry group discovery
                                % =====================================
                                %
        findSymTol = 1e-12;     % tolerance used to check for invariant objective value
        findSymHist = [];       % histogram of differences below tolerance
                                %
                                % Settings for numerical block diagonalization
                                % ============================================
                                %
        blockDiagRefine = true; % Whether to use eigenspace refinement
        blockDiagOrbits = true; % Whether to consider group orbits in the basis computation
        blockDiagEigsOpts = []; % Options to pass to "eigs" when performing refinement
        blockDiagEigTol = 1e-9; % tolerance used to compare eigenvalues
        blockDiagEigHist = [];  % histogram of eigenvalue magnitudes below tolerance
        blockDiagMatTol = 1e-9; % tolerance used to compare singular values of blocks to discover structure
        blockDiagMatHist = [];  % histogram of element magnitudes below tolerance
    end

    methods
        
        function self = NVSettings(varargin)
            if nargin > 0
                set(self, varargin{:});
            end
        end
        
        function log(self, message, level)
            if nargin < 3
                level = 1;
            end
            if self.verbosityLevel >= level
                disp(message);
            end
        end
        
        function h = plotHistograms(self)
            figure
            hists = {self.monomialHist self.blockDiagEigHist self.blockDiagMatHist};
            nHists = length(hists);
            allEqual = @(x) all(x(1:end-1) == x(2:end));
            alignAxes = allEqual(cellfun(@(th) th.logMin, hists)) && allEqual(cellfun(@(th) th.logMax, hists)) && ...
                allEqual(cellfun(@(th) th.nBins, hists));
            if alignAxes
                I = cellfun(@(th) sum(th.counts) > 0, hists);
                minBin = min(cellfun(@(th) min(find(th.counts > 0)), hists(I)));
                maxBin = max(cellfun(@(th) max(find(th.counts > 0)), hists(I)));
                range = minBin:maxBin;
            end
            for i = 1:nHists
                ax = subplot(nHists,1,i);
                if alignAxes
                    hists{i}.plot(ax, range);
                else
                    hists{i}.plot(ax);
                end
            end
        end
        
    end
    
    methods (Static)
        
        function yalmipSettings = yalmipMOSEK(tolerance, varargin)
            yalmipSettings = sdpsettings('solver', 'mosek', ...
                                         'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS', tolerance, ...
                                         'mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS', tolerance, ...
                                         'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP', tolerance, ...
                                         'mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS', tolerance, ...
                                         'mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED', tolerance, ...
                                         'mosek.MSK_DPAR_INTPNT_TOL_PFEAS', tolerance, ...
                                         'mosek.MSK_DPAR_INTPNT_TOL_DFEAS', tolerance, ...
                                         'mosek.MSK_DPAR_INTPNT_TOL_INFEAS', tolerance, ...
                                         'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP', max( 1e-14, tolerance), ...
                                         varargin{:});
        end
        
        function yalmipSettings = yalmipSeDuMi(tolerance, varargin)
            yalmipSettings = sdpsettings('solver', 'sedumi', ...
                                         'sedumi.eps', tolerance, ...
                                         'sedumi.cg.qprec', 1, ...
                                         'sedumi.cg.maxiter', 49, ...
                                         'sedumi.stepdif', 2, ...
                                         varargin{:});
        end
        
        function yalmipSettings = yalmipSDPT3(tolerance, varargin)
            yalmipSettings = sdpsettings('solver', 'sdpt3', ...
                                         'sdpt3.vers', 0, ... % restore default
                                         'sdpt3.inftol', 1e-8, ... % was 1e-7
                                         'sdpt3.maxit', 100, ... %
                                         'sdpt3.spdensity', 0.4, ... % was 0.5
                                         'sdpt3.smallblkdim', 50, ... % was 30
                                         'sdpt3.gaptol', tolerance, ...
                                         varargin{:});
        end
        
        function yalmipSettings = yalmipSDPNAL(tolerance, varargin)
            yalmipSettings = sdpsettings('solver', 'sdpnal', ...
                                         'sdpnal.maxiter', 20000, ...
                                         'sdpnal.AAtsolve', 'direct', ...
                                         'sdpnal.tol', tolerance, ...
                                         varargin{:}); % default is 1e-6
                                                       % 'sdpnal.maxitersub', 20, ... does not matter
                                                       % 'sdpnal.stopoption', 0, ... % do not stop prematurely if stagnation
                                                       % 'sdpnal.tolADM', 1e-4, ...
                                                       % 'sdpnal.maxiterADM', 200, ... % from 200 to 400
        end

        function scsSettings = yalmipSCS(tolerance, varargin)
            scs_settings = @(e) sdpsettings('solver', 'scs', ...
                                            'scs.alpha', 1.8, ... % restore
                                            'scs.max_iters', 10000, ...
                                            'scs.eps', tolerance, ...
                                            varargin{:});
        end
        
        function settings = withHistograms(varargin)
            import qdimsum.*                        
            settings = NVSettings(varargin{:});
            settings.monomialHist = ToleranceHistogram('monomialHist: monomial comparison');
            settings.blockDiagEigHist = ToleranceHistogram('blockDiagEigHist: block diagonalization eigenvalues');
            settings.blockDiagMatHist = ToleranceHistogram('blockDiagMatHist: block diagonalization matrix coefficients');
            settings.findSymHist = ToleranceHistogram('findSymHist: invariance of objective value under symmetry');
        end
    
    end
    
end
