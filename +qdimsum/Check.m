classdef Check

    methods (Static)
        
        function problem(problem, settings)
            import qdimsum.*
            Check.samples(problem, settings);
            Check.symmetryGroup(problem, settings);
            Check.ambientGroup(problem, settings);
        end
        
        function samples(problem, settings)
            import qdimsum.*
            X = problem.sampleOperators;
            K = problem.sampleStateKraus;
            Check.sampleObeysConstraints(problem, X, K, settings, '');
        end

        function symmetryGroup(problem, settings)
            import qdimsum.*
            Check.generatorsLeaveObjectiveInvariant(problem, problem.symmetryGroupGenerators, settings);
            Check.generatorsLeaveOperatorsFeasible(problem, problem.symmetryGroupGenerators, settings);
        end
        
        function ambientGroup(problem, settings)
            import qdimsum.*
            Check.generatorsLeaveOperatorsFeasible(problem, problem.ambientGroupGenerators, settings);
        end

        function generatorsLeaveObjectiveInvariant(problem, generators, settings)
            import qdimsum.*
            X = problem.sampleOperators;
            K = problem.sampleStateKraus;
            obj = problem.computeObjective(X, K);
            for i = 1:size(generators, 1)
                g = generators(i,:);
                X1 = GenPerm.operatorsImage(g, X);
                obj1 = problem.computeObjective(X1, K);
                diff = abs(obj1 - obj);
                assert(diff < settings.checksTol, ...
                       sprintf('The generator %s does not preserve the group', num2str(g)));
            end
        end
        
        function generatorsLeaveOperatorsFeasible(problem, generators, settings)
            import qdimsum.*
            X = problem.sampleOperators;
            K = problem.sampleStateKraus;
            for i = 1:size(generators, 1)
                g = generators(i, :);
                message = sprintf('For generator %s:', num2str(g));
                Check.sampleObeysConstraints(problem, GenPerm.operatorsImage(g, X), K, settings, message);
            end
        end
                
        function sampleObeysConstraints(problem, X, K, settings, message)
            import qdimsum.*
            Csdp = problem.operatorSDPConstraints(X);
            CopEq = problem.operatorEqualityConstraints(X);
            CscEq = problem.scalarEqualityConstraints(X);
            for i = 1:length(Csdp)
                mineig = min(eig(Csdp{i}));
                assert(mineig > -settings.checksEigTol, ...
                       sprintf('%s%d th SDP constraint violated, minimal eigenvalue %f', message, i, mineig));
            end
            for i = 1:length(CopEq)
                maxAbsSingVal = max(abs(eig(CopEq{i})));
                assert(maxAbsSingVal < settings.checksEigTol, ...
                       '%s%d th operator equality constraint violated, maximal singular value %f', message, i, maxAbsSingVal);
            end
            for i = 1:length(CscEq)
                diff = abs(trace(K'*CscEq{i}*K));
                assert(diff < settings.checksTol, ...
                       '%s%d th scalar equality constraint violated, difference %f', message, i, diff);
            end
        end
                
    end
    
end
