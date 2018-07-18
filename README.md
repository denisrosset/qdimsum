# QDimSum

Symmetric SDP relaxations for qudits systems

The code requires YALMIP and a semidefinite solver of your choice (for example Mosek or SeDuMi).

See [RAC22.m](https://github.com/denisrosset/qdimsum/blob/master/RAC22.m) and [TestRAC22.m](https://github.com/denisrosset/qdimsum/blob/master/TestRAC22.m) for a simple example.


## Limitations

- Code should work on complex moment matrices, but this has not been tested. (Please contact us with a use case!).
- Only decomposition of representations when all irreps are real (complex and quaternionic irreps work in progress).

## Future tutorial

2. In the meantime, unpack the ZIP file attached with the latest version.

3. To use the package with Matlab, you need to add the base path (say /qdimsum-master) to Matlab, but not the subfolders (please check what happens if you do: does this fail? please document it in a way that normal Matlab users can understand -- I have a different environment).
You can achieve the steps below by cut'n'pasting stuff from the existing "RAC22.m" and "TestRAC22.m". The idea is to have a tutorial that introduces step by step our conventions, in a way that the user can test the code as soon as possible.

Please check that you can follow these steps.

a. Include only the methods "sampleOperators", "sampleStateKraus", "computeObjective", and set "forceReal = true". Optimize over that first example with monomials = {'npa' 2} (because "operatorTypes" is not defined), and "method='none'" because no symmetry group is defined.

b. Discuss how to split operators into types using "operatorTypes", and how to use the {'families' ...} notation.

c. Then discuss the symmetry group notation, and add the "symmetryGroupGenerators" method, discuss the different methods "reynolds", "isotypic", "irreps", "blocks".

d. Implement "ambientGroupGenerators", and call "findSymmetryGroupGenerators(RAC22, NVSettings)" to discover automatically the symmetries.

e. Implement the constraints "operatorSDPConstraints", "operatorEqualityConstraints", "scalarEqualityConstraints", introduce a deliberate mistake in the sampling and check whether the library detects it.

To write a recap of that, simply use Markdown, the text format that Github understands. A markdown file is a text file with additional tricks. Use the editor https://stackedit.io/app , and the syntax below for code and $x+1$ (inline) or $$x+1$$ (own line) for LaTeX.
```matlab
classdef RAC22 < NVProblem

    properties
       forceReal = true;
    end
   
    methods

        function X = sampleOperators(self, rank)
        % the first 4 operators are the states
        % the next 2 represent a projective measurement
        % the last 2 represent another projective measurement
            dim = 2;
            X = cell(1, 8);
            for i = 1:4
                X{i} = qdimsum.Random.pureNormalizedDensityMatrix(dim);
            end
            U = qdimsum.Random.unitary(2);
            X{5} = U*[1 0; 0 0]*U';
            X{6} = U*[0 0; 0 1]*U';
            U = qdimsum.Random.unitary(2);
            X{7} = U*[1 0; 0 0]*U';
            X{8} = U*[0 0; 0 1]*U';
        end
       
        function K = sampleStateKraus(self)
            dim = 2;
            K = eye(dim);
        end
       
        function obj = computeObjective(self, X, K)
            obj = 0;
            for x1 = 1:2
                for x2 = 1:2
                    for y = 1:2
                        if y == 1
                            b = x1;
                        else
                            b = x2;
                        end
                        rho = X{x1+(x2-1)*2};
                        M = X{4+b+(y-1)*2};
                        obj = obj + trace(M * rho); % K is always identity
                    end
                end
            end
            obj = real(obj / 8);
        end

    end
end
```
