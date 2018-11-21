import qdimsum.*
settings = NVSettings;
% Dihedral group D8
G = Group([2 3 4 1; 3 2 1 4]);
I = IsoDec.forGroup(G, settings);
assert(isequal(sort(I.compDims), [1 1 2]));
I.refine.check;

% Cyclic group of order 4
G = Group([2:16 1]);
I = IsoDec.forGroup(G, settings);
assert(isequal(sort(I.compDims), [1 1 2 2 2 2 2 2 2]));
I.refine.check;

% Klein four group
G = Group([2 1 3 4; 1 2 4 3]);
I = IsoDec.forGroup(G, settings);
assert(isequal(sort(I.compDims), [1 1 2]));
I.refine.check;

% Mathieu group M9
g1 = GenPerm.fromCycles(9, [1,4,9,8], [2,5,3,6]);
g2 = GenPerm.fromCycles(9, [1,6,5,2], [3,7,9,8]);
G = Group([g1;g2]);
I = IsoDec.forGroup(G, settings);
assert(isequal(sort(I.compDims), [1 8]));
I.refine.check;

% Rubik cube group

g1 = GenPerm.fromCycles(48, [1, 3, 8, 6], [2, 5, 7, 4], [9,33,25,17], [10,34,26,18], [11,35,27,19]);
g2 = GenPerm.fromCycles(48, [9,11,16,14], [10,13,15,12], [1,17,41,40], [4,20,44,37], [6,22,46,35]);
g3 = GenPerm.fromCycles(48, [17,19,24,22], [18,21,23,20], [6,25,43,16], [7,28,42,13], [8,30,41,11]);
g4 = GenPerm.fromCycles(48, [25,27,32,30], [26,29,31,28], [3,38,43,19], [5,36,45,21], [8,33,48,24]);
g5 = GenPerm.fromCycles(48, [33,35,40,38], [34,37,39,36], [3, 9,46,32], [2,12,47,29], [1,14,48,27]);
g6 = GenPerm.fromCycles(48, [41,43,48,46], [42,45,47,44], [14,22,30,38], [15,23,31,39], [16,24,32,40]);
G = Group([g1;g2;g3;g4;g5;g6]);
I = IsoDec.forGroup(G, settings);
assert(isequal(sort(I.compDims), [2 7 11 12 16]));
I.refine.check;

G = Group.binaryTetrahedralGroup;
% Representations
% 1 trivial
% d=1 complex and its conjugate
% d=2 quaternionic
% d=2 complex and its conjugate
% d=3 real


%\hat{n}_j is the isotypic component size = \overline{m}_j \overline{n}_j


g1 = GenPerm.fromCycles(55, [1,37], [3,39], [4,46], [5,55], [6,48], [7,16], [8,43], [9,18], [10,52], [12,42], [14,51], ...
                       [15,21], [17,41], [19,50], [22,29], [23,45], [24,31], [25,54], [26,35], [27,34], [28,33], ...
                       [30,44], [32,53], [40,47]);
g2 = GenPerm.fromCycles(55, [1,10,13,36], [2,25], [3,26,11,32], [4,6,24,22], [5,9,27,29], [7,14,31,23], [8,15,34,30], ...
                        [12,33], [16,19,18,28], [17,20,21,35], [37,40,39,44], [38,41,42,45], [46,49,48,53], [47,50,51,54]);
G = Group([g1;g2]);