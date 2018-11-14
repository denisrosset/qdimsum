import qdimsum.*
settings = NVSettings;
% Dihedral group D8
G = Group([2 3 4 1; 3 2 1 4]);
I = IsoDec.forGroup(G, settings);
assert(isequal(sort(I.compDims), [1 1 2]));
% Cyclic group of order 4
G = Group([2:16 1]);
I = IsoDec.forGroup(G, settings);
assert(isequal(sort(I.compDims), [1 1 2 2 2 2 2 2 2]));
assert(isequal(sort(I.repTypes), [1 1 2 2 2 2 2 2 2]));
% Klein four group
G = Group([2 1 3 4; 1 2 4 3]);
I = IsoDec.forGroup(G, settings);
assert(isequal(sort(I.compDims), [1 1 2]));
assert(isequal(sort(I.repTypes), [1 1 1]));

% Mathieu group M9
g1 = GenPerm.fromCycles(9, [1,4,9,8], [2,5,3,6]);
g2 = GenPerm.fromCycles(9, [1,6,5,2], [3,7,9,8]);
G = Group([g1;g2]);
I = IsoDec.forGroup(G, settings);
assert(isequal(sort(I.compDims), [1 8]));
assert(isequal(sort(I.repTypes), [1 1]));

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
assert(isequal(sort(I.repTypes), [1 1 1 1 2]));
