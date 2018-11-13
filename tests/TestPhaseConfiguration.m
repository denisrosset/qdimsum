import qdimsum.*
C = PhaseConfiguration.fromGenPerm(4, [2 3 -4 1; 2 1 4 3]);
[M P] = C.toPhasePartition;
Mexp = [1 2 0 2
        2 1 2 0
        0 2 1 2
        2 0 2 1];
Pexp = [0 0 1 0
        0 0 0 1
        1 0 0 1
        0 1 1 0];
assert(isequal(M, uint32(Mexp)));
assert(isequal(P, uint32(Pexp)));

% from I3322 npa level 2
G = [1 5 6 7 2 3 4 23 24 10 15 20 25 26 11 16 21 27 28 12 17 22 8 9 13 14 18 19
     1 3 2 4 5 6 -7 13 14 15 16 -17 8 9 10 11 -12 19 18 20 21 -22 23 -24 25 -26 -27 -28
     1 2 3 -4 6 5 7 8 -9 11 10 12 13 -14 16 15 17 -18 -19 -21 -20 -22 25 26 23 24 28 27];
C = PhaseConfiguration.fromGenPerm(28, G);
[M P] = C.toPhasePartition;
Mexp = [1 2 2 0 2 2 0 3 0 4 4 5 3 0 4 4 5 0 0 5 5 0 3 0 3 0 0 0
        6 7 8 0 9 9 10 11 0 12 12 13 14 0 15 15 16 0 0 17 17 0 18 19 18 19 20 20
        6 8 7 0 9 9 10 14 0 15 15 16 11 0 12 12 13 0 0 17 17 0 18 19 18 19 20 20
        0 0 0 21 22 22 0 0 23 24 24 0 0 23 24 24 0 25 25 26 26 0 27 0 27 0 0 0
        6 9 9 10 7 8 0 18 19 12 15 17 18 19 12 15 17 20 20 13 16 0 11 0 14 0 0 0
        6 9 9 10 8 7 0 18 19 15 12 17 18 19 15 12 17 20 20 16 13 0 14 0 11 0 0 0
        0 22 22 0 0 0 21 27 0 24 24 26 27 0 24 24 26 0 0 0 0 0 0 23 0 23 25 25
        28 29 30 0 31 31 32 33 0 34 34 35 36 0 37 37 38 0 0 39 39 0 40 41 40 41 42 42
        0 0 0 43 44 44 0 0 45 46 46 0 0 47 48 48 0 49 50 51 51 52 53 54 53 54 55 55
        56 57 58 59 57 58 59 60 61 62 63 64 65 66 63 67 68 69 70 64 68 71 60 61 65 66 69 70
        56 57 58 59 58 57 59 60 61 63 62 64 65 66 67 63 68 69 70 68 64 71 65 66 60 61 70 69
        72 73 74 0 75 75 76 77 0 78 78 79 80 0 81 81 82 0 0 83 83 0 84 85 84 85 86 86
        28 30 29 0 31 31 32 36 0 37 37 38 33 0 34 34 35 0 0 39 39 0 40 41 40 41 42 42
        0 0 0 43 44 44 0 0 47 48 48 0 0 45 46 46 0 50 49 51 51 52 53 54 53 54 55 55
        56 58 57 59 57 58 59 65 66 63 67 68 60 61 62 63 64 70 69 64 68 71 60 61 65 66 69 70
        56 58 57 59 58 57 59 65 66 67 63 68 60 61 63 62 64 70 69 68 64 71 65 66 60 61 70 69
        72 74 73 0 75 75 76 80 0 81 81 82 77 0 78 78 79 0 0 83 83 0 84 85 84 85 86 86
        0 0 0 87 88 88 0 0 89 90 90 0 0 91 92 92 0 93 94 95 95 96 97 98 97 98 99 99
        0 0 0 87 88 88 0 0 91 92 92 0 0 89 90 90 0 94 93 95 95 96 97 98 97 98 99 99
        72 75 75 76 73 74 0 84 85 78 81 83 84 85 78 81 83 86 86 79 82 0 77 0 80 0 0 0
        72 75 75 76 74 73 0 84 85 81 78 83 84 85 81 78 83 86 86 82 79 0 80 0 77 0 0 0
        0 0 0 0 0 0 0 0 100 101 101 0 0 100 101 101 0 102 102 0 0 103 0 100 0 100 102 102
        28 31 31 32 29 30 0 40 41 34 37 39 40 41 34 37 39 42 42 35 38 0 33 0 36 0 0 0
        0 44 44 0 0 0 43 53 54 46 48 51 53 54 46 48 51 55 55 0 0 52 0 45 0 47 49 50
        28 31 31 32 30 29 0 40 41 37 34 39 40 41 37 34 39 42 42 38 35 0 36 0 33 0 0 0
        0 44 44 0 0 0 43 53 54 48 46 51 53 54 48 46 51 55 55 0 0 52 0 47 0 45 50 49
        0 88 88 0 0 0 87 97 98 90 92 95 97 98 90 92 95 99 99 0 0 96 0 89 0 91 93 94
        0 88 88 0 0 0 87 97 98 92 90 95 97 98 92 90 95 99 99 0 0 96 0 91 0 89 94 93];

Pexp = [ 0 0 0 1 0 0 1 0 1 0 0 0 0 1 0 0 1 1 1 0 1 1 0 1 0 1 1 1
         0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 1 0 1 1 0 0 0 0 0 0
         0 0 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 1 1 0 1 1 0 1 0 1 1 1
         1 1 1 0 0 1 1 1 0 0 1 1 1 0 0 1 1 0 0 0 0 1 0 1 1 1 1 1
         0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 0 1 1 1
         0 0 0 1 0 0 1 0 1 0 0 0 0 1 0 0 1 1 1 1 1 1 0 1 0 1 1 1
         1 0 1 1 1 1 0 0 1 0 0 0 1 1 1 1 0 1 1 1 1 1 1 0 1 0 0 0
         0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 1 0 1 1 0 0 0 0 0 0
         1 1 1 0 0 1 1 1 0 0 1 1 1 0 0 1 1 0 0 0 0 0 0 0 1 1 0 1
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
         0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 1 1 1 1 0 0 0 0 0 0
         0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 1 0 1 1 0 0 0 0 0 0
         0 0 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 1 1 0 1 1 0 1 0 1 1 1
         1 1 1 0 0 1 1 1 0 0 1 1 1 0 0 1 1 0 0 0 0 1 0 1 1 0 1 0
         0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 1 0 1 1 1
         0 0 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 1 1 1 1 0 0 1 0 1 1 1
         1 1 1 1 1 1 0 1 1 1 1 0 1 1 1 1 0 1 1 1 0 1 1 0 1 0 0 0
         1 1 1 0 0 1 1 1 0 0 1 1 1 0 0 1 1 0 0 0 0 0 0 0 1 1 0 1
         1 1 1 0 0 1 1 1 0 0 1 1 1 0 0 1 1 0 0 0 0 1 0 1 1 0 1 0
         0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 0 1 1 1
         1 1 1 0 1 1 1 1 0 1 1 1 1 0 1 1 0 0 0 0 0 1 1 1 1 1 1 1
         1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 0 1 0 1 1 1 0 1 0 1 1 0 1
         0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 1 0 1 1 1
         1 0 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 1 1 1 0 1 0 1 0 0 0
         0 0 0 1 0 0 1 0 1 0 0 0 0 1 0 0 1 1 1 1 1 1 0 1 0 1 1 1
         1 0 1 1 1 1 0 0 1 0 0 0 1 0 1 1 0 1 0 1 1 1 1 0 1 0 0 0
         1 0 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 1 1 1 0 1 0 1 0 0 0
         1 0 1 1 1 1 0 0 1 0 0 0 1 0 1 1 0 1 0 1 1 1 1 0 1 0 0 0];
assert(isequal(M, uint32(Mexp)));
assert(isequal(P, uint32(Pexp)));