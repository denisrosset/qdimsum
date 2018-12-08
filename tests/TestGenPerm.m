import qdimsum.*
n = 20;
for i = 1:100
    g1 = GenPerm.random(n);
    g2 = GenPerm.random(n);
    g3 = GenPerm.random(n);
    % test associativity
    g12 = GenPerm.compose(g1, g2);
    g23 = GenPerm.compose(g2, g3);
    comp1 = GenPerm.compose(g1, g23);
    comp2 = GenPerm.compose(g12, g3);
    assert(isequal(comp1, comp2));
    % test inverse law
    ginv1 = GenPerm.inverse(g1);
    id1 = GenPerm.compose(g1, ginv1);
    id2 = GenPerm.compose(ginv1, g1);
    assert(isequal(id1, 1:n));
    assert(isequal(id2, 1:n));
    % test isomorphism to (unsigned) permutations
    u1 = GenPerm.toUnsigned(g1);
    u2 = GenPerm.toUnsigned(g2);
    u12 = GenPerm.toUnsigned(g12);
    assert(isequal(GenPerm.compose(u1, u2), u12));
    v = rand(n, 1);
    img1 = GenPerm.vectorImage(g12, v);
    img2 = GenPerm.vectorImage(g1, GenPerm.vectorImage(g2, v));
    assert(isequal(img1, img2));
    img3 = GenPerm.orthogonalMatrix(g12) * v;
    assert(isequal(img1, img3));
    M1 = GenPerm.orthogonalMatrix(g1);
    M2 = GenPerm.slowOrthogonalMatrix(g1);
    assert(isequal(M1, M2));
end
