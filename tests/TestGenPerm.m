import qdimsum.*
n = 20;
for i = 1:100
    g1 = GenPerm.random(n);
    g2 = GenPerm.random(n);
    g3 = GenPerm.random(n);
    g12 = GenPerm.compose(g1, g2);
    g23 = GenPerm.compose(g2, g3);
    comp1 = GenPerm.compose(g1, g23);
    comp2 = GenPerm.compose(g12, g3);
    assert(isequal(comp1, comp2));
    ginv1 = GenPerm.inverse(g1);
    id1 = GenPerm.compose(g1, ginv1);
    id2 = GenPerm.compose(ginv1, g1);
    assert(isequal(id1, 1:n));
    assert(isequal(id2, 1:n));
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
