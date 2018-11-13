package com.faacets.qdimsum;

/** Disjoint set forest data structure to compute the orbits of a permutation group
 *
 * See the corresponding Matlab code at qdimsum.group.OrbitsBuilder.
 */
public class OrbitsBuilder {

    public int n;
    public int[] parent;
    public int[] size;
    public int[] index;
    public int nOrbits;
    public int[][] orbits;

    public int[] getOrbit(int o) {
        return orbits[o];
    }

    public static OrbitsBuilder fromPerm(int n, int[][] generators) {
        OrbitsBuilder b = new OrbitsBuilder(n);
        for(int i = 0; i < generators.length; i ++) {
            for (int x = 0; x < n; x ++) {
                int xi = generators[i][x];
                b.union(x, xi);
            }
        }
        b.computeOrbits();
        return b;
    }

    public OrbitsBuilder(int _n) {
        n = _n;
        parent = new int[n];
        size = new int[n];
        index = new int[n];
        nOrbits = n;
        for(int x = 0; x < n; x ++) {
            parent[x] = x;
            size[x] = 1;
            index[x] = -1;
        }
    }

    public int find(int x) {
        while (parent[x] != x) {
            int nxt = parent[x];
            parent[x] = parent[nxt];
            x = nxt;
        }
        return x;
    }

    public void union(int x, int y) {
        int xRoot = find(x);
        int yRoot = find(y);
        if (xRoot < yRoot) {
            parent[yRoot] = xRoot;
            nOrbits --;
            size[xRoot] += size[yRoot];
            size[yRoot] = 0;
        } else if (xRoot > yRoot) {
            parent[xRoot] = yRoot;
            nOrbits --;
            size[yRoot] += size[xRoot];
            size[xRoot] = 0;
        }
    }

    public void computeOrbits() {
        orbits = new int[nOrbits][];
        int[] orbitIndex = new int[nOrbits];
        int orbit = 0;
        for(int x = 0; x < n; x ++) {
            int x0 = find(x);
            if (index[x0] == -1) {
                // new orbit
                index[x0] = orbit;
                orbits[orbit] = new int[size[x0]];
                orbit ++;
            }
            int o = index[x0];
            index[x] = o;
            orbits[o][orbitIndex[o]] = x;
            orbitIndex[o] ++;
        }
    }

}
