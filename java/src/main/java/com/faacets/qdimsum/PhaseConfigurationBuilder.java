package com.faacets.qdimsum;

/** See Matlab version in package qdimsum.group.PhaseConfigurationBuilder
 *
 * This version uses 0-based indices, and encodes a sign flip by storing
 * ~permImage in the generalized permutation image, instead of -permImage
 * as in the Matlab version.
 * */
public class PhaseConfigurationBuilder {

    public int n;
    public int[][] shift;
    public int order;
    public int[][] parentRow;
    public int[][] parentCol;
    public int[][] size;
    public int[][] index;
    public int nOrbits;
    public int[] orbitStart;
    public int[] orbitRow;
    public int[] orbitCol;
    public int zr;
    public int zc;

    public static PhaseConfigurationBuilder fromGenPerm(int n, int[][] generators) {
        PhaseConfigurationBuilder b = new PhaseConfigurationBuilder(n, 2);
        for(int i = 0; i < generators.length; i ++) {
            for (int xr = 0; xr < n; xr ++) {
                for (int xc = 0; xc < n; xc ++) {
                    int yr = generators[i][xr];
                    int yc = generators[i][xc];
                    int p = 0;
                    if (yr < 0) {
                        p += 1;
                        yr = ~yr;
                    }
                    if (yc < 0) {
                        p += 1;
                        yc = ~yc;
                    }
                    b.union(xr, xc, yr, yc, p % 2);
                }
            }
        }
        b.computeOrbits();
        return b;
    }

    public void computeOrbits() {
        int [][][] orbits = new int[nOrbits][][];
        int[] orbitIndex = new int[nOrbits + 1];
        int orbit = 0;
        int[] x0 = new int[3];
        for (int xr = 0; xr < n; xr ++) {
            for (int xc = 0; xc < n; xc ++) {
                find(xr, xc, x0);
                int xr0 = x0[0];
                int xc0 = x0[1];
                if (xr0 != zr || xc0 != zc) {
                    // nonzero element
                    if (index[xr0][xc0] == -1) {
                        // new orbit
                        index[xr0][xc0] = orbit;
                        orbits[orbit] = new int[][] { new int[size[xr0][xc0]], new int[size[xr0][xc0]] };
                        orbit ++;
                    }
                    int o = index[xr0][xc0];
                    index[xr][xc] = o;
                    orbits[o][0][orbitIndex[o]] = xr;
                    orbits[o][1][orbitIndex[o]] = xc;
                    orbitIndex[o] ++;
                }
            }
        }
        orbitStart = new int[nOrbits+1];
        orbitRow = new int[n*n];
        orbitCol = new int[n*n];
        int i = 0;
        for (int o = 0; o < nOrbits; o ++) {
            orbitStart[o] = i;
            int orbitSize = orbits[o][0].length;
            System.arraycopy(orbits[o][0], 0, orbitRow, i, orbitSize);
            System.arraycopy(orbits[o][1], 0, orbitCol, i, orbitSize);
            i += orbitSize;
        }
        orbitStart[nOrbits] = n*n;
    }

    public PhaseConfigurationBuilder(int _n, int _maxPhase) {
        n = _n;
        order = _maxPhase;
        parentRow = new int[n][];
        parentCol = new int[n][];
        shift = new int[n][];
        size = new int[n][];
        index = new int[n][];
        nOrbits = n * n;
        zr = -1;
        zc = -1;
        for(int r = 0; r < n; r ++) {
            parentRow[r] = new int[n];
            parentCol[r] = new int[n];
            shift[r] = new int[n];
            size[r] = new int[n];
            index[r] = new int[n];
            for(int c = 0; c < n; c ++) {
                parentRow[r][c] = r;
                parentCol[r][c] = c;
                size[r][c] = 1;
                index[r][c] = -1;
            }
        }
    }

    public int phaseAdd(int x, int y) {
        int r = x + y;
        if (r >= order)
            r -= order;
        return r;
    }

    public int phaseSub(int x, int y) {
        int r = x - y;
        if (r < 0)
            r += order;
        return r;
    }

    public int phaseNegate(int p) {
        if (p == 0)
            return p;
        else
            return order - p;
    }

    /** Finds the representative of the cell to which (xr, xc) belongs
     * and returns the representative (yr, yc) and the multiplicative sign
     * such that M(xr, xc) = M(yr, yc) * rootOfUnity(yp, order)
     */
    public void find(int xr, int xc, int[] out) {
        int p = 0;
        while (true) {
            int yr = parentRow[xr][xc];
            int yc = parentCol[xr][xc];
            int yp = shift[xr][xc];
            if (xr == yr && xc == yc) {
                out[0] = yr;
                out[1] = yc;
                out[2] = p;
                return;
            }
            parentRow[xr][xc] = parentRow[yr][yc];
            parentCol[xr][xc] = parentCol[yr][yc];
            shift[xr][xc] = phaseAdd(yp, shift[yr][yc]);
            xr = yr;
            xc = yc;
            p = phaseAdd(p, yp);
        }
    }

    public int nZeros() {
        if (zr == -1 || zc == -1)
            return 0;
        else
            return size[zr][zc];
    }

    /** Compares indices (xr, xc) and (yr, yc) according to the lexicographic order */
    public int compareIndices(int xr, int xc, int yr, int yc) {
        if (xr < yr || (xr == yr && xc < yc))
            return -1;
        else if (xr > yr || (xr == yr && xc > yc))
            return 1;
        else
            return 0;
    }

    /** Merges the root (xr0, xc0) to (yr0, yc0)
     * so that M(xr0, xc0) = M(yr0, yc0) * rootOfUnity(p, order)
     * By convention, we require (xr0, xc0) > (yr0, yc0)
     */
     public void mergeTo(int xr0, int xc0, int yr0, int yc0, int p) {
         assert(xr0 != yr0 || xc0 != yc0);
         parentRow[xr0][xc0] = yr0;
         parentCol[xr0][xc0] = yc0;
         shift[xr0][xc0] = p;
         size[yr0][yc0] += size[xr0][xc0];
         size[xr0][xc0] = 0;
         nOrbits -= 1;
     }

    /** Merges the sets from which (xr, xc) and (yr, yc) are members
     * with a shift difference such that
     * M(xr, xc) = M(yr, yc) * rootOfUnity(p, order)
     */
    public void union(int xr, int xc, int yr, int yc, int p) {
        int[] x0 = new int[3];
        int[] y0 = new int[3];
        find(xr, xc, x0);
        find(yr, yc, y0);
        int xr0 = x0[0];
        int xc0 = x0[1];
        int xp0 = x0[2];
        int yr0 = y0[0];
        int yc0 = y0[1];
        int yp0 = y0[2];
        // shift difference between roots, i.e.
        // M(xr0, xc0) = M(yr0, yc0) * rootOfUnity(p0, order)
        int p0 = phaseSub(phaseAdd(p, yp0), xp0);

        // is our cell zero?
        if (xr0 == yr0 && xc0 == yc0 && p0 != 0) {
            // new zero element discovered
            if (zr == -1 || zc == -1) {
                // no zero element known until now
                zr = xr0;
                zc = xc0;
                shift[zr][zc] = 0;
                nOrbits --; // the zero orbit does not count
            } else {
                // zero element present
                shift[xr0][xc0] = 0;
                switch (compareIndices(zr, zc, xr0, xc0)) {
                    case -1:
                        // the zero representative stays the same
                        mergeTo(xr0, xc0, zr, zc, 0);
                        break;
                    case 1:
                        // the new zero cell has the zero representative
                        mergeTo(zr, zc, xr0, xc0, 0);
                        zr = xr0;
                        zc = xc0;
                        break;
                }
            }
        } else {
            switch (compareIndices(xr0, xc0, yr0, yc0)) {
                case -1:
                    mergeTo(yr0, yc0, xr0, xc0, phaseNegate(p0));
                    break;
                case 1:
                    mergeTo(xr0, xc0, yr0, yc0, p0);
                    break;
            }
        }
    }

}
