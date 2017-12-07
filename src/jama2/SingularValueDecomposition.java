package jama2;

import java.io.Serializable;

import static jama2.util.Maths.eps;
import static jama2.util.Maths.hypot;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

/**
 * Singular Value Decomposition.
 * <P>
 * For an m-by-n matrix A with m &gt;= n, the singular value decomposition is an
 * m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and an n-by-n
 * orthogonal matrix V so that A = U*S*V'.
 * </P>
 * <P>
 * The singular values, sigma[k] = S[k][k], are ordered so that sigma[0] &gt;=
 * sigma[1] &gt;= ... &gt;= sigma[n-1].
 * </P>
 * <P>
 * The singular value decompostion always exists, so the constructor will never
 * fail. The matrix condition number and the effective numerical rank can be
 * computed from this decomposition.
 * </P>
 *
 * @author The MathWorks, Inc. and the National Institute of Standards and
 *         Technology.
 * @version 2.0
 * @see <a href="http://tweimer.github.io/java-matrix/">java-matrix</a>
 */
public class SingularValueDecomposition implements Serializable {
    /**
     * For the Serializeable interface
     */
    private static final long serialVersionUID = 1;

    final static double tiny = Math.pow(2.0, -966.0);

    /**
     * Arrays for internal storage of U and V.
     *
     * @serial internal storage of U.
     * @serial internal storage of V.
     */
    private final double[][] U, V;

    /**
     * Array for internal storage of singular values.
     *
     * @serial internal storage of singular values.
     */
    private final double[] s;

    /**
     * Row and column dimensions.
     *
     * @serial row dimension.
     * @serial column dimension.
     */
    private final int m, n;

    /**
     * Construct the singular value decomposition Structure to access U, S and V.
     * 
     * <p>This is a package-private constructor.
     * Use {@link Matrix#svd()} to create a cholesky decomposition of a given matrix.</p>
     *
     * @param Arg
     *            Rectangular matrix
     * @see Matrix#svd()
     */
    SingularValueDecomposition(final Matrix Arg) {
        // Derived from LINPACK code.
        // Initialize.
        final double[][] A = Arg.getArrayCopy();
        this.m = Arg.getRowDimension();
        this.n = Arg.getColumnDimension();
        final int nu = min(this.m, this.n);

        /*
         * Apparently the failing cases are only a proper subset of (m<n), so
         * let's not throw error. Correct fix to come later?
         */
        // if (m<n)
        // {
        // throw new IllegalArgumentException("Jama SVD only works for m >= n");
        // }

        this.s = new double[min(this.m + 1, this.n)];
        this.U = new double[this.m][nu];
        this.V = new double[this.n][this.n];
        final double[] e = new double[this.n], work = new double[this.m];
        final boolean wantu = true, wantv = true;

        // Reduce A to bidiagonal form, storing the diagonal elements
        // in s and the super-diagonal elements in e.
        final int nct = min(this.m - 1, this.n),
                nrt = max(0, min(this.n - 2, this.m));
        for (int k = 0; k < max(nct, nrt); k++) {
            if (k < nct) {
                // Compute 2-norm of k-th column without under/overflow.
                this.s[k] = 0;
                for (int i = k; i < this.m; i++) {
                    this.s[k] = hypot(this.s[k], A[i][k]);
                }

                // Compute the transformation for the k-th column
                if (this.s[k] != 0D) {
                    if (A[k][k] < 0D) {
                        this.s[k] = -this.s[k];
                    }
                    for (int i = k; i < this.m; i++) {
                        A[i][k] /= this.s[k];
                    }
                    A[k][k]++;
                }

                // place the k-th diagonal in s[k].
                this.s[k] = -this.s[k];
            }

            for (int j = k + 1; j < this.n; j++) {
                if ((k < nct) && (this.s[k] != 0D)) {
                    // Apply the transformation.
                    double t = 0D;
                    for (int i = k; i < this.m; i++) {
                        t += A[i][k] * A[i][j];
                    }
                    t /= -A[k][k];

                    for (int i = k; i < this.m; i++) {
                        A[i][j] += t * A[i][k];
                    }
                }

                // Place the k-th row of A into e for the
                // subsequent calculation of the row transformation.
                e[j] = A[k][j];
            }

            if (wantu && (k < nct)) {
                // Place the transformation in U for subsequent back
                // multiplication.
                for (int i = k; i < this.m; i++) {
                    this.U[i][k] = A[i][k];
                }
            }

            if (k < nrt) {
                // Compute the k-th row transformation and place the
                // k-th super-diagonal in e[k].
                // Compute 2-norm without under/overflow.
                e[k] = 0D;
                for (int i = k + 1; i < this.n; i++) {
                    e[k] = hypot(e[k], e[i]);
                }

                if (e[k] != 0D) {
                    if (e[k + 1] < 0D) {
                        e[k] = -e[k];
                    }
                    for (int i = k + 1; i < this.n; i++) {
                        e[i] /= e[k];
                    }
                    e[k + 1]++;
                }
                e[k] = -e[k];

                if (((k + 1) < this.m) && (e[k] != 0.0)) {
                    // Apply the transformation.
                    for (int i = k + 1; i < this.m; i++) {
                        work[i] = 0D;
                    }
                    for (int j = k + 1; j < this.n; j++) {
                        for (int i = k + 1; i < this.m; i++) {
                            work[i] += (e[j] * A[i][j]);
                        }
                    }
                    for (int j = k + 1; j < this.n; j++) {
                        final double t = -e[j] / e[k + 1];
                        for (int i = k + 1; i < this.m; i++) {
                            A[i][j] += (t * work[i]);
                        }
                    }
                }

                if (wantv) {
                    // Place the transformation in V for subsequent back
                    // multiplication.
                    for (int i = k + 1; i < this.n; i++) {
                        this.V[i][k] = e[i];
                    }
                }
            }
        }

        // Set up the final bidiagonal matrix or order p.
        int p = min(this.n, this.m + 1);
        if (nct < this.n) {
            this.s[nct] = A[nct][nct];
        }

        if (this.m < p) {
            this.s[p - 1] = 0D;
        }

        if ((nrt + 1) < p) {
            e[nrt] = A[nrt][p - 1];
        }

        e[p - 1] = 0D;

        // If required, generate U.
        if (wantu) {
            for (int j = nct; j < nu; j++) {
                for (final double[] rowU : this.U) {
                    rowU[j] = 0D;
                }
                this.U[j][j] = 1D;
            }

            for (int k = nct - 1; k >= 0; k--) {
                if (this.s[k] != 0D) {
                    for (int j = k + 1; j < nu; j++) {
                        double t = 0D;
                        for (int i = k; i < this.m; i++) {
                            t += this.U[i][k] * this.U[i][j];
                        }
                        t /= -this.U[k][k];

                        for (int i = k; i < this.m; i++) {
                            this.U[i][j] += t * this.U[i][k];
                        }
                    }
                    for (int i = k; i < this.m; i++) {
                        this.U[i][k] = -this.U[i][k];
                    }
                    this.U[k][k]++;

                    for (int i = 0; i < (k - 1); i++) {
                        this.U[i][k] = 0D;
                    }
                } else {
                    for (final double[] rowU : this.U) {
                        rowU[k] = 0D;
                    }
                    this.U[k][k] = 1D;
                }
            }
        }

        // If required, generate V.
        if (wantv) {
            for (int k = this.n - 1; k >= 0; k--) {
                if ((k < nrt) && (e[k] != 0D)) {
                    for (int j = k + 1; j < nu; j++) {
                        double t = 0;
                        for (int i = k + 1; i < this.n; i++) {
                            t += this.V[i][k] * this.V[i][j];
                        }
                        t /= -this.V[k + 1][k];

                        for (int i = k + 1; i < this.n; i++) {
                            this.V[i][j] += t * this.V[i][k];
                        }
                    }
                }
                for (final double[] rowV : this.V) {
                    rowV[k] = 0D;
                }
                this.V[k][k] = 1D;
            }
        }

        // Main iteration loop for the singular values.
        final int pp = p - 1;
        while (p > 0) {
            int k = p - 2;
            for (; k >= -1; k--) {
                if (k == -1) {
                    break;
                } else if (abs(e[k]) <= (tiny + (eps * (abs(this.s[k]) + abs(this.s[k + 1]))))) {
                    e[k] = 0D;
                    break;
                }
            }

            // Here is where a test for too many iterations would go.
            //
            // This section of the program inspects for negligible elements in
            // the s and e arrays.
            //
            // On completion the variables kase and k are set as follows.
            //
            // kase = 1 if s(p) and e[k-1] are negligible and k<p
            // kase = 2 if s(k) is negligible and k<p
            // kase = 3 if e[k-1] is negligible, k<p, and s(k), ..., s(p) are
            // not negligible (qr step).
            // kase = 4 if e(p-1) is negligible (convergence).
            final byte kase;
            if (k == (p - 2)) {
                // e(p-1) is negligible (convergence).
                kase = 4;
            } else {
                int ks = p - 1;
                while (ks > k) {
                    final double t = (ks != p ? abs(e[ks]) : 0D) + (ks != (k + 1) ? abs(e[ks - 1]) : 0D);
                    if (abs(this.s[ks]) <= (tiny + (eps * t))) {
                        this.s[ks] = 0D;
                        break;
                    } else {
                        ks--;
                    }
                }

                if (ks == k) {
                    // e[k-1] is negligible, k<p, and s(k), ..., s(p) are not
                    // negligible (qr step).
                    kase = 3;
                } else if (ks == (p - 1)) {
                    // s(p) and e[k-1] are negligible and k<p
                    kase = 1;
                } else {
                    // s(k) is negligible and k<p
                    kase = 2;
                    k = ks;
                }
            }
            k++;

            // Perform the task indicated by kase.
            switch (kase) {
            // Deflate negligible s(p).
            case 1: {
                // Remember e[p - 2] in f and reset it
                double f = e[p - 2];
                e[p - 2] = 0D;

                for (int j = p - 2; j >= k; j--) {
                    double t = hypot(this.s[j], f);
                    final double cs = (this.s[j] / t), sn = (f / t);
                    this.s[j] = t;
                    if (j != k) {
                        f = -sn * e[j - 1];
                        e[j - 1] *= cs;
                    }

                    if (wantv) {
                        for (final double rowV[] : this.V) {
                            // remember V[i][p - 1]
                            t = rowV[p - 1];
                            rowV[p - 1] = (cs * t) - (sn * rowV[j]);
                            rowV[j] = (cs * rowV[j]) + (sn * t);
                        }
                    }
                }
            }
                break;

            // Split at negligible s(k).
            case 2: {
                // Remember e[k - 1] in f and reset it
                double f = e[k - 1];
                e[k - 1] = 0D;

                for (int j = k; j < p; j++) {
                    double t = hypot(this.s[j], f);
                    final double cs = (this.s[j] / t), sn = (f / t);
                    this.s[j] = t;
                    f = -sn * e[j];
                    e[j] *= cs;

                    if (wantu) {
                        for (final double[] rowU : this.U) {
                            // Update U[i][k - 1] and U[i][j]
                            t = rowU[k - 1]; // remember U[i][k - 1]
                            rowU[k - 1] = (cs * t) - (sn * rowU[j]);
                            rowU[j] = (cs * rowU[j]) + (sn * t);
                        }
                    }
                }
            }
                break;

            // Perform one qr step.
            case 3: {
                final double scale = max(
                        max(max(max(abs(this.s[p - 1]),
                                abs(this.s[p - 2])), abs(e[p - 2])), abs(this.s[k])),
                                    abs(e[k]));
                final double sp = this.s[p - 1] / scale;
                final double spm1 = this.s[p - 2] / scale,
                        epm1 = e[p - 2] / scale;
                final double sk = this.s[k] / scale, ek = e[k] / scale;
                final double b = (((spm1 + sp) * (spm1 - sp)) + (epm1 * epm1)) / 2.0,
                        c = (sp * epm1) * (sp * epm1);

                // Calculate the shift.
                double shift;
                if ((b != 0D) || (c != 0D)) {
                    shift = sqrt((b * b) + c);
                    if (b < 0D) {
                        shift = -shift;
                    }
                    shift = c / (b + shift);
                } else {
                    shift = 0D;
                }

                // Chase zeros.
                double f = (((sk + sp) * (sk - sp)) + shift), g = (sk * ek);
                for (int j = k; j < (p - 1); j++) {
                    double t = hypot(f, g), cs = f / t, sn = g / t;
                    if (j != k) {
                        e[j - 1] = t;
                    }
                    f = (cs * this.s[j]) + (sn * e[j]);
                    e[j] = (cs * e[j]) - (sn * this.s[j]);
                    g = sn * this.s[j + 1];
                    this.s[j + 1] *= cs;

                    if (wantv) {
                        for (final double[] rowV : this.V) {
                            t = rowV[j + 1]; // remember V[i][j + 1]
                            rowV[j + 1] = (cs * t) - (sn * rowV[j]);
                            rowV[j] = (cs * rowV[j]) + (sn * t);
                        }
                    }

                    t = hypot(f, g);
                    cs = f / t;
                    sn = g / t;
                    this.s[j] = t;
                    f = ((cs * e[j]) + (sn * this.s[j + 1]));
                    this.s[j + 1] = ((cs * this.s[j + 1]) - (sn * e[j]));
                    g = sn * e[j + 1];
                    e[j + 1] *= cs;

                    if (wantu && (j < (this.m - 1))) {
                        for (final double[] rowU : this.U) {
                            t = rowU[j + 1];
                            rowU[j + 1] = ((cs * t) - (sn * rowU[j]));
                            rowU[j] = ((cs * rowU[j]) + (sn * t));
                        }
                    }
                }
                e[p - 2] = f;
            }
                break;

            // Convergence.
            case 4: {
                // Make the singular values positive.
                if (this.s[k] <= 0D) {
                    this.s[k] = -this.s[k];
                    if (wantv) {
                        for (int i = 0; i <= pp; i++) {
                            this.V[i][k] = -this.V[i][k];
                        }
                    }
                }

                // Order the singular values.
                while (k < pp) {
                    if (this.s[k] < this.s[k + 1]) {
                        // swap s[k] and s[k + 1]
                        double t = this.s[k];
                        this.s[k] = this.s[k + 1];
                        this.s[k + 1] = t;

                        if (wantv && (k < (this.n - 1))) {
                            for (final double[] rowV : this.V) {
                                // swap rowV[k + 1] and rowV[k]
                                t = rowV[k + 1];
                                rowV[k + 1] = rowV[k];
                                rowV[k] = t;
                            }
                        }

                        if (wantu && (k < (this.m - 1))) {
                            for (final double[] rowU : this.U) {
                                // swap rowU[k + 1] and rowU[k]
                                t = rowU[k + 1];
                                rowU[k + 1] = rowU[k];
                                rowU[k] = t;
                            }
                        }
                        k++;
                    } else {
                        break;
                    }
                }
                p--;
            }
                break;
            }
        }
    }

    /**
     * Two norm condition number.
     *
     * @return max(S)/min(S)
     */
    public double cond() {
        return this.s[0] / this.s[min(this.m, this.n) - 1];
    }

    /**
     * Return the diagonal matrix of singular values.
     *
     * @return S
     */
    public Matrix getS() {
        final Matrix X = new Matrix(this.n);
        final double[][] S = X.getArray();
        for (int i = 0; i < this.n; i++) {
            S[i][i] = this.s[i];
        }
        return X;
    }

    /**
     * Return the one-dimensional array of singular values.
     *
     * @return diagonal of S.
     */
    public double[] getSingularValues() {
        return this.s;
    }

    /**
     * Return the left singular vectors.
     *
     * @return U
     */
    public Matrix getU() {
        return new Matrix(this.m, min(this.m + 1, this.n), this.U);
    }

    /**
     * Return the right singular vectors.
     *
     * @return V
     */
    public Matrix getV() {
        return new Matrix(this.n, this.V);
    }

    /**
     * Two norm.
     *
     * @return max(S)
     */
    public double norm2() {
        return this.s[0];
    }

    /**
     * Effective numerical matrix rank.
     *
     * @return Number of nonnegligible singular values.
     */
    public int rank() {
        final double tol = max(this.m, this.n) * this.s[0] * eps;
        int r = 0;
        for (final double element : this.s) {
            if (element > tol) {
                r++;
            }
        }
        return r;
    }
}
