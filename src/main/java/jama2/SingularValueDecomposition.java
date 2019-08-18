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
     * <p>
     * This is a package-private constructor. Use {@link Matrix#svd()} to create a
     * cholesky decomposition of a given matrix.
     * </p>
     *
     * @param Arg Rectangular matrix
     * @see Matrix#svd()
     */
    SingularValueDecomposition(final Matrix Arg) {
        // Derived from LINPACK code.
        // Initialize.
        final var A = Arg.getArrayCopy();
        m = Arg.getRowDimension();
        n = Arg.getColumnDimension();
        final var nu = min(m, n);

        /*
         * Apparently the failing cases are only a proper subset of (m<n), so let's not
         * throw error. Correct fix to come later?
         */
        // if (m<n)
        // {
        // throw new IllegalArgumentException("Jama SVD only works for m >= n");
        // }

        s = new double[min(m + 1, n)];
        U = new double[m][nu];
        V = new double[n][n];
        final var e = new double[n];
        final var work = new double[m];
        final var wantu = true;
        final var wantv = true;

        // Reduce A to bidiagonal form, storing the diagonal elements
        // in s and the super-diagonal elements in e.
        final var nct = min(m - 1, n);
        final var nrt = max(0, min(n - 2, m));
        for (var k = 0; k < max(nct, nrt); k++) {
            if (k < nct) {
                // Compute 2-norm of k-th column without under/overflow.
                s[k] = 0.0;
                for (var i = k; i < m; i++) {
                    s[k] = hypot(s[k], A[i][k]);
                }

                // Compute the transformation for the k-th column
                if (s[k] != 0.0) {
                    if (A[k][k] < 0.0) {
                        s[k] = -s[k];
                    }
                    for (int i = k; i < m; i++) {
                        A[i][k] /= s[k];
                    }
                    A[k][k]++;
                }

                // place the k-th diagonal in s[k].
                s[k] = -s[k];
            }

            for (int j = k + 1; j < n; j++) {
                if (k < nct && s[k] != 0.0) {
                    // Apply the transformation.
                    var t = 0.0;
                    for (var i = k; i < m; i++) {
                        t += A[i][k] * A[i][j];
                    }
                    t /= -A[k][k];

                    for (var i = k; i < m; i++) {
                        A[i][j] += t * A[i][k];
                    }
                }

                // Place the k-th row of A into e for the
                // subsequent calculation of the row transformation.
                e[j] = A[k][j];
            }

            if (wantu && k < nct) {
                // Place the transformation in U for subsequent back multiplication.
                for (var i = k; i < m; i++) {
                    U[i][k] = A[i][k];
                }
            }

            if (k < nrt) {
                // Compute the k-th row transformation and place the
                // k-th super-diagonal in e[k].
                // Compute 2-norm without under/overflow.
                e[k] = 0.0;
                for (var i = k + 1; i < n; i++) {
                    e[k] = hypot(e[k], e[i]);
                }

                if (e[k] != 0.0) {
                    if (e[k + 1] < 0.0) {
                        e[k] = -e[k];
                    }
                    for (var i = k + 1; i < n; i++) {
                        e[i] /= e[k];
                    }
                    e[k + 1]++;
                }
                e[k] = -e[k];

                if (k + 1 < m && e[k] != 0.0) {
                    // Apply the transformation.
                    for (var i = k + 1; i < m; i++) {
                        work[i] = 0.0;
                    }
                    for (var j = k + 1; j < n; j++) {
                        for (var i = k + 1; i < m; i++) {
                            work[i] += e[j] * A[i][j];
                        }
                    }
                    for (var j = k + 1; j < n; j++) {
                        final var t = -e[j] / e[k + 1];
                        for (var i = k + 1; i < m; i++) {
                            A[i][j] += t * work[i];
                        }
                    }
                }

                if (wantv) {
                    // Place the transformation in V for subsequent back multiplication.
                    for (var i = k + 1; i < n; i++) {
                        V[i][k] = e[i];
                    }
                }
            }
        }

        // Set up the final bidiagonal matrix or order p.
        var p = min(n, m + 1);
        if (nct < n) {
            s[nct] = A[nct][nct];
        }

        if (m < p) {
            s[p - 1] = 0.0;
        }

        if (nrt + 1 < p) {
            e[nrt] = A[nrt][p - 1];
        }

        e[p - 1] = 0.0;

        // If required, generate U.
        if (wantu) {
            for (var j = nct; j < nu; j++) {
                for (final var rowU : U) {
                    rowU[j] = 0.0;
                }
                U[j][j] = 1.0;
            }

            for (int k = nct - 1; k >= 0; k--) {
                if (s[k] != 0.0) {
                    for (int j = k + 1; j < nu; j++) {
                        var t = 0.0;
                        for (var i = k; i < m; i++) {
                            t += U[i][k] * U[i][j];
                        }
                        t /= -U[k][k];

                        for (var i = k; i < m; i++) {
                            U[i][j] += t * U[i][k];
                        }
                    }
                    for (var i = k; i < m; i++) {
                        U[i][k] = -U[i][k];
                    }
                    U[k][k]++;

                    for (var i = 0; i < k - 1; i++) {
                        U[i][k] = 0.0;
                    }
                } else {
                    for (final var rowU : U) {
                        rowU[k] = 0.0;
                    }
                    U[k][k] = 1.0;
                }
            }
        }

        // If required, generate V.
        if (wantv) {
            for (int k = n - 1; k >= 0; k--) {
                if (k < nrt && e[k] != 0.0) {
                    for (int j = k + 1; j < nu; j++) {
                        var t = 0.0;
                        for (var i = k + 1; i < n; i++) {
                            t += V[i][k] * V[i][j];
                        }
                        t /= -V[k + 1][k];

                        for (var i = k + 1; i < n; i++) {
                            V[i][j] += t * V[i][k];
                        }
                    }
                }
                for (final double[] rowV : V) {
                    rowV[k] = 0.0;
                }
                V[k][k] = 1.0;
            }
        }

        // Main iteration loop for the singular values.
        final var pp = p - 1;
        while (p > 0) {
            var k = p - 2;
            while (k >= -1) {
                if (k == -1) {
                    break;
                }
                if (abs(e[k]) <= tiny + eps * (abs(s[k]) + abs(s[k + 1]))) {
                    e[k] = 0.0;
                    break;
                }
                k--;
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
            if (k == p - 2) {
                // e(p-1) is negligible (convergence).
                kase = 4;
            } else {
                var ks = p - 1;
                while (ks >= k) {
                    if (ks == k) {
                        break;
                    }
                    
                    final var t = (ks != p ? abs(e[ks]) : 0.0) +
                            (ks != k + 1 ? abs(e[ks - 1]) : 0.0);
                    if (abs(s[ks]) <= tiny + eps * t) {
                        s[ks] = 0.0;
                        break;
                    }
                    
                    ks--;
                }

                if (ks == k) {
                    // e[k-1] is negligible, k<p, and s(k), ..., s(p) are not
                    // negligible (qr step).
                    kase = 3;
                } else if (ks == p - 1) {
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
                var f = e[p - 2];
                e[p - 2] = 0.0;

                for (var j = p - 2; j >= k; j--) {
                    final var h = hypot(s[j], f);
                    final var cs = s[j] / h;
                    final var sn = f / h;
                    s[j] = h;
                    if (j != k) {
                        f = -sn * e[j - 1];
                        e[j - 1] *= cs;
                    }

                    if (wantv) {
                        for (var i = 0; i < n; i++) {
                            final var t = cs*V[i][j] + sn*V[i][p-1];
                            V[i][p-1] = -sn*V[i][j] + cs*V[i][p-1];
                            V[i][j] = t;
                         }
                    }
                }
            }
                break;

            // Split at negligible s(k).
            case 2: {
                // Remember e[k - 1] in f and reset it
                double f = e[k - 1];
                e[k - 1] = 0.0;

                for (int j = k; j < p; j++) {
                    var t = hypot(s[j], f);
                    final double cs = s[j] / t, sn = f / t;
                    s[j] = t;
                    f = -sn * e[j];
                    e[j] *= cs;

                    if (wantu) {
                        for (final var rowU : U) {
                            // Update U[i][k - 1] and U[i][j]
                            t = rowU[k - 1]; // remember U[i][k - 1]
                            rowU[k - 1] = cs * t - sn * rowU[j];
                            rowU[j] = cs * rowU[j] + sn * t;
                        }
                    }
                }
            }
                break;

            // Perform one qr step.
            case 3: {
                final var scale = max(
                        max(max(max(abs(s[p - 1]), abs(s[p - 2])), abs(e[p - 2])), abs(s[k])),
                        abs(e[k]));
                final var sp = s[p - 1] / scale;
                final var spm1 = s[p - 2] / scale;
                final var epm1 = e[p - 2] / scale;
                final var sk = s[k] / scale;
                final var ek = e[k] / scale;
                final var b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
                final var c = (sp * epm1) * (sp * epm1);

                // Calculate the shift.
                double shift;
                if (b != 0.0 || c != 0.0) {
                    shift = sqrt(b * b + c);
                    if (b < 0.0) {
                        shift = -shift;
                    }
                    shift = c / (b + shift);
                } else {
                    shift = 0.0;
                }

                // Chase zeros.
                var f = (sk + sp) * (sk - sp) + shift;
                var g = sk * ek;
                for (var j = k; j < p - 1; j++) {
                    var h = hypot(f, g);
                    var cs = f / h;
                    var sn = g / h;
                    if (j != k) {
                        e[j - 1] = h;
                    }
                    f = cs * s[j] + sn * e[j];
                    e[j] = cs * e[j] - sn * s[j];
                    g = sn * s[j + 1];
                    s[j + 1] *= cs;

                    if (wantv) {
                        for (var i = 0; i < n; i++) {
                            final var t = cs*V[i][j] + sn*V[i][j+1];
                            V[i][j+1] = -sn*V[i][j] + cs*V[i][j+1];
                            V[i][j] = t;
                        }
                    }

                    h = hypot(f, g);
                    cs = f / h;
                    sn = g / h;
                    s[j] = h;
                    f = cs * e[j] + sn * s[j + 1];
                    s[j + 1] = cs * s[j + 1] - sn * e[j];
                    g = sn * e[j + 1];
                    e[j + 1] *= cs;

                    if (wantu && j < m - 1) {
                        for (var  i = 0; i < m; i++) {
                            final var t = cs*U[i][j] + sn*U[i][j+1];
                            U[i][j+1] = -sn*U[i][j] + cs*U[i][j+1];
                            U[i][j] = t;
                         }
                    }
                }
                e[p - 2] = f;
            }
                break;

            // Convergence.
            case 4: {
                // Make the singular values positive.
                if (s[k] <= 0.0) {
                    s[k] = -s[k];
                    if (wantv) {
                        for (var i = 0; i <= pp; i++) {
                            V[i][k] = -V[i][k];
                        }
                    }
                }

                // Order the singular values.
                while (k < pp) {
                    if (s[k] < s[k + 1]) {
                        // swap s[k] and s[k + 1]
                        double t = s[k];
                        s[k] = s[k + 1];
                        s[k + 1] = t;

                        if (wantv && k < n - 1) {
                            for (final double[] rowV : V) {
                                // swap rowV[k + 1] and rowV[k]
                                t = rowV[k + 1];
                                rowV[k + 1] = rowV[k];
                                rowV[k] = t;
                            }
                        }

                        if (wantu && k < m - 1) {
                            for (final double[] rowU : U) {
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
        return s[0] / s[min(m, n) - 1];
    }

    /**
     * Return the diagonal matrix of singular values.
     *
     * @return S
     */
    public Matrix getS() {
        return Matrix.diag(s);
    }

    /**
     * Return the one-dimensional array of singular values.
     *
     * @return diagonal of S.
     */
    public double[] getSingularValues() {
        return s;
    }

    /**
     * Return the left singular vectors.
     *
     * @return U
     */
    public Matrix getU() {
        return new Matrix(m, min(m, n), U);
    }

    /**
     * Return the right singular vectors.
     *
     * @return V
     */
    public Matrix getV() {
        return new Matrix(n, V);
    }

    /**
     * Two norm.
     *
     * @return max(S)
     */
    public double norm2() {
        return s[0];
    }

    /**
     * Effective numerical matrix rank.
     *
     * @return Number of nonnegligible singular values.
     */
    public int rank() {
        final double tol = max(m, n) * s[0] * eps;
        int r = 0;
        for (final double element : s) {
            if (element > tol) {
                r++;
            }
        }
        return r;
    }
}
