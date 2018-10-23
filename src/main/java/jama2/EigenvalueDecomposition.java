package jama2;

import java.io.*;

import jama2.util.Maths;

/**
 * Eigenvalues and eigenvectors of a real matrix.
 * <P>
 * If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is diagonal
 * and the eigenvector matrix V is orthogonal. I.e. A =
 * V.times(D.times(V.transpose())) and V.times(V.transpose()) equals the
 * identity matrix.
 * </P>
 * <P>
 * If A is not symmetric, then the eigenvalue matrix D is block diagonal with
 * the real eigenvalues in 1-by-1 blocks and any complex eigenvalues, lambda +
 * i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda]. The columns of V represent
 * the eigenvectors in the sense that A*V = V*D, i.e. A.times(V) equals
 * V.times(D). The matrix V may be badly conditioned, or even singular, so the
 * validity of the equation A = V*D*inverse(V) depends upon V.cond().
 * </P>
 * 
 * @author The MathWorks, Inc. and the National Institute of Standards and
 *         Technology.
 * @version 2.0
 * @see <a href="http://tweimer.github.io/java-matrix/">java-matrix</a>
 **/
public class EigenvalueDecomposition implements Serializable {

    private static final long serialVersionUID = 1;

    /**
     * Row and column dimension (square matrix).
     * 
     * @serial matrix dimension.
     */
    private final int n;

    /**
     * Symmetry flag.
     * 
     * @serial internal symmetry flag.
     */
    private boolean issymmetric;

    /**
     * Arrays for internal storage of eigenvalues.
     * 
     * @serial internal storage of eigenvalues.
     */
    private final double[] d, e;

    /**
     * Array for internal storage of eigenvectors (n-times-n).
     * 
     * @serial internal storage of eigenvectors.
     */
    private final double[][] V;

    /**
     * Array for internal storage of nonsymmetric Hessenberg form.
     * 
     * @serial internal storage of nonsymmetric Hessenberg form.
     */
    private double[][] H;

    /**
     * Working storage for nonsymmetric algorithm.
     * 
     * @serial working storage for nonsymmetric algorithm.
     */
    private double[] ort;

    /**
     * Real and imaginary result of a complex division.
     */
    private transient double cdivr, cdivi;

    /**
     * Check for symmetry, then construct the eigenvalue decomposition Structure
     * to access D and V.
     * 
     * <p>This is a package-private constructor.
     * Use {@link Matrix#eig()} to create a cholesky decomposition of a given matrix.</p>
     * 
     * @param Arg
     *            Square matrix
     * @see Matrix#eig()
     */
    EigenvalueDecomposition(final Matrix Arg) {
        final double[][] A = Arg.getArray();
        this.n = Arg.getColumnDimension();
        this.V = new double[this.n][this.n];
        this.d = new double[this.n];
        this.e = new double[this.n];

        this.issymmetric = true;
        for (int j = 0; (j < this.n) && this.issymmetric; j++) {
            for (int i = 0; i < this.n; i++) {
                if (A[i][j] != A[j][i]) {
                    this.issymmetric = false;
                    break;
                }
            }
        }

        if (this.issymmetric) {
            for (int i = 0; i < this.n; i++) {
                for (int j = 0; j < this.n; j++) {
                    this.V[i][j] = A[i][j];
                }
            }

            // Tridiagonalize.
            this.tred2();

            // Diagonalize.
            this.tql2();
        } else {
            this.H = new double[this.n][this.n];
            for (int j = 0; j < this.n; j++) {
                for (int i = 0; i < this.n; i++) {
                    this.H[i][j] = A[i][j];
                }
            }

            // Reduce to Hessenberg form.
            this.orthes();

            // Reduce Hessenberg to real Schur form.
            this.hqr2();
        }
    }

    /**
     * Complex scalar division. Result is stored in cdivr+cdivi*i
     * 
     * @param xr
     *            x (real part)
     * @param xi
     *            x (imaginary part)
     * @param yr
     *            y (real part)
     * @param yi
     *            y (imaginary part)
     */
    private void cdiv(final double xr, final double xi, final double yr, final double yi) {
        final double r, d1;
        if (Math.abs(yr) > Math.abs(yi)) {
            r = yi / yr;
            d1 = yr + (r * yi);
            this.cdivr = (xr + (r * xi)) / d1;
            this.cdivi = (xi - (r * xr)) / d1;
        } else {
            r = yr / yi;
            d1 = yi + (r * yr);
            this.cdivr = ((r * xr) + xi) / d1;
            this.cdivi = ((r * xi) - xr) / d1;
        }
    }

    /**
     * Return the block diagonal eigenvalue matrix
     * 
     * @return D
     */
    public Matrix getD() {
        final Matrix X = new Matrix(this.n, this.n);
        final double[][] D = X.getArray();
        for (int i = 0; i < this.n; i++) {
            D[i][i] = this.d[i];
            if (this.e[i] > 0) {
                D[i][i + 1] = this.e[i];
            } else if (this.e[i] < 0) {
                D[i][i - 1] = this.e[i];
            }
        }
        return X;
    }

    /**
     * Return the imaginary parts of the eigenvalues
     * 
     * @return imag(diag(D))
     */
    public double[] getImagEigenvalues() {
        return this.e;
    }

    /**
     * Return the real parts of the eigenvalues
     * 
     * @return real(diag(D))
     */
    public double[] getRealEigenvalues() {
        return this.d;
    }

    /**
     * Return the eigenvector matrix
     * 
     * @return V
     */
    public Matrix getV() {
        return new Matrix(this.n, this.n, this.V);
    }

    /**
     * Nonsymmetric reduction from Hessenberg to real Schur form.
     * 
     * This is derived from the Algol procedure hqr2, by Martin and Wilkinson,
     * Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the correspondin
     * Fortran subroutine in EISPACK.
     */
    private void hqr2() {
        // Initialize
        final int nn = this.n;
        int n1 = nn - 1;
        final int low = 0, high = nn - 1;
        double exshift = 0.0, p = 0, q = 0, r = 0, s = 0, z = 0;

        // Store roots isolated by balanc and compute matrix norm
        double norm = 0.0;
        for (int i = 0; i < nn; i++) {
            if ((i < low) || (i > high)) {
                this.d[i] = this.H[i][i];
                this.e[i] = 0.0;
            }
            for (int j = Math.max(i - 1, 0); j < nn; j++) {
                norm += Math.abs(this.H[i][j]);
            }
        }

        // Outer loop over eigenvalue index
        int iter = 0;
        while (n1 >= low) {
            // Look for single small sub-diagonal element
            int l = n1;
            while (l > low) {
                s = Math.abs(this.H[l - 1][l - 1]) + Math.abs(this.H[l][l]);
                if (s == 0D) {
                    s = norm;
                }
                if (Math.abs(this.H[l][l - 1]) < (Maths.eps * s)) {
                    break;
                } else {
                    l--;
                }
            }

            // Check for convergence

            // One root found
            if (l == n1) {
                this.d[n1] = (this.H[n1][n1] += exshift);
                this.e[n1--] = 0D;
                iter = 0;
            }
            // Two roots found
            else if (l == (n1 - 1)) {
                var w = this.H[n1][n1 - 1] * this.H[n1 - 1][n1];
                p = (this.H[n1 - 1][n1 - 1] - this.H[n1][n1]) / 2.0;
                q = (p * p) + w;
                z = Math.sqrt(Math.abs(q));
                this.H[n1 - 1][n1 - 1] += exshift;
                var x = (this.H[n1][n1] += exshift);

                // Real pair
                if (q >= 0) {
                    z = ((p >= 0) ? (p + z) : (p - z));
                    this.d[n1] = this.d[n1 - 1] = x + z;
                    if (z != 0.0) {
                        this.d[n1] = x - (w / z);
                    }
                    this.e[n1 - 1] = this.e[n1] = 0D;
                    x = this.H[n1][n1 - 1];
                    s = Math.abs(x) + Math.abs(z);
                    p = x / s;
                    q = z / s;
                    r = Math.sqrt((p * p) + (q * q));
                    p /= r;
                    q /= r;

                    // Row modification
                    for (int j = n1 - 1; j < nn; j++) {
                        z = this.H[n1 - 1][j];
                        this.H[n1 - 1][j] = (q * z) + (p * this.H[n1][j]);
                        this.H[n1][j] = (q * this.H[n1][j]) - (p * z);
                    }

                    // Column modification
                    for (int i = 0; i <= n1; i++) {
                        z = this.H[i][n1 - 1];
                        this.H[i][n1 - 1] = (q * z) + (p * this.H[i][n1]);
                        this.H[i][n1] = (q * this.H[i][n1]) - (p * z);
                    }

                    // Accumulate transformations
                    for (int i = low; i <= high; i++) {
                        z = this.V[i][n1 - 1];
                        this.V[i][n1 - 1] = (q * z) + (p * this.V[i][n1]);
                        this.V[i][n1] = (q * this.V[i][n1]) - (p * z);
                    }
                }
                // Complex pair
                else {
                    this.d[n1 - 1] = this.d[n1] = x + p;
                    this.e[n1 - 1] = z;
                    this.e[n1] = -z;
                }
                n1 -= 2;
                iter = 0;
                // No convergence yet
            }
            // Form shift
            else {
                var x = this.H[n1][n1];
                var y = 0D;
                var w = 0D;
                if (l < n1) {
                    y = this.H[n1 - 1][n1 - 1];
                    w = this.H[n1][n1 - 1] * this.H[n1 - 1][n1];
                }

                // Wilkinson's original ad hoc shift
                if (iter == 10) {
                    exshift += x;
                    for (int i = low; i <= n1; i++) {
                        this.H[i][i] -= x;
                    }
                    s = Math.abs(this.H[n1][n1 - 1]) + Math.abs(this.H[n1 - 1][n1 - 2]);
                    x = y = 0.75 * s;
                    w = -0.4375 * s * s;
                }

                // MATLAB's new ad hoc shift
                if (iter == 30) {
                    s = (y - x) / 2.0;
                    s = (s * s) + w;
                    if (s > 0) {
                        s = Math.sqrt(s);
                        if (y < x) {
                            s = -s;
                        }
                        s = x - (w / (((y - x) / 2.0) + s));
                        for (int i = low; i <= n1; i++) {
                            this.H[i][i] -= s;
                        }
                        exshift += s;
                        x = y = w = 0.964;
                    }
                }

                iter++; // (Could check iteration count here.)

                // Look for two consecutive small sub-diagonal elements
                int m = n1 - 2;
                while (m >= l) {
                    z = this.H[m][m];
                    r = x - z;
                    s = y - z;
                    p = (((r * s) - w) / this.H[m + 1][m]) + this.H[m][m + 1];
                    q = this.H[m + 1][m + 1] - z - r - s;
                    r = this.H[m + 2][m + 1];
                    s = Maths.norm1(p, q, r);
                    p /= s;
                    q /= s;
                    r /= s;
                    if (m == l) {
                        break;
                    }
                    if ((Math.abs(this.H[m][m - 1]) * (Math.abs(q) + Math.abs(r))) < (Maths.eps
                            * (Math.abs(p) * (Math.abs(this.H[m - 1][m - 1]) + Math.abs(z) + Math.abs(this.H[m + 1][m + 1]))))) {
                        break;
                    }
                    m--;
                }

                for (int i = m + 2; i <= n1; i++) {
                    this.H[i][i - 2] = 0D;
                    if (i > (m + 2)) {
                        this.H[i][i - 3] = 0D;
                    }
                }

                // Double QR step involving rows l:n and columns m:n
                for (int k = m; k <= (n1 - 1); k++) {
                    final boolean notlast = (k < (n1 - 1));
                    if (k != m) {
                        p = this.H[k][k - 1];
                        q = this.H[k + 1][k - 1];
                        r = (notlast ? this.H[k + 2][k - 1] : 0.0);
                        x = Maths.norm1(p, q, r);
                        if (x == 0D) {
                            continue;
                        } else {
                            p /= x;
                            q /= x;
                            r /= x;
                        }
                    }

                    s = Maths.hypot(p, q, r);
                    if (p < 0) {
                        s = -s;
                    }
                    if (s != 0) {
                        if (k != m) {
                            this.H[k][k - 1] = -s * x;
                        } else if (l != m) {
                            this.H[k][k - 1] = -this.H[k][k - 1];
                        }
                        p += s;
                        x = p / s;
                        y = q / s;
                        z = r / s;
                        q /= p;
                        r /= p;

                        // Row modification
                        for (int j = k; j < nn; j++) {
                            p = this.H[k][j] + (q * this.H[k + 1][j]);
                            if (notlast) {
                                p += (r * this.H[k + 2][j]);
                                this.H[k + 2][j] -= (p * z);
                            }
                            this.H[k][j] -= (p * x);
                            this.H[k + 1][j] -= (p * y);
                        }

                        // Column modification
                        for (int i = 0; i <= Math.min(n1, k + 3); i++) {
                            p = (x * this.H[i][k]) + (y * this.H[i][k + 1]);
                            if (notlast) {
                                p += (z * this.H[i][k + 2]);
                                this.H[i][k + 2] -= (p * r);
                            }
                            this.H[i][k] -= p;
                            this.H[i][k + 1] -= (p * q);
                        }

                        // Accumulate transformations
                        for (int i = low; i <= high; i++) {
                            p = (x * this.V[i][k]) + (y * this.V[i][k + 1]);
                            if (notlast) {
                                p += (z * this.V[i][k + 2]);
                                this.V[i][k + 2] -= (p * r);
                            }
                            this.V[i][k] -= p;
                            this.V[i][k + 1] -= (p * q);
                        }
                    } // (s != 0)
                } // k loop
            } // check convergence
        } // while (n >= low)

        // Backsubstitute to find vectors of upper triangular form
        if (norm == 0D) {
            return;
        }

        for (n1 = nn - 1; n1 >= 0; n1--) {
            p = this.d[n1];
            q = this.e[n1];

            // Real vector
            if (q == 0) {
                int l = n1;
                this.H[n1][n1] = 1.0;
                for (int i = n1 - 1; i >= 0; i--) {
                    var w = this.H[i][i] - p;
                    r = 0.0;
                    for (int j = l; j <= n1; j++) {
                        r += (this.H[i][j] * this.H[j][n1]);
                    }
                    if (this.e[i] < 0D) {
                        z = w;
                        s = r;
                    } else {
                        l = i;
                        if (this.e[i] == 0D) {
                            this.H[i][n1] = -r / ((w != 0D) ? w : (Maths.eps * norm));
                            // Solve real equations
                        } else {
                            var x = this.H[i][i + 1];
                            var y = this.H[i + 1][i];
                            q = ((this.d[i] - p) * (this.d[i] - p)) + (this.e[i] * this.e[i]);
                            final var t = ((x * s) - (z * r)) / q;
                            this.H[i][n1] = t;
                            this.H[i + 1][n1] = (Math.abs(x) > Math.abs(z)) ? ((-r - (w * t)) / x) : ((-s - (y * t)) / z);
                        }

                        // Overflow control
                        var t = Math.abs(this.H[i][n1]);
                        if ((Maths.eps * (t * t)) > 1) {
                            for (int j = i; j <= n1; j++) {
                                this.H[j][n1] /= t;
                            }
                        }
                    }
                }
            }
            // Complex vector
            else if (q < 0) {
                int l = n1 - 1;

                // Last vector component imaginary so matrix is triangular
                if (Math.abs(this.H[n1][n1 - 1]) > Math.abs(this.H[n1 - 1][n1])) {
                    this.H[n1 - 1][n1 - 1] = q / this.H[n1][n1 - 1];
                    this.H[n1 - 1][n1] = -(this.H[n1][n1] - p) / this.H[n1][n1 - 1];
                } else {
                    this.cdiv(0.0, -this.H[n1 - 1][n1], this.H[n1 - 1][n1 - 1] - p, q);
                    this.H[n1 - 1][n1 - 1] = this.cdivr;
                    this.H[n1 - 1][n1] = this.cdivi;
                }
                this.H[n1][n1 - 1] = 0D;
                this.H[n1][n1] = 1D;
                for (int i = n1 - 2; i >= 0; i--) {
                    final double[] rowHi = this.H[i], rowHi2 = this.H[i + 1];
                    double ra = 0D, sa = 0D;
                    for (int j = l; j <= n1; j++) {
                        ra += (rowHi[j] * this.H[j][n1 - 1]);
                        sa += (rowHi[j] * this.H[j][n1]);
                    }
                    var w = this.H[i][i] - p;

                    if (this.e[i] < 0.0) {
                        z = w;
                        r = ra;
                        s = sa;
                    } else {
                        l = i;
                        if (this.e[i] == 0) {
                            this.cdiv(-ra, -sa, w, q);
                            rowHi[n1 - 1] = this.cdivr;
                            rowHi[n1] = this.cdivi;
                        } else {
                            // Solve complex equations
                            var x = rowHi[i + 1];
                            var y = rowHi2[i];
                            var vr = (((this.d[i] - p) * (this.d[i] - p)) + (this.e[i] * this.e[i])) - (q * q);
                            var vi = (this.d[i] - p) * 2D * q;
                            if ((vr == 0D) && (vi == 0D)) {
                                vr = Maths.eps * norm * (Math.abs(w) + Math.abs(q) + Math.abs(x) + Math.abs(y) + Math.abs(z));
                            }
                            this.cdiv(((x * r) - (z * ra)) + (q * sa), (x * s) - (z * sa) - (q * ra), vr, vi);
                            rowHi[n1 - 1] = this.cdivr;
                            rowHi[n1] = this.cdivi;
                            if (Math.abs(x) > Math.abs(z) + Math.abs(q)) {
                                rowHi2[n1 - 1] = ((-ra - (w * rowHi[n1 - 1])) + (q * rowHi[n1])) / x;
                                rowHi2[n1] = (-sa - (w * rowHi[n1]) - (q * rowHi[n1 - 1])) / x;
                            } else {
                                this.cdiv(-r - (y * rowHi[n1 - 1]), -s - (y * rowHi[n1]), z, q);
                                rowHi2[n1 - 1] = this.cdivr;
                                rowHi2[n1] = this.cdivi;
                            }
                        }

                        // Overflow control
                        final var t = Math.max(Math.abs(rowHi[n1 - 1]), Math.abs(rowHi[n1]));
                        if ((Maths.eps * (t * t)) > 1) {
                            for (int j = i; j <= n1; j++) {
                                this.H[j][n1 - 1] /= t;
                                this.H[j][n1] /= t;
                            }
                        }
                    }
                }
            }
        }

        // Vectors of isolated roots
        for (int i = 0; i < nn; i++) {
            if ((i < low) || (i > high)) {
                for (int j = i; j < nn; j++) {
                    this.V[i][j] = this.H[i][j];
                }
            }
        }

        // Back transformation to get eigenvectors of original matrix
        for (int j = nn - 1; j >= low; j--) {
            for (int i = low; i <= high; i++) {
                var tmp = 0D;
                for (int k = low; k <= Math.min(j, high); k++) {
                    tmp += (this.V[i][k] * this.H[k][j]);
                }
                this.V[i][j] = tmp;
            }
        }
    }

    /**
     * Nonsymmetric reduction to Hessenberg form.
     * 
     * This is derived from the Algol procedures orthes and ortran, by Martin
     * and Wilkinson, Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the
     * corresponding Fortran subroutines in EISPACK.
     */
    private void orthes() {
        this.ort = new double[this.n];
        final int low = 0, high = this.n - 1;

        for (int m = low + 1; m <= (high - 1); m++) {
            // Scale column.
            var scale = 0D;
            for (int i = m; i <= high; i++) {
                scale += Math.abs(this.H[i][m - 1]);
            }
            if (scale != 0D) {
                // Compute Householder transformation.
                double h = 0D;
                for (int i = high; i >= m; i--) {
                    this.ort[i] = (this.H[i][m - 1] / scale);
                    h += (this.ort[i] * this.ort[i]);
                }
                var g = Math.sqrt(h);
                if (this.ort[m] > 0) {
                    g = -g;
                }
                
                h -= (this.ort[m] * g);
                this.ort[m] -= g;

                // Apply Householder similarity transformation
                // H = (I-u*u'/h)*H*(I-u*u')/h)
                for (int j = m; j < this.n; j++) {
                    double f = 0D;
                    for (int i = high; i >= m; i--) {
                        f += this.ort[i] * this.H[i][j];
                    }
                    f /= h;
                    for (int i = m; i <= high; i++) {
                        this.H[i][j] -= f * this.ort[i];
                    }
                }

                for (int i = 0; i <= high; i++) {
                    var f = 0D;
                    for (int j = high; j >= m; j--) {
                        f += this.ort[j] * this.H[i][j];
                    }
                    f /= h;
                    for (int j = m; j <= high; j++) {
                        this.H[i][j] -= f * this.ort[j];
                    }
                }
                this.ort[m] *= scale;
                this.H[m][m - 1] = scale * g;
            }
        }

        // Accumulate transformations (Algol's ortran).
        for (int i = 0; i < this.n; i++) {
            for (int j = 0; j < this.n; j++) {
                this.V[i][j] = (i == j ? 1.0 : 0.0);
            }
        }

        for (int m = high - 1; m > low; m--) {
            if (this.H[m][m - 1] != 0D) {
                for (int i = m + 1; i <= high; i++) {
                    this.ort[i] = this.H[i][m - 1];
                }
                for (int j = m; j <= high; j++) {
                    double g = 0D;
                    for (int i = m; i <= high; i++) {
                        g += this.ort[i] * this.V[i][j];
                    }
                    // Double division avoids possible underflow
                    g /= this.ort[m];
                    g /= this.H[m][m - 1];
                    for (int i = m; i <= high; i++) {
                        this.V[i][j] += g * this.ort[i];
                    }
                }
            }
        }
    }

    /**
     * Symmetric tridiagonal QL algorithm.
     * 
     * This is derived from the Algol procedures tql2, by Bowdler, Martin,
     * Reinsch, and Wilkinson, Handbook for Auto. Comp., Vol.ii-Linear Algebra,
     * and the corresponding Fortran subroutine in EISPACK.
     */
    private void tql2() {
        for (int i = 1; i < this.n; i++) {
            this.e[i - 1] = this.e[i];
        }
        this.e[this.n - 1] = 0D;

        double f = 0D, tst1 = 0D;
        for (int l = 0; l < this.n; l++) {
            // Find small subdiagonal element
            final double tst2 = Math.abs(this.d[l]) + Math.abs(this.e[l]);
            if (tst1 < tst2) {
                tst1 = tst2;
            }
            int m = l;
            while (m < this.n) {
                if (Math.abs(this.e[m]) <= (Maths.eps * tst1)) {
                    break;
                } else {
                    m++;
                }
            }

            // If m == l, d[l] is an eigenvalue,
            // otherwise, iterate.
            if (m > l) {
                do {
                    // Compute implicit shift
                    double g = this.d[l], p = (this.d[l + 1] - g) / (2D * this.e[l]),
                            r = Maths.hypot(p, 1D);
                    if (p < 0) {
                        r = -r;
                    }
                    final double dl1 = this.d[l + 1] = (this.e[l] * (p + r));
                    double h = g - (this.d[l] = (this.e[l] / (p + r)));
                    for (int i = l + 2; i < this.n; i++) {
                        this.d[i] -= h;
                    }
                    f += h;

                    // Implicit QL transformation.
                    p = this.d[m];
                    double c = 1D, c2 = c, c3 = c, s = 0D, s2 = 0D;
                    final double el1 = this.e[l + 1];
                    for (int i = m - 1; i >= l; i--) {
                        c3 = c2;
                        c2 = c;
                        s2 = s;
                        g = c * this.e[i];
                        h = c * p;
                        r = Maths.hypot(p, this.e[i]);
                        this.e[i + 1] = s * r;
                        s = this.e[i] / r;
                        c = p / r;
                        p = (c * this.d[i]) - (s * g);
                        this.d[i + 1] = h + (s * ((c * g) + (s * this.d[i])));

                        // Accumulate transformation.
                        for (final double[] rowV : this.V) {
                            h = rowV[i + 1];
                            rowV[i + 1] = (s * rowV[i]) + (c * h);
                            rowV[i] = (c * rowV[i]) - (s * h);
                        }
                    }
                    p = -(s * s2 * c3 * el1 * this.e[l]) / dl1;
                    this.e[l] = s * p;
                    this.d[l] = c * p;
                    // Check for convergence.
                } while (Math.abs(this.e[l]) > (Maths.eps * tst1));
            }
            this.d[l] += f;
            this.e[l] = 0D;
        }

        // Sort eigenvalues and corresponding vectors.
        for (int i = 0; i < (this.n - 1); i++) {
            double p = this.d[i];

            int k = i;
            for (int j = i + 1; j < this.n; j++) {
                if (this.d[j] < p) {
                    k = j;
                    p = this.d[j];
                }
            }

            if (k != i) {
                this.d[k] = this.d[i];
                this.d[i] = p;
                for (final double[] rowV : this.V) {
                    // swap rowV[i] and rowV[k]
                    p = rowV[i];
                    rowV[i] = rowV[k];
                    rowV[k] = p;
                }
            }
        }
    }

    /**
     * Symmetric Householder reduction to tridiagonal form.
     * 
     * This is derived from the Algol procedures tred2 by Bowdler, Martin,
     * Reinsch, and Wilkinson, Handbook for Auto. Comp., Vol.ii-Linear Algebra,
     * and the corresponding Fortran subroutine in EISPACK.
     * 
     */
    private void tred2() {
        for (int j = 0; j < this.n; j++) {
            this.d[j] = this.V[this.n - 1][j];
        }

        // Householder reduction to tridiagonal form.
        for (int i = this.n - 1; i > 0; i--) {
            // Scale to avoid under/overflow.
            double scale = 0.0, h = 0.0;
            for (int k = 0; k < i; k++) {
                scale += Math.abs(this.d[k]);
            }
            if (scale == 0.0) {
                this.e[i] = this.d[i - 1];
                for (int j = 0; j < i; j++) {
                    this.d[j] = this.V[i - 1][j];
                    this.V[i][j] = this.V[j][i] = 0D;
                }
            } else {
                // Generate Householder vector.
                for (int k = 0; k < i; k++) {
                    this.d[k] /= scale;
                    h += this.d[k] * this.d[k];
                }
                double f = this.d[i - 1], g = Math.sqrt(h);
                if (f > 0D) {
                    g = -g;
                }
                this.e[i] = scale * g;
                h -= (f * g);
                this.d[i - 1] = f - g;
                for (int j = 0; j < i; j++) {
                    this.e[j] = 0D;
                }

                // Apply similarity transformation to remaining columns.
                for (int j = 0; j < i; j++) {
                    f = this.d[j];
                    this.V[j][i] = f;
                    g = this.e[j] + (this.V[j][j] * f);
                    for (int k = j + 1; k <= (i - 1); k++) {
                        final double[] rowV = this.V[k];
                        g += rowV[j] * this.d[k];
                        this.e[k] += rowV[j] * f;
                    }
                    this.e[j] = g;
                }
                f = 0D;
                for (int j = 0; j < i; j++) {
                    f += (this.e[j] /= h) * this.d[j];
                }
                final double hh = f / (h + h);
                for (int j = 0; j < i; j++) {
                    this.e[j] -= hh * this.d[j];
                }
                for (int j = 0; j < i; j++) {
                    f = this.d[j];
                    g = this.e[j];
                    for (int k = j; k <= (i - 1); k++) {
                        this.V[k][j] -= ((f * this.e[k]) + (g * this.d[k]));
                    }
                    this.d[j] = this.V[i - 1][j];
                    this.V[i][j] = 0.0;
                }
            }
            this.d[i] = h;
        }

        // Accumulate transformations.
        for (int i = 0; i < (this.n - 1); i++) {
            this.V[this.n - 1][i] = this.V[i][i];
            this.V[i][i] = 1.0;
            final double h = this.d[i + 1];
            if (h != 0.0) {
                for (int k = 0; k <= i; k++) {
                    this.d[k] = this.V[k][i + 1] / h;
                }
                for (int j = 0; j <= i; j++) {
                    double g = 0.0;
                    for (int k = 0; k <= i; k++) {
                        final double[] rowV = this.V[k];
                        g += rowV[i + 1] * rowV[j];
                    }
                    for (int k = 0; k <= i; k++) {
                        this.V[k][j] -= g * this.d[k];
                    }
                }
            }
            for (int k = 0; k <= i; k++) {
                this.V[k][i + 1] = 0.0;
            }
        }
        for (int j = 0; j < this.n; j++) {
            this.d[j] = this.V[this.n - 1][j];
            this.V[this.n - 1][j] = 0.0;
        }
        this.V[this.n - 1][this.n - 1] = 1.0;
        this.e[0] = 0.0;
    }
}
