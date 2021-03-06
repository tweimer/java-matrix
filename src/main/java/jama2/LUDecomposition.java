package jama2;

import java.io.Serializable;
import java.util.Arrays;

/**
 * LU Decomposition.
 * <P>
 * For an m-by-n matrix A with m &gt;= n, the LU decomposition is an m-by-n unit
 * lower triangular matrix L, an n-by-n upper triangular matrix U, and a
 * permutation vector piv of length m so that A(piv,:) = L*U. If m &lt; n, then
 * L is m-by-m and U is m-by-n.
 * </P>
 * <P>
 * The LU decompostion with pivoting always exists, even if the matrix is
 * singular, so the constructor will never fail. The primary use of the LU
 * decomposition is in the solution of square systems of simultaneous linear
 * equations. This will fail if isNonsingular() returns false.
 * </P>
 * 
 * @author The MathWorks, Inc. and the National Institute of Standards and
 *         Technology.
 * @version 2.0
 * @see <a href="http://tweimer.github.io/java-matrix/">java-matrix</a>
 */
public class LUDecomposition implements ISolver, FunctionalMatrix, Serializable {
    /**
     * For the Serializeable interface
     */
    private static final long serialVersionUID = 1;

    /**
     * Array for internal storage of decomposition.
     * 
     * @serial internal array storage.
     */
    private final double[][] LU;

    /**
     * Row and column dimensions and pivot sign.
     * 
     * @serial column dimension.
     * @serial row dimension.
     * @serial pivot sign.
     */
    private final int m, n, pivsign;

    /**
     * Internal storage of pivot vector.
     * 
     * @serial pivot vector.
     */
    private final int piv[];


    /**
     * LU Decomposition, computed by Gaussian elimination.
     * <P>
     * This constructor computes L and U with the "daxpy"-based
     * elimination algorithm used in LINPACK and MATLAB. In
     * Java, we suspect the dot-product, Crout algorithm will be
     * faster. We have temporarily included this constructor
     * until timing experiments confirm this suspicion.
     * </P>
     * Structure to access L, U and piv.
     *
     * @param A
     *            Rectangular matrix
     * @param linpackflag
     *            If true, use Gaussian elimination,
     *            otherwise Crout algorithm.
     */
    LUDecomposition(final Matrix A, final boolean linpackflag) {
        // Initialize.
        LU = A.getArrayCopy();
        m = A.getRowDimension();
        n = A.getColumnDimension();
        piv = new int[this.m];
        
        for (var i = 1; i < m; i++) {
            piv[i] = i;
        }
        
        var pivsign = 1;
        if (linpackflag) {
            // Main loop.
            for (var k = 0; k < n; k++) {
                // Find pivot.
                var p = k;
                for (var i = k + 1; i < m; i++) {
                    if (Math.abs(LU[i][k]) > Math.abs(LU[p][k])) {
                        p = i;
                    }
                }
                
                // Exchange if necessary.
                if (p > k) {
                    for (var j = 0; j < n; j++) {
                        final var t = LU[p][j];
                        LU[p][j] = LU[k][j];
                        LU[k][j] = t;
                    }
                    
                    // Swap piv[p] and piv[k]
                    final var t = piv[p];
                    piv[p] = piv[k];
                    piv[k] = t;
                    
                    // alternate pivsign
                    pivsign = -pivsign;
                }
                
                // Compute multipliers and eliminate k-th column.
                if (LU[k][k] != 0D) {
                    for (var i = k + 1; i < m; i++) {
                        LU[i][k] /= LU[k][k];
                        for (var j = k + 1; j < n; j++) {
                            LU[i][j] -= LU[i][k] * LU[k][j];
                        }
                    }
                }
            }
        } else {
            // Outer loop.
            for (var j = 0; j < n; j++) {
                final var LUcolj = new double[this.m];

                // Make a copy of the j-th column to localize references.
                for (var i = 0; i < this.m; i++) {
                    LUcolj[i] = this.LU[i][j];
                }

                // Apply previous transformations.
                for (var i = 0; i < m; i++) {
                    final var LUrowi = LU[i];

                    // Most of the time is spent in the following dot product.
                    final var kmax = Math.min(i, j);
                    var s = 0D;
                    for (var k = 0; k < kmax; k++) {
                        s += LUrowi[k] * LUcolj[k];
                    }
                    LUrowi[j] = LUcolj[i] -= s;
                }

                // Find pivot and exchange if necessary.
                var p = j;
                for (var i = j + 1; i < m; i++) {
                    if (Math.abs(LUcolj[i]) > Math.abs(LUcolj[p])) {
                        p = i;
                    }
                }

                if (p > j) {
                    // Swap this.LU[p][k] and this.LU[j][k]
                    final var row_p = LU[p];
                    LU[p] = LU[j];
                    LU[j] = row_p;

                    // swap piv[p] and piv[j]
                    final var k = piv[p];
                    piv[p] = piv[j];
                    piv[j] = k;

                    // alternate pivsign
                    pivsign = -pivsign;
                }

                // Compute multipliers.
                if ((j < m) && (LU[j][j] != 0D)) {
                    for (var i = j + 1; i < m; i++) {
                        LU[i][j] /= LU[j][j];
                    }
                }
            }
        }
        this.pivsign = pivsign;
    }


    /**
     * LU Decomposition Structure to access L, U and piv.
     * 
     * <p>This is a package-private constructor.
     * Use {@link Matrix#lu()} to create a cholesky decomposition of a given matrix.</p>
     * 
     * @param A
     *            Rectangular matrix
     * @see Matrix#lu()
     */
    LUDecomposition(final Matrix A) {
        this(A, false);
    }

    /**
     * Determinant of a square matrix
     * 
     * @return det(A) for square matrices, NaN otherwise. Iff det equals 0,
     *         matrix is singular.
     */
    public double det() {
        if (m == n) {
            // go through the LU decomposition
            // For A = P*L*U, det(A)=det(P)*det(L)*det(U), where
            // det(P) is pivsign
            // det(L) is equal to 1
            // det(U) is the product of the diagonal entries
            // So lets start with pivsign and multiply with the diagonal entries.
            double det = pivsign;
            for (var j = 0; j < this.n; j++) {
                det *= this.LU[j][j];
            }
            return det;
        } else {
            return Double.NaN;
        }
    }

    /**
     * Return pivot permutation vector as a one-dimensional double array
     * 
     * @return (double) piv
     */
    public double[] getDoublePivot() {
        final var vals = new double[this.m];
        for (var i = 0; i < this.m; i++) {
            vals[i] = this.piv[i];
        }
        return vals;
    }

    /**
     * Return lower triangular factor
     * 
     * @return L lower triangular factor
     */
    public Matrix getL() {
        final var L = new double[this.m][this.n];
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < i; j++) {
                L[i][j] = this.LU[i][j];
            }

            L[i][i] = 1D;
        }
        return new Matrix(this.m, this.n, L);
    }

    /**
     * Return pivot permutation vector
     * 
     * @return piv pivot permutation vector
     */
    public int[] getPivot() {
        return Arrays.copyOf(this.piv, this.m);
    }

    /**
     * Return upper triangular factor
     * 
     * @return U upper triangular factor
     */
    public Matrix getU() {
        final var U = new double[this.n][this.n];
        for (var i = 0; i < this.n; i++) {
            for (var j = i; j < this.n; j++) {
                U[i][j] = this.LU[i][j];
            }
        }
        return new Matrix(this.n, U);
    }

    /**
     * Is the matrix nonsingular?
     * 
     * @return true if U, and hence A, is nonsingular.
     */
    public boolean isNonsingular() {
        for (var j = 0; j < this.n; j++) {
            if (this.LU[j][j] == 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * Solve A*X = B
     * 
     * @param B
     *            A Matrix with as many rows as A and any number of columns.
     * @return Returns null if matrix row dimensions don't agree or A is
     *         singular. Returns X so that <code>L*U*X = B(piv,:)</code> otherwise.
     * @throws NullPointerException Iff <code>B == null</code>
     */
    @Override
    public Matrix solve(final Matrix B) {
        if ((B.getRowDimension() == this.m) && this.isNonsingular()) {
            // Copy right hand side with pivoting
            final int nx = B.getColumnDimension();
            final var Xmat = B.getMatrix(this.piv, 0, nx - 1);
            final var X = Xmat.getArray();

            // Solve L*Y = B(piv,:)
            for (int k = 0; k < this.n; k++) {
                for (int i = k + 1; i < this.n; i++) {
                    for (int j = 0; j < nx; j++) {
                        X[i][j] -= X[k][j] * this.LU[i][k];
                    }
                }
            }
            // Solve U*X = Y;
            for (int k = this.n - 1; k >= 0; k--) {
                for (int j = 0; j < nx; j++) {
                    // LU[k][k] != 0
                    X[k][j] /= this.LU[k][k];
                }
                for (int i = 0; i < k; i++) {
                    for (int j = 0; j < nx; j++) {
                        X[i][j] -= X[k][j] * this.LU[i][k];
                    }
                }
            }
            return Xmat;
        } else {
            return null;
        }
    }


    /**
     * Returns the elements of the LUDecomposition
     * @param i row index
     * @param j column index
     * @return value
     */
	@Override
	public double get(final int i, final int j) {
		return this.LU[i][j];
	}
}
