package jama2;

import jama2.util.Maths;

import java.io.Serializable;

/**
 * QR Decomposition.
 * <P>
 * For an m-by-n matrix A with m &gt;= n, the QR decomposition is an m-by-n
 * orthogonal matrix Q and an n-by-n upper triangular matrix R so that A = Q*R.
 * <P>
 * The QR decompostion always exists, even if the matrix does not have full
 * rank, so the constructor will never fail. The primary use of the QR
 * decomposition is in the least squares solution of nonsquare systems of
 * simultaneous linear equations. This will fail if isFullRank() returns false.
 * 
 * @author The MathWorks, Inc. and the National Institute of Standards and
 *         Technology.
 * @version 2.0
 * @see <a href="http://tweimer.github.io/java-matrix/">java-matrix</a>
 */
public class QRDecomposition implements FunctionalMatrix, Serializable {
    /**
     * For the Serializable interface.
     */
    private static final long serialVersionUID = 1;

    /**
     * Array for internal storage of decomposition.
     *
     * @serial internal array storage.
     */
    private final double[][] QR;

    /**
     * Row and column dimensions.
     *
     * @serial column dimension.
     * @serial row dimension.
     */
    private final int nRows, nColumns;

    /**
     * Array for internal storage of diagonal of R.
     *
     * @serial diagonal of R.
     */
    private final double[] Rdiag;

    /**
     * QR Decomposition, computed by Householder reflections. Structure to
     * access R and the Householder vectors and compute Q.
     *
     * <p>This is a package-private constructor.
     * Use {@link Matrix#qr()} to create a cholesky decomposition of a given matrix.</p>
     * 
     * @param A
     *            Rectangular matrix
     *@see Matrix#qr()
     */
    QRDecomposition(final SizedMatrix A) {
        // Initialize.
        QR = A.getArrayCopy();
        nRows = A.getRowDimension();
        nColumns = A.getColumnDimension();
        Rdiag = new double[nColumns];

        // Main loop.
        for (var k = 0; k < nColumns; k++) {
            // Compute 2-norm of k-th column without under/overflow.
            var nrm = 0.0;
            for (var i = k; i < nRows; i++) {
                nrm = Maths.hypot(nrm, this.QR[i][k]);
            }

            if (nrm != 0D) {
                // Form k-th Householder vector.
                if (QR[k][k] < 0.0) {
                    nrm = -nrm;
                }
                for (var i = k; i < nRows; i++) {
                    QR[i][k] /= nrm;
                }
                QR[k][k]++;

                // Apply transformation to remaining columns.
                for (var j = k + 1; j < nColumns; j++) {
                    var s = 0.0;
                    for (var i = k; i < nRows; i++) {
                        s += QR[i][k] * QR[i][j];
                    }
                    s /= -QR[k][k];
                    for (var i = k; i < nRows; i++) {
                        QR[i][j] += s * QR[i][k];
                    }
                }
            }
            Rdiag[k] = -nrm;
        }
    }

    /**
     * Return the Householder vectors
     *
     * @return m &times; n lower trapezoidal matrix
     * whose columns define the reflections
     */
    public Matrix getH() {
        final var H = new double[nRows][nColumns];
        for (var i = 0; i < nRows; i++) {
            for (var j = 0; j <= i; j++) {
                H[i][j] = QR[i][j];
            }
        }
        return new Matrix(nRows, nColumns, H);
    }

    /**
     * Generate and return the (economy-sized) orthogonal factor
     *
     * @return m &times; n (economy-sized) orthogonal factor Q
     */
    public Matrix getQ() {
        final var Q = new double[nRows][nColumns];
        for (var k = nColumns - 1; k >= 0; k--) {
            for (var i = 0; i < nRows; i++) {
                Q[i][k] = 0D;
            }
            Q[k][k] = 1D;
            for (var j = k; j < nColumns; j++) {
                if (QR[k][k] != 0) {
                    var s = 0D;
                    for (var i = k; i < nRows; i++) {
                        s += QR[i][k] * Q[i][j];
                    }
                    s /= -QR[k][k];
                    for (var i = k; i < nRows; i++) {
                        Q[i][j] += s * QR[i][k];
                    }
                }
            }
        }
        return new Matrix(nRows, nColumns, Q);
    }

    /**
     * Return the upper triangular factor
     *
     * @return n &times; n upper triangular factor R
     */
    public Matrix getR() {
        final var R = new double[nColumns][nColumns];
        for (var i = 0; i < nColumns; i++) {
            for (var j = i; j < nColumns; j++) {
                R[i][j] = ((i < j) ? QR[i][j] : Rdiag[i]);
            }
        }
        return new Matrix(nColumns, R);
    }

    /**
     * Is the matrix full rank?
     * 
     * @return true if R, and hence A, has full rank.
     */
    public boolean isFullRank() {
        for (final var r : this.Rdiag) {
            if (r == 0D) {
                return false;
            }
        }
        return true;
    }

    /**
     * Least squares solution of A*X = B.
     *
     * @param B
     *            A Matrix with as many rows as A and any number of columns.
     * @return Returns null if row dimensions don't agree or matrix is rang
     *         deficient. Returns X that minimizes the two norm of Q*R*X-B
     *         Matrix row otherwise.
     */
    public Matrix solve(final Matrix B) {
        if (B.getRowDimension() != this.nRows) {
            return null;
        } else if (!this.isFullRank()) {
            return null;
        }

        // Copy right hand side
        final var M = new Matrix(B);
        final var X = M.getArray();

        // Compute Y = transpose(Q)*B
        final int nx = B.getColumnDimension();
        for (var k = 0; k < this.nColumns; k++) {
            for (var j = 0; j < nx; j++) {
                var s = 0D;
                for (var i = k; i < this.nRows; i++) {
                    s += this.QR[i][k] * X[i][j];
                }
                s /= -this.QR[k][k];
                for (var i = k; i < this.nRows; i++) {
                    X[i][j] += s * this.QR[i][k];
                }
            }
        }
        
        // Solve R*X = Y;
        for (var k = this.nColumns - 1; k >= 0; k--) {
            for (var j = 0; j < nx; j++) {
                X[k][j] /= this.Rdiag[k];
            }
            for (var i = 0; i < k; i++) {
                for (var j = 0; j < nx; j++) {
                    X[i][j] -= X[k][j] * this.QR[i][k];
                }
            }
        }
        return M.getMatrix(0, this.nColumns - 1, 0, nx - 1);
    }

    /**
     * Returns the elements of the QRDecomposition
     * @param i row index
     * @param j column index
     * @return value
     */
	@Override
	public double get(final int i, final int j) {
		return this.QR[i][j] ;
	}
}
