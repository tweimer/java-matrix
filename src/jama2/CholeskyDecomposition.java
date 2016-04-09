package jama2;

import java.io.Serializable;

/**
 * Cholesky Decomposition.
 * <P>
 * For a symmetric, positive definite matrix A, the Cholesky decomposition is an
 * lower triangular matrix L so that A = L*L'.
 * <P>
 * If the matrix is not symmetric or positive definite, the constructor returns
 * a partial decomposition and sets an internal flag that may be queried by the
 * isSPD() method.
 * 
 * 
 * @author The MathWorks, Inc. and the National Institute of Standards and
 *         Technology.
 * @version 2.0
 * @see http://tweimer.github.io/jama/
 */
public class CholeskyDecomposition implements Serializable
{
    /**
     * Array for internal storage of decomposition.
     * 
     * @serial internal array storage.
     */
    private final double[][] L;

    /**
     * Row and column dimension (square matrix).
     * 
     * @serial matrix dimension.
     */
    private final int n;

    /**
     * Symmetric and positive definite flag.
     * 
     * @serial is symmetric and positive definite flag.
     */
    private boolean isspd;

    private static final long serialVersionUID = 1;

    //    /* ------------------------
    //       Temporary, experimental code.
    //     * ------------------------ */
    //
    //    /**
    //     * Right Triangular Cholesky Decomposition.
    //     * <P>
    //     * For a symmetric, positive definite matrix A, the Right Cholesky
    //     * decomposition is an upper triangular matrix R so that A = R'*R. This
    //     * constructor computes R with the Fortran inspired column oriented
    //     * algorithm used in LINPACK and MATLAB. In Java, we suspect a row oriented,
    //     * lower triangular decomposition is faster. We have temporarily included
    //     * this constructor here until timing experiments confirm this suspicion.
    //     */
    //
    //    /** Array for internal storage of right triangular decomposition. **/
    //    private transient double[][] R;
    //
    //    /**
    //     * Cholesky algorithm for symmetric and positive definite matrix.
    //     * 
    //     * @param Arg
    //     *            Square, symmetric matrix.
    //     * @param rightflag
    //     *            Actual value ignored.
    //     * @return Structure to access R and isspd flag.
    //     */
    //    public CholeskyDecomposition(final Matrix Arg, final int rightflag)
    //    {
    //        this.L = null;
    //        // Initialize.
    //        final double[][] A = Arg.getArray();
    //        this.n = Arg.getColumnDimension();
    //        this.R = new double[this.n][this.n];
    //        this.isspd = (Arg.getColumnDimension() == this.n);
    //        // Main loop.
    //        for (int j = 0; j < this.n; j++)
    //        {
    //            double d = 0.0;
    //            for (int k = 0; k < j; k++)
    //            {
    //                double s = A[k][j];
    //                for (int i = 0; i < k; i++)
    //                {
    //                    s -= this.R[i][k] * this.R[i][j];
    //                }
    //                s /= this.R[k][k];
    //                this.R[k][j] = s;
    //                d += s * s;
    //                this.isspd = this.isspd && (A[k][j] == A[j][k]);
    //            }
    //            d = A[j][j] - d;
    //            this.isspd = this.isspd && (d > 0.0);
    //            this.R[j][j] = Math.sqrt(Math.max(d, 0.0));
    //            for (int k = j + 1; k < this.n; k++)
    //            {
    //                this.R[k][j] = 0.0;
    //            }
    //        }
    //    }
    //
    //    /**
    //     * Return upper triangular factor.
    //     * 
    //     * @return R
    //     */
    //    public Matrix getR()
    //    {
    //        return new Matrix(this.R, this.n, this.n);
    //    }
    //
    //    /* ------------------------
    //       End of temporary code.
    //     * ------------------------ */

    /**
     * Cholesky algorithm for symmetric and positive definite matrix. Structure
     * to access L and isspd flag.
     * 
     * @param Arg
     *            Square, symmetric matrix.
     */
    public CholeskyDecomposition(final Matrix Arg)
    {
        // Initialize.
        final double[][] A = Arg.getArray();
        this.n = Arg.getRowDimension();
        this.L = new double[this.n][this.n];
        this.isspd = (Arg.getColumnDimension() == this.n);
        // Main loop.
        for (int j = 0; j < this.n; j++)
        {
            final double[] Lrowj = this.L[j];
            double d = 0D;
            for (int k = 0; k < j; k++)
            {
                final double[] Lrowk = this.L[k];
                double s = 0D;
                for (int i = 0; i < k; i++)
                {
                    s += Lrowk[i] * Lrowj[i];
                }
                Lrowj[k] = s = (A[j][k] - s) / this.L[k][k];
                d += (s * s);
                if (A[k][j] != A[j][k])
                {
                    this.isspd = false;
                }
            }
            d = A[j][j] - d;
            if (d <= 0D)
            {
                this.isspd = false;
                this.L[j][j] = 0D;
            }
            else
            {
                this.L[j][j] = Math.sqrt(d);
            }
            for (int k = j + 1; k < this.n; k++)
            {
                this.L[j][k] = 0D;
            }
        }
    }

    /**
     * Return triangular factor.
     * 
     * @return L
     */
    public Matrix getL()
    {
        return new Matrix(this.n, this.n, this.L);
    }

    /**
     * Is the matrix symmetric and positive definite?
     * 
     * @return true if A is symmetric and positive definite.
     */
    public boolean isSPD()
    {
        return this.isspd;
    }

    /**
     * Solve A*X = B
     * 
     * @param B
     *            A Matrix with as many rows as A and any number of columns.
     * @return X so that L*L'*X = B
     * @exception IllegalArgumentException
     *                Matrix row dimensions must agree.
     * @exception RuntimeException
     *                Matrix is not symmetric positive definite.
     */
    public Matrix solve(final Matrix B)
    {
        if (B.getRowDimension() != this.n)
        {
            throw new IllegalArgumentException("Matrix row dimensions must agree."); //$NON-NLS-1$
        }
        else if (!this.isspd)
        {
            throw new RuntimeException("Matrix is not symmetric positive definite."); //$NON-NLS-1$
        }

        // Copy right hand side.
        final Matrix M = new Matrix(B);
        final double[][] X = M.getArray();
        final int nx = B.getColumnDimension();

        // Solve L*Y = B;
        for (int k = 0; k < this.n; k++)
        {
            for (int j = 0; j < nx; j++)
            {
                for (int i = 0; i < k; i++)
                {
                    X[k][j] -= X[i][j] * this.L[k][i];
                }
                X[k][j] /= this.L[k][k];
            }
        }

        // Solve L'*X = Y;
        for (int k = this.n - 1; k >= 0; k--)
        {
            for (int j = 0; j < nx; j++)
            {
                for (int i = k + 1; i < this.n; i++)
                {
                    X[k][j] -= X[i][j] * this.L[i][k];
                }
                X[k][j] /= this.L[k][k];
            }
        }

        return M;
    }

}
