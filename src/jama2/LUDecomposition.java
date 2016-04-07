package jama2;

import java.io.Serializable;

/**
 * LU Decomposition.
 * <P>
 * For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n unit
 * lower triangular matrix L, an n-by-n upper triangular matrix U, and a
 * permutation vector piv of length m so that A(piv,:) = L*U. If m < n, then L
 * is m-by-m and U is m-by-n.
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
 * @version 8. Februar 2016
 * @see http://math.nist.gov/javanumerics/jama/
 */
public class LUDecomposition implements Serializable
{

    /**
     * Array for internal storage of decomposition.
     * 
     * @serial internal array storage.
     */
    private final double[][] LU;

    /**
     * Row and column dimensions, and pivot sign.
     * 
     * @serial column dimension.
     * @serial row dimension.
     * @serial pivot sign.
     */
    private final int m, n;

    /**
     * 
     */
    private int pivsign;

    /**
     * Internal storage of pivot vector.
     * 
     * @serial pivot vector.
     */
    private final int piv[];

    /**
     * 
     */
    private static final long serialVersionUID = 1;

    /* ------------------------
       Temporary, experimental code.
       ------------------------ */

    //    /**
    //     * LU Decomposition, computed by Gaussian elimination.
    //     * <P>
    //     * This constructor computes L and U with the "daxpy"-based elimination
    //     * algorithm used in LINPACK and MATLAB. In Java, we suspect the
    //     * dot-product, Crout algorithm will be faster. We have temporarily included
    //     * this constructor until timing experiments confirm this suspicion.
    //     * </P>
    //     * Structure to access L, U and piv.
    //     * 
    //     * @param A
    //     *            Rectangular matrix
    //     * @param linpackflag
    //     *            Use Gaussian elimination. Actual value ignored.
    //     */
    //    public LUDecomposition(final Matrix A, final int linpackflag)
    //    {
    //        // Initialize.
    //        this.LU = A.getArrayCopy();
    //        this.m = A.getRowDimension();
    //        this.n = A.getColumnDimension();
    //        this.piv = new int[this.m];
    //        for (int i = 0; i < this.m; i++)
    //        {
    //            this.piv[i] = i;
    //        }
    //        this.pivsign = 1;
    //        // Main loop.
    //        for (int k = 0; k < this.n; k++)
    //        {
    //            // Find pivot.
    //            int p = k;
    //            for (int i = k + 1; i < this.m; i++)
    //            {
    //                if (Math.abs(this.LU[i][k]) > Math.abs(this.LU[p][k]))
    //                {
    //                    p = i;
    //                }
    //            }
    //            // Exchange if necessary.
    //            if (p != k)
    //            {
    //                for (int j = 0; j < this.n; j++)
    //                {
    //                    final double t = this.LU[p][j];
    //                    this.LU[p][j] = this.LU[k][j];
    //                    this.LU[k][j] = t;
    //                }
    //                final int t = this.piv[p];
    //                this.piv[p] = this.piv[k];
    //                this.piv[k] = t;
    //                this.pivsign = -this.pivsign;
    //            }
    //            // Compute multipliers and eliminate k-th column.
    //            if (this.LU[k][k] != 0.0)
    //            {
    //                for (int i = k + 1; i < this.m; i++)
    //                {
    //                    this.LU[i][k] /= this.LU[k][k];
    //                    for (int j = k + 1; j < this.n; j++)
    //                    {
    //                        this.LU[i][j] -= this.LU[i][k] * this.LU[k][j];
    //                    }
    //                }
    //            }
    //        }
    //    }

    /* ------------------------
       End of temporary code.
     * ------------------------ */

    /* ------------------------
       Public Methods
     * ------------------------ */

    /**
     * LU Decomposition Structure to access L, U and piv.
     * 
     * @param A
     *            Rectangular matrix
     */
    public LUDecomposition(final Matrix A)
    {
        // Use a "left-looking", dot-product, Crout/Doolittle algorithm.
        this.LU = A.getArrayCopy();
        this.m = A.getRowDimension();
        this.n = A.getColumnDimension();
        this.piv = new int[this.m];
        for (int i = 0; i < this.m; i++)
        {
            this.piv[i] = i;
        }
        this.pivsign = 1;

        // Outer loop.
        for (int j = 0; j < this.n; j++)
        {
            final double[] LUcolj = new double[this.m];

            // Make a copy of the j-th column to localize references.
            for (int i = 0; i < this.m; i++)
            {
                LUcolj[i] = this.LU[i][j];
            }

            // Apply previous transformations.
            for (int i = 0; i < this.m; i++)
            {
                final double[] LUrowi = this.LU[i];

                // Most of the time is spent in the following dot product.
                final int kmax = Math.min(i, j);
                double s = 0D;
                for (int k = 0; k < kmax; k++)
                {
                    s += LUrowi[k] * LUcolj[k];
                }
                LUrowi[j] = LUcolj[i] -= s;
            }

            // Find pivot and exchange if necessary.
            int p = j;
            for (int i = j + 1; i < this.m; i++)
            {
                if (Math.abs(LUcolj[i]) > Math.abs(LUcolj[p]))
                {
                    p = i;
                }
            }

            if (p != j)
            {
                for (int k = 0; k < this.n; k++)
                {
                    // Swap this.LU[p][k] and this.LU[j][k]
                    final double t = this.LU[p][k];
                    this.LU[p][k] = this.LU[j][k];
                    this.LU[j][k] = t;
                }

                // swap piv[p] and piv[j]
                final int k = this.piv[p];
                this.piv[p] = this.piv[j];
                this.piv[j] = k;

                // alternate pivsign
                this.pivsign = -this.pivsign;
            }

            // Compute multipliers.
            if ((j < this.m) && (this.LU[j][j] != 0D))
            {
                for (int i = j + 1; i < this.m; i++)
                {
                    this.LU[i][j] /= this.LU[j][j];
                }
            }
        }
    }

    /**
     * Determinant
     * 
     * @return det(A)
     * @exception IllegalArgumentException
     *                Matrix must be square
     */
    public double det()
    {
        if (this.m == this.n)
        {
            // go through the LU decomposition
            double d = this.pivsign;
            for (int j = 0; j < this.n; j++)
            {
                d *= this.LU[j][j];
            }
            return d;
        }
        else
        {
            throw new IllegalArgumentException("Matrix must be square."); //$NON-NLS-1$
        }
    }

    /**
     * Return pivot permutation vector as a one-dimensional double array
     * 
     * @return (double) piv
     */
    public double[] getDoublePivot()
    {
        final double vals[] = new double[this.m];
        for (int i = 0; i < this.m; i++)
        {
            vals[i] = this.piv[i];
        }
        return vals;
    }

    /**
     * Return lower triangular factor
     * 
     * @return L
     */
    public Matrix getL()
    {
        final Matrix X = new Matrix(this.m, this.n);
        final double L[][] = X.getArray();
        for (int i = 0; i < this.m; i++)
        {
            for (int j = 0; j < this.n; j++)
            {
                L[i][j] = ((i > j) ? this.LU[i][j] : ((i == j) ? 1D : 0D));
            }
        }
        return X;
    }

    /**
     * Return pivot permutation vector
     * 
     * @return piv
     */
    public int[] getPivot()
    {
        final int p[] = new int[this.m];
        for (int i = 0; i < this.m; i++)
        {
            p[i] = this.piv[i];
        }
        return p;
    }

    /**
     * Return upper triangular factor
     * 
     * @return U
     */
    public Matrix getU()
    {
        final Matrix X = new Matrix(this.n, this.n);
        final double[][] U = X.getArray();
        for (int i = 0; i < this.n; i++)
        {
            for (int j = i; j < this.n; j++)
            {
                U[i][j] = this.LU[i][j];
            }
        }
        return X;
    }

    /**
     * Is the matrix nonsingular?
     * 
     * @return true if U, and hence A, is nonsingular.
     */
    public boolean isNonsingular()
    {
        for (int j = 0; j < this.n; j++)
        {
            if (this.LU[j][j] == 0)
            {
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
     * @return Returns null if matrix row dimensions don't agree or matrix is
     *         singular. Returns X so that L*U*X = B(piv,:) otherwise.
     */
    public Matrix solve(final Matrix B)
    {
        if (B.getRowDimension() != this.m)
        {
            return null;
            //throw new IllegalArgumentException("Matrix row dimensions must agree."); //$NON-NLS-1$
        }
        else if (!this.isNonsingular())
        {
            return null;
            //throw new RuntimeException("Matrix is singular."); //$NON-NLS-1$
        }
        else
        {
            // Copy right hand side with pivoting
            final int nx = B.getColumnDimension();
            final Matrix Xmat = B.getMatrix(this.piv, 0, nx - 1);
            final double X[][] = Xmat.getArray();

            // Solve L*Y = B(piv,:)
            for (int k = 0; k < this.n; k++)
            {
                for (int i = k + 1; i < this.n; i++)
                {
                    for (int j = 0; j < nx; j++)
                    {
                        X[i][j] -= X[k][j] * this.LU[i][k];
                    }
                }
            }
            // Solve U*X = Y;
            for (int k = this.n - 1; k >= 0; k--)
            {
                for (int j = 0; j < nx; j++)
                {
                    X[k][j] /= this.LU[k][k];
                }
                for (int i = 0; i < k; i++)
                {
                    for (int j = 0; j < nx; j++)
                    {
                        X[i][j] -= X[k][j] * this.LU[i][k];
                    }
                }
            }
            return Xmat;
        }
    }
}
