package jama2;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.io.StreamTokenizer;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;
import java.util.Vector;

import static jama2.util.Maths.hypot;
import static java.lang.Math.abs;

/**
 * Jama = Java Matrix class.
 * <P>
 * The Java Matrix Class provides the fundamental operations of numerical linear
 * algebra. Various constructors create Matrices from two dimensional arrays of
 * double precision floating point numbers. Various "gets" and "sets" provide
 * access to submatrices and matrix elements. Several methods implement basic
 * matrix arithmetic, including matrix addition and multiplication, matrix
 * norms, and element-by-element array operations. Methods for reading and
 * printing matrices are also included. All the operations in this version of
 * the Matrix Class involve real matrices. Complex matrices may be handled in a
 * future version.
 * </P>
 * <P>
 * Five fundamental matrix decompositions, which consist of pairs or triples of
 * matrices, permutation vectors, and the like, produce results in five
 * decomposition classes. These decompositions are accessed by the Matrix class
 * to compute solutions of simultaneous linear equations, determinants, inverses
 * and other matrix functions. The five decompositions are:
 * </P>
 * <UL>
 * <LI>Cholesky Decomposition of symmetric, positive definite matrices.</LI>
 * <LI>LU Decomposition of rectangular matrices.</LI>
 * <LI>QR Decomposition of rectangular matrices.</LI>
 * <LI>Singular Value Decomposition of rectangular matrices.</LI>
 * <LI>Eigenvalue Decomposition of both symmetric and nonsymmetric square
 * matrices.</LI>
 * </UL>
 * <DL>
 * <DT><B>Example of use:</B></DT>
 * <DD>Solve a linear system A x = b and compute the residual norm, ||b - A x||.
 * <P>
 * 
 * <PRE>
 * double[][] vals = { { 1., 2., 3 }, { 4., 5., 6. }, { 7., 8., 10. } };
 * Matrix A = new Matrix(vals);
 * Matrix b = Matrix.random(3, 1);
 * Matrix x = A.solve(b);
 * Matrix r = A.times(x).minus(b);
 * double rnorm = r.normInf();
 * </PRE>
 * 
 * </DD>
 * </DL>
 * 
 * @author The MathWorks, Inc. and the National Institute of Standards and
 *         Technology.
 * @version 2.0
 * @see <a href="http://tweimer.github.io/java-matrix/">java-matrix</a>
 */
public class Matrix implements Serializable
{
    /**
     * For the Serializeable interface
     */
    private static final long serialVersionUID = 1;

    /**
     * Construct a matrix from a copy of a 2-D array.
     * 
     * @param A
     *            Two-dimensional array of doubles.
     * @return Matrix with copied array
     * @exception IllegalArgumentException
     *                All rows must have the same length
     */
    public static Matrix constructWithCopy(final double A[][])
    {
        final int m = A.length;
        final int n = A[0].length;
        final double C[][] = new double[m][];
        for (int i = 0; i < C.length; i++)
        {
			if (A[i].length != n)
			{
				throw new IllegalArgumentException("All rows must have the same length."); //$NON-NLS-1$
			}
			else
			{
				C[i] = Arrays.copyOf(A[i], n);
			}
        }
        return new Matrix(m, n, C);
    }

    /**
     * Creates a diagonal Matrix
     * 
     * @param d
     *            diagonal elements
     * @return diagonal Matrix
     */
    public static Matrix diag(final double... d)
    {
        final Matrix D = new Matrix(d.length);
        for (int i = 0; i < d.length; i++)
        {
            D.set(i, i, d[i]);
        }
        return D;
    }

    /**
     * Generate identity matrix
     * 
     * @param m
     *            Number of rows and colums.
     * @return An m-by-n matrix with ones on the diagonal and zeros elsewhere.
     */
    public static Matrix identity(final int m)
    {
        return Matrix.identity(m, m);
    }

    /**
     * Generate identity matrix
     * 
     * @param m
     *            Number of rows.
     * @param n
     *            Number of colums.
     * @return An m-by-n matrix with ones on the diagonal and zeros elsewhere.
     */
    public static Matrix identity(final int m, final int n)
    {
        final Matrix I = new Matrix(m, n);
        I.identity();
        return I;
    }

    /**
     * Generate matrix with random elements
     * 
     * @param n
     *            Number of rows and colums.
     * @return An n-by-n matrix with uniformly distributed random elements.
     */
    public static Matrix random(final int n)
    {
        return Matrix.random(n, n);
    }

    /**
     * Generate matrix with random elements
     * 
     * @param m
     *            Number of rows.
     * @param n
     *            Number of colums.
     * @return An m-by-n matrix with uniformly distributed random elements.
     */
    public static Matrix random(final int m, final int n)
    {
        final Matrix R = new Matrix(m, n);
        R.random();
        return R;
    }

    /**
     * Generate matrix with random elements
     * 
     * @param n
     *            Number of rows and colums.
     * @return An n-by-n matrix with uniformly distributed random elements.
     */
    public static Matrix randomInt(final int n)
    {
        return Matrix.randomInt(n, n);
    }

    /**
     * Generate matrix with random elements
     * 
     * @param m
     *            Number of rows.
     * @param n
     *            Number of colums.
     * @return An m-by-n matrix with uniformly distributed random elements.
     */
    public static Matrix randomInt(final int m, final int n)
    {
        final Matrix A = new Matrix(m, n);
        A.randomInt();
        return A;
    }

    /**
     * Read a matrix from a stream. The format is the same the print method, so
     * printed matrices can be read back in (provided they were printed using US
     * Locale). Elements are separated by whitespace, all the elements for each
     * row appear on a single line, the last row is followed by a blank line.
     * 
     * @param input
     *            the input stream.
     * @return Matrix
     * @throws IOException
     */
    public static Matrix read(final BufferedReader input) throws IOException
    {
        final StreamTokenizer tokenizer = new StreamTokenizer(input);

        // Although StreamTokenizer will parse numbers, it doesn't recognize
        // scientific notation (E or D); however, Double.valueOf does.
        // The strategy here is to disable StreamTokenizer's number parsing.
        // We'll only get whitespace delimited words, EOL's and EOF's.
        // These words should all be numbers, for Double.valueOf to parse.

        tokenizer.resetSyntax();
        tokenizer.wordChars(0, 255);
        tokenizer.whitespaceChars(0, ' ');
        tokenizer.eolIsSignificant(true);

        while (tokenizer.nextToken() == StreamTokenizer.TT_EOL)
        {
            // Ignore initial empty lines
        }

        if (tokenizer.ttype == StreamTokenizer.TT_EOF)
        {
            throw new IOException("Unexpected EOF on matrix read."); //$NON-NLS-1$
        }

        final Vector<Double> vD = new Vector<>();
        do
        {
            // Read & store 1st row.
            vD.addElement(new Double(tokenizer.sval));
        }
        while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);

        // Now we've got the number of columns!
        final int n = vD.size();
        double row[] = new double[n];
        int j = 0;
        for (final Double d : vD)
        {
            // extract the elements of the 1st row.
            row[j++] = d.doubleValue();
        }

        // Start storing rows instead of columns.
        final Vector<double[]> v = new Vector<>();
        v.addElement(row);
        while (tokenizer.nextToken() == StreamTokenizer.TT_WORD)
        {
            // While non-empty lines
            v.addElement(row = new double[n]);
            j = 0;
            do
            {
                if (j >= n)
                {
                    throw new IOException("Row " + v.size() + " is too long."); //$NON-NLS-1$ //$NON-NLS-2$
                }
                row[j++] = Double.valueOf(tokenizer.sval).doubleValue();
            }
            while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
            if (j < n)
            {
                throw new IOException("Row " + v.size() + " is too short."); //$NON-NLS-1$ //$NON-NLS-2$
            }
        }

        // Now we've got the number of rows.
        final int m = v.size();
        final double[][] A = new double[m][];

        // copy the rows out of the vector
        v.copyInto(A);
        return new Matrix(A);
    }

    /**
     * Array for internal storage of elements.
     * 
     * @serial internal array storage.
     */
    private final double A[][];

    /**
     * Row and column dimensions.
     * 
     * @serial row dimension.
     * @serial column dimension.
     */
    private final int m, n;

    /**
     * Construct a matrix from a one-dimensional packed array
     * 
     * @param vals
     *            One-dimensional array of doubles, packed by columns (ala
     *            Fortran).
     * @param m
     *            Number of rows.
     * @throws IllegalArgumentException
     *             Array length must be a multiple of m.
     */
    public Matrix(final double vals[], final int m)
    {
        this.n = (m != 0 ? vals.length / m : 0);
        if ((m * this.n) != vals.length)
        {
            throw new IllegalArgumentException("Array length must be a multiple of m."); //$NON-NLS-1$
        }
        else
        {
            this.A = new double[this.m = m][this.n];
            for (int i = 0; i < this.A.length; i++)
            {
                for (int j = 0; j < this.A[i].length; j++)
                {
                    this.A[i][j] = vals[i + (j * m)];
                }
            }
        }
    }

    /**
     * Construct a matrix from a 2-D array. This does <b>not</b> copy the given
     * array, it just stores it's reference.
     * 
     * @param A
     *            Two-dimensional array of doubles.
     * @exception IllegalArgumentException
     *                All rows must have the same length
     * @see #constructWithCopy
     */
    public Matrix(final double A[][])
    {
        // TODO: how to handle null here?
        this(A.length, A[0].length, A);

        // check if each row has the same length
        for (final double r[] : this.A)
        {
            if (r.length != this.n)
            {
                throw new IllegalArgumentException("All rows must have the same length."); //$NON-NLS-1$
            }
        }
    }

    /**
     * Construct a matrix quickly without checking arguments.
     * 
     * @param m
     *            Number of rows.
     * @param n
     *            Number of colums.
     * @param A
     *            Two-dimensional array of doubles.
     */
    public Matrix(final int m, final int n, final double A[][])
    {
        this.A = A;
        this.m = m;
        this.n = n;
    }

    /**
     * Construct a matrix quickly without checking arguments.
     * 
     * @param n
     *            Number of rows and colums.
     * @param A
     *            Two-dimensional array of doubles.
     */
    public Matrix(final int n, final double A[][])
    {
        this(n, n, A);
    }

    /**
     * Construct an n-by-n matrix of zeros.
     * 
     * @param n
     *            Number of rows and colums.
     */
    public Matrix(final int n)
    {
        this(n, n);
    }

    /**
     * Construct an m-by-n matrix of zeros.
     * 
     * @param m
     *            Number of rows.
     * @param n
     *            Number of colums
     */
    public Matrix(final int m, final int n)
    {
        this.A = new double[this.m = m][this.n = n];
    }

    /**
     * Construct an m-by-n constant matrix.
     * 
     * @param m
     *            Number of rows.
     * @param n
     *            Number of colums.
     * @param s
     *            Fill the matrix with this scalar value.
     */
    public Matrix(final int m, final int n, final double s)
    {
        this(m, n);
        this.set(s);
    }

    /**
     * Copies Matrix X
     * 
     * @param X
     *            Matrix to be copied
     */
    public Matrix(final Matrix X)
    {
        this(X.m, X.n, X.getArrayCopy());
    }

    /**
     * Element-by-element left division, C = A.\B
     * 
     * @param B
     *            another matrix
     * @return A.\B
     */
    public Matrix arrayLeftDivide(final Matrix B)
    {
        this.checkMatrixDimensions(B);
        final Matrix M = new Matrix(this.m, this.n);
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                M.A[i][j] = B.A[i][j] / this.A[i][j];
            }
        }
        return M;
    }

    /**
     * Element-by-element left division in place, A = A.\B
     * 
     * @param B
     *            another matrix
     */
    public void arrayLeftDivideEquals(final Matrix B)
    {
        this.checkMatrixDimensions(B);
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                this.A[i][j] = B.A[i][j] / this.A[i][j];
            }
        }
    }

    /**
     * Element-by-element right division, C = A./B
     * 
     * @param B
     *            another matrix
     * @return A./B
     */
    public Matrix arrayRightDivide(final Matrix B)
    {
        this.checkMatrixDimensions(B);
        final Matrix M = new Matrix(this.m, this.n);
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                M.A[i][j] = this.A[i][j] / B.A[i][j];
            }
        }
        return M;
    }

    /**
     * Element-by-element right division in place, A = A./B
     * 
     * @param B
     *            another matrix
     */
    public void arrayRightDivideEquals(final Matrix B)
    {
        this.checkMatrixDimensions(B);
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                this.A[i][j] /= B.A[i][j];
            }
        }
    }

    /**
     * Element-by-element multiplication, C = A.*B
     * 
     * @param B
     *            another matrix
     * @return A.*B
     * @see #arrayTimesEquals
     */
    public Matrix arrayTimes(final Matrix B)
    {
        this.checkMatrixDimensions(B);
        final Matrix M = new Matrix(this.m, this.n);
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                M.A[i][j] = this.A[i][j] * B.A[i][j];
            }
        }
        return M;
    }

    /**
     * Element-by-element multiplication in place, A = A.*B
     * 
     * @param B
     *            another matrix
     * @see #arrayTimes
     */
    public void arrayTimesEquals(final Matrix B)
    {
        this.checkMatrixDimensions(B);
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                this.A[i][j] *= B.A[i][j];
            }
        }
    }

    /**
     * Check if size(A) == size(B)
     * 
     * @param B
     * @throws IllegalArgumentException
     *             if they don't agree.
     */
    private void checkMatrixDimensions(final Matrix B)
    {
        if (!this.equalDimensions(B))
        {
            throw new IllegalArgumentException("Matrix dimensions must agree."); //$NON-NLS-1$
        }
    }

    /**
     * Compares dimension of both matrices
     * 
     * @param B
     *            another Matrix
     * @return true if both have the same dimension, false otherwise, returns
     *         false if B==null.
     */
    public boolean equalDimensions(final Matrix B)
    {
        return (B != null) && (B.m == this.m) && (B.n == this.n);
    }

    /**
     * Cholesky Decomposition
     * 
     * @return CholeskyDecomposition
     * @see CholeskyDecomposition
     */
    public CholeskyDecomposition chol()
    {
        return new CholeskyDecomposition(this);
    }

    /**
     * Matrix condition (2 norm)
     * 
     * @return ratio of largest to smallest singular value.
     * @see SingularValueDecomposition#cond()
     */
    public double cond()
    {
        return this.svd().cond();
    }

    /**
     * Returns a copy of the matrix
     * 
     * @return copy of this Matrix
     * @deprecated use {@link #Matrix(Matrix)}
     */
    @Deprecated
    public Matrix copy()
    {
        return new Matrix(this);
    }

    /**
     * Matrix determinant
     * 
     * @return determinant
     */
    public double det()
    {
        return this.lu().det();
    }

    /**
     * Eigenvalue Decomposition
     * 
     * @return EigenvalueDecomposition
     * @see EigenvalueDecomposition
     */
    public EigenvalueDecomposition eig()
    {
        return new EigenvalueDecomposition(this);
    }

    /**
     * Compares a Matrix to another Matrix
     * 
     * @param other
     *            another Matrix
     * @return true if other equals A
     */
    public boolean equals(final Matrix other)
    {
        return (other == this) || ((other != null) && ((this.m == other.m) && (this.n == other.n) && Arrays.deepEquals(this.A, other.A)));
    }

    /**
     * Overloads
     * 
     * @param obj
     *            another object
     * @return true if other equals A
     */
    @Override
    public boolean equals(final Object obj)
    {
        return (obj instanceof Matrix) && this.equals((Matrix) obj);
    }

    /**
     * Returns index at index i, row-by-row (MATLAB like)
     * 
     * @param i
     *            index
     * @return entry
     */
    public double get(final int i)
    {
        return this.A[i / this.m][i % this.n];
    }

    /**
     * Get a single element.
     * 
     * @param i
     *            Row index.
     * @param j
     *            Column index.
     * @return A(i,j)
     * @throws ArrayIndexOutOfBoundsException
     *             if indices are invalid
     */
    public double get(final int i, final int j)
    {
        return this.A[i][j];
    }

    /**
     * Access the internal two-dimensional array.
     * 
     * @return Pointer to the two-dimensional array of matrix elements.
     */
    public double[][] getArray()
    {
        return this.A;
    }

    /**
     * Copy the internal two-dimensional array.
     * 
     * @return Two-dimensional array copy of matrix elements.
     */
    public double[][] getArrayCopy()
    {
        final double C[][] = new double[this.m][];
        for (int i = 0; i < C.length; i++)
        {
            C[i] = Arrays.copyOf(this.A[i], this.n);
        }
        return C;
    }

    /**
     * Get column dimension.
     * 
     * @return n, the number of columns.
     */
    public int getColumnDimension()
    {
        return this.n;
    }

    /**
     * Make a one-dimensional column packed copy of the internal array.
     * 
     * @return Matrix elements packed in a one-dimensional array by columns.
     */
    public double[] getColumnPackedCopy()
    {
        final double vals[] = new double[this.m * this.n];
        for (int i = 0; i < this.m; i++)
        {
            for (int j = 0; j < this.n; j++)
            {
                vals[i + (j * this.m)] = this.A[i][j];
            }
        }
        return vals;
    }

    /**
     * Make a one-dimensional row packed copy of the internal array.
     * 
     * @return Matrix elements packed in a one-dimensional array by rows.
     */
    public double[] getRowPackedCopy()
    {
        final double vals[] = new double[this.m * this.n];
        int i = 0;
        for (final double row[] : this.A)
        {
            for (final double d : row)
            {
                vals[i++] = d;
            }
        }
        return vals;
    }

    /**
     * Get a submatrix.
     * 
     * @param i0
     *            Initial row index (inclusive)
     * @param i1
     *            Final row index (inclusive)
     * @param c
     *            Array of column indices.
     * @return A(i0:i1,c(:))
     * @throws ArrayIndexOutOfBoundsException
     *             Submatrix indices
     */
    public Matrix getMatrix(final int i0, final int i1, final int c[])
    {
        final int m1 = (i1 - i0) + 1, n1 = c.length;
        final Matrix M = new Matrix(m1, n1);
        for (int i = i0; i <= i1; i++)
        {
            for (int j = 0; j < n1; j++)
            {
                M.A[i - i0][j] = this.A[i][c[j]];
            }
        }
        return M;
    }

    /**
     * Get a submatrix.
     * 
     * @param i0
     *            Initial row index (inclusive)
     * @param i1
     *            Final row index (inclusive)
     * @param j0
     *            Initial column index (inclusive)
     * @param j1
     *            Final column index (inclusive)
     * @return A(i0:i1,j0:j1)
     * @exception ArrayIndexOutOfBoundsException
     *                Submatrix indices
     */
    public Matrix getMatrix(final int i0, final int i1, final int j0, final int j1)
    {
        final int m1 = (i1 - i0) + 1, n1 = (j1 - j0) + 1;
        final Matrix M = new Matrix(m1, n1);
        for (int i = i0; i <= i1; i++)
        {
            for (int j = j0; j <= j1; j++)
            {
                M.A[i - i0][j - j0] = this.A[i][j];
            }
        }
        return M;
    }

    /**
     * Get a submatrix.
     * 
     * @param r
     *            Array of row indices
     * @param j0
     *            Initial column index (inclusive)
     * @param j1
     *            Final column index (inclusive)
     * @return A(r(:),j0:j1)
     * @exception ArrayIndexOutOfBoundsException
     *                Submatrix indices
     */
    public Matrix getMatrix(final int r[], final int j0, final int j1)
    {
        final int m1 = r.length, n1 = (j1 - j0) + 1;
        final Matrix M = new Matrix(m1, n1);
        for (int i = 0; i < m1; i++)
        {
            for (int j = j0; j <= j1; j++)
            {
                M.A[i][j - j0] = this.A[r[i]][j];
            }
        }
        return M;
    }

    /**
     * Get a submatrix.
     * 
     * @param r
     *            Array of row indices.
     * @param c
     *            Array of column indices.
     * @return A(r(:),c(:))
     * @exception ArrayIndexOutOfBoundsException
     *                Submatrix indices
     */
    public Matrix getMatrix(final int r[], final int c[])
    {
        final int m1 = r.length, n1 = c.length;
        final Matrix M = new Matrix(m1, n1);
        for (int i = 0; i < m1; i++)
        {
            for (int j = 0; j < n1; j++)
            {
                M.A[i][j] = this.A[r[i]][c[j]];
            }
        }
        return M;
    }

    /**
     * Get row dimension.
     * 
     * @return m, the number of rows.
     */
    public int getRowDimension()
    {
        return this.m;
    }

    @Override
    public int hashCode()
    {
        final int prime = 31;
        int result = 1;
        result = (prime * result) + Arrays.hashCode(this.A);
        result = (prime * result) + this.m;
        result = (prime * result) + this.n;
        return result;
    }

    /**
     * Makes the identity matrix
     */
    public void identity()
    {
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                this.A[i][j] = (i == j ? 1D : 0D);
            }
        }
    }

    /**
     * Matrix inverse or pseudoinverse
     * 
     * @return inverse(A) if A is square, pseudoinverse otherwise.
     * @see #solve
     */
    public Matrix inverse()
    {
        final Matrix I = new Matrix(this.m);
        I.identity();
        return this.solve(I);
    }

    /**
     * LU Decomposition
     * 
     * @return LUDecomposition
     * @see LUDecomposition
     */
    public LUDecomposition lu()
    {
        return new LUDecomposition(this);
    }

    /**
     * C = A - B
     * 
     * @param B
     *            another matrix
     * @return new Matrix A - B
     * @see #minusEquals
     */
    public Matrix minus(final Matrix B)
    {
        this.checkMatrixDimensions(B);
        final Matrix M = new Matrix(this.m, this.n);
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                M.A[i][j] = this.A[i][j] - B.A[i][j];
            }
        }
        return M;
    }

    /**
     * A = A - B
     * 
     * @param B
     *            another matrix
     * @see #minus
     */
    public void minusEquals(final Matrix B)
    {
        this.checkMatrixDimensions(B);
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                this.A[i][j] -= B.A[i][j];
            }
        }
    }

    /**
     * One norm
     * 
     * @return maximum column sum.
     */
    public double norm1()
    {
        double f = 0;
        for (int j = 0; j < this.n; j++)
        {
            double s = 0D;
            for (int i = 0; i < this.m; i++)
            {
                s += abs(this.A[i][j]);
            }
            if (f < s)
            {
                f = s;
            }
        }
        return f;
    }

    /**
     * Two norm
     * 
     * @return maximum singular value.
     */
    public double norm2()
    {
        return this.svd().norm2();
    }

    /**
     * Frobenius norm
     * 
     * @return sqrt of sum of squares of all elements.
     */
    public double normF()
    {
        double f = 0;
        for (final double[] row : this.A)
        {
            for (final double d : row)
            {
                f = hypot(f, d);
            }
        }
        return f;
    }

    /**
     * Infinity norm
     * 
     * @return maximum row sum.
     */
    public double normInf()
    {
        double f = 0;
        for (final double[] row : this.A)
        {
            double s = 0;
            for (final double d : row)
            {
                s += abs(d);
            }
            if (f < s)
            {
                f = s;
            }
        }
        return f;
    }

    /**
     * C = A + B
     * 
     * @param B
     *            another matrix
     * @return new Matrix A + B
     * @see #plusEquals
     */
    public Matrix plus(final Matrix B)
    {
        this.checkMatrixDimensions(B);
        final Matrix M = new Matrix(this.m, this.n);
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                M.A[i][j] = this.A[i][j] + B.A[i][j];
            }
        }
        return M;
    }

    /**
     * A = A + B
     * 
     * @param B
     *            another matrix
     * @see #plus
     */
    public void plusEquals(final Matrix B)
    {
        this.checkMatrixDimensions(B);
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                this.A[i][j] += B.A[i][j];
            }
        }
    }

    /**
     * Print the matrix to stdout. Line the elements up in columns. Use the
     * format object, and right justify within columns of width characters. Note
     * that is the matrix is to be read back in, you probably will want to use a
     * NumberFormat that is set to US Locale.
     * 
     * @param format
     *            A Formatting object for individual elements.
     * @param width
     *            Field width for each column.
     * @see java.text.DecimalFormat#setDecimalFormatSymbols
     */

    public void print(final NumberFormat format, final int width)
    {
        this.print(new PrintWriter(System.out, true), format, width);
    }

    /**
     * Print the matrix to the output stream. Line the elements up in columns.
     * Use the format object, and right justify within columns of width
     * characters. Note that is the matrix is to be read back in, you probably
     * will want to use a NumberFormat that is set to US Locale.
     * 
     * @param output
     *            the output stream.
     * @param format
     *            A formatting object to format the matrix elements
     * @param width
     *            Column width.
     * @see java.text.DecimalFormat#setDecimalFormatSymbols
     */
    public void print(final PrintWriter output, final NumberFormat format, final int width)
    {
        // start on new line.
        output.println();

        for (final double row[] : this.A)
        {
            for (final double d : row)
            {
                // format the number
                final String s = format.format(d);
                // At _least_ 1 space
                final int padding = (width > s.length()) ? width - s.length() : 1;
                for (int k = 0; k < padding; k++)
                {
                    output.print(' ');
                }
                output.print(s);
            }
            output.println();
        }

        // end with blank line.
        output.println();
    }

    /**
     * Print the matrix to the output stream. Line the elements up in columns
     * with a Fortran-like 'Fw.d' style format.
     * 
     * @param output
     *            Output stream.
     * @param w
     *            Column width.
     * @param d
     *            Number of digits after the decimal.
     */
    public void print(final PrintWriter output, final int w, final int d)
    {
        // DecimalFormat is a little disappointing coming from Fortran or C's printf.
        // Since it doesn't pad on the left, the elements will come out different
        // widths.  Consequently, we'll pass the desired column width in as an
        // argument and do the extra padding ourselves.
        final DecimalFormat format = new DecimalFormat();
        format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
        format.setMinimumIntegerDigits(1);
        format.setMaximumFractionDigits(d);
        format.setMinimumFractionDigits(d);
        format.setGroupingUsed(false);
        this.print(output, format, w + 2);
    }

    /**
     * Print the matrix to stdout. Line the elements up in columns with a
     * Fortran-like 'Fw.d' style format.
     * 
     * @param w
     *            Column width.
     * @param d
     *            Number of digits after the decimal.
     */
    public void print(final int w, final int d)
    {
        this.print(new PrintWriter(System.out, true), w, d);
    }

    /**
     * QR Decomposition
     * 
     * @return QRDecomposition
     * @see QRDecomposition
     */
    public QRDecomposition qr()
    {
        return new QRDecomposition(this);
    }

    /**
     * Creates a random Matrix
     */
    public void random()
    {
        final Random rnd = new Random();
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                this.A[i][j] = rnd.nextDouble();
            }
        }
    }

    /**
     * Creates a random Matrix
     */
    public void randomInt()
    {
        final Random rnd = new Random();
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                this.A[i][j] = rnd.nextInt();
            }
        }
    }

    /**
     * Matrix rank
     * 
     * @return effective numerical rank, obtained from SVD.
     * @see SingularValueDecomposition#rank
     */
    public int rank()
    {
        return this.svd().rank();
    }

    /**
     * 
     * @param s
     */
    public void set(final double s)
    {
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                this.A[i][j] = s;
            }
        }
    }

    /**
     * 
     * @param i
     *            index
     * @param s
     */
    public void set(final int i, final double s)
    {
        this.A[i / this.m][i % this.n] = s;
    }

    /**
     * Set a single element.
     * 
     * @param i
     *            Row index.
     * @param j
     *            Column index.
     * @param s
     *            A(i,j).
     */
    public void set(final int i, final int j, final double s)
    {
        this.A[i][j] = s;
    }

    /**
     * Set a submatrix.
     * 
     * @param i0
     *            Initial row index
     * @param i1
     *            Final row index
     * @param j0
     *            Initial column index
     * @param j1
     *            Final column index
     * @param X
     *            A(i0:i1,j0:j1)
     */
    public void setMatrix(final int i0, final int i1, final int j0, final int j1, final Matrix X)
    {
        for (int i = i0; i <= i1; i++)
        {
            for (int j = j0; j <= j1; j++)
            {
                this.A[i][j] = X.A[i - i0][j - j0];
            }
        }
    }

    /**
     * Set a submatrix.
     * 
     * @param i0
     *            Initial row index
     * @param i1
     *            Final row index
     * @param c
     *            Array of column indices.
     * @param X
     *            A(i0:i1,c(:))
     */
    public void setMatrix(final int i0, final int i1, final int c[], final Matrix X)
    {
        for (int i = i0; i <= i1; i++)
        {
            for (int j = 0; j < c.length; j++)
            {
                this.A[i][c[j]] = X.A[i - i0][j];
            }
        }
    }

    /**
     * Set a submatrix.
     * 
     * @param r
     *            Array of row indices.
     * @param j0
     *            Initial column index
     * @param j1
     *            Final column index
     * @param X
     *            A(r(:),j0:j1)
     */
    public void setMatrix(final int r[], final int j0, final int j1, final Matrix X)
    {
        for (int i = 0; i < r.length; i++)
        {
            for (int j = j0; j <= j1; j++)
            {
                this.A[r[i]][j] = X.A[i][j - j0];
            }
        }
    }

    /**
     * Set a submatrix.
     * 
     * @param r
     *            Array of row indices.
     * @param c
     *            Array of column indices.
     * @param X
     *            A(r(:),c(:))
     */
    public void setMatrix(final int r[], final int c[], final Matrix X)
    {
        for (int i = 0; i < r.length; i++)
        {
            for (int j = 0; j < c.length; j++)
            {
                this.A[r[i]][c[j]] = X.get(i, j);
            }
        }
    }

    /**
     * Solve A*X = B
     * 
     * @param B
     *            right hand side
     * @return solution if A is square, least squares solution otherwise.
     *         Returns null if no such solution exists (matrix is singular or
     *         rank deficient).
     */
    public Matrix solve(final Matrix B)
    {
        return (this.m == this.n ? this.lu().solve(B) : this.qr().solve(B));
    }

    /**
     * Solve X*A = B, which is also A'*X' = B'
     * 
     * @param B
     *            right hand side
     * @return solution if A is square, least squares solution otherwise.
     *         Returns null if no such solution exists (matrix is singular or
     *         rank deficient).
     */
    public Matrix solveTranspose(final Matrix B)
    {
        return this.transpose().solve(B.transpose());
    }

    /**
     * Singular Value Decomposition
     * 
     * @return SingularValueDecomposition
     * @see SingularValueDecomposition
     */
    public SingularValueDecomposition svd()
    {
        return new SingularValueDecomposition(this);
    }

    /**
     * Multiply a matrix by a scalar, C = s*A.
     * 
     * @param s
     *            scalar
     * @return new Matrix s*A
     */
    public Matrix times(final double s)
    {
        final Matrix M = new Matrix(this.m, this.n);
        for (int i = 0; i < this.m; i++)
        {
            for (int j = 0; j < this.n; j++)
            {
                M.A[i][j] = s * this.A[i][j];
            }
        }
        return M;
    }

    /**
     * Linear algebraic matrix multiplication, A * B
     * 
     * @param B
     *            another matrix
     * @return Matrix product, A * B, null if matrix dimensions don't agree.
     */
    public Matrix times(final Matrix B)
    {
        if (B.m != this.n)
        {
            return null;
            //throw new IllegalArgumentException("Matrix inner dimensions must agree."); //$NON-NLS-1$
        }
        else
        {
            final Matrix X = new Matrix(this.m, B.n);
            for (int j = 0; j < B.n; j++)
            {
                for (int i = 0; i < this.m; i++)
                {
                    double s = 0;
                    for (int k = 0; k < this.n; k++)
                    {
                        s += this.A[i][k] * B.A[k][j];
                    }
                    X.A[i][j] = s;
                }
            }
            return X;
        }
    }

    /**
     * Multiply a matrix by a scalar in place, A = s*A
     * 
     * @param s
     *            scalar
     */
    public void timesEquals(final double s)
    {
        for (int i = 0; i < this.m; i++)
        {
            for (int j = 0; j < this.n; j++)
            {
                this.A[i][j] *= s;
            }
        }
    }

    @Override
    public String toString()
    {
        final StringBuffer st = new StringBuffer();
        for (final double row[] : this.A)
        {
            for (final double d : row)
            {
                st.append(d);
                st.append(" ");
            }
            st.append("\n");
        }
        return st.toString();
    }

    /**
     * Matrix trace.
     * 
     * @return sum of the diagonal elements.
     */
    public double trace()
    {
        final int dim = ((this.m < this.n) ? this.m : this.n);
        double t = 0;
        for (int i = 0; i < dim; i++)
        {
            t += this.A[i][i];
        }
        return t;
    }

    /**
     * Matrix transpose.
     * 
     * @return new Matrix A'
     */
    public Matrix transpose()
    {
        final Matrix M = new Matrix(this.n, this.m);
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                M.A[j][i] = this.A[i][j];
            }
        }
        return M;
    }

    /**
     * Transposes this matrix.
     */
    public void transposeThis()
    {
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < i; j++)
            {
                final double t = this.A[i][j];
                this.A[i][j] = this.A[j][i];
                this.A[j][i] = t;
            }
        }
    }

    /**
     * Unary minus
     * 
     * @return new Matrix -A
     */
    public Matrix uminus()
    {
        final Matrix M = new Matrix(this.m, this.n);
        for (int i = 0; i < this.A.length; i++)
        {
            for (int j = 0; j < this.A[i].length; j++)
            {
                M.A[i][j] = -this.A[i][j];
            }
        }
        return M;
    }
}
