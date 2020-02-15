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
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoublePredicate;
import java.util.function.DoubleUnaryOperator;

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
 * <P>
 * Note that this class is <b>mutable</b> and none of it's methods are synchronized.
 * Be careful using it from different threads.
 * </P>
 *
 * @author The MathWorks, Inc. and the National Institute of Standards and
 *         Technology.
 * @version 2.0
 * @see <a href="http://tweimer.github.io/java-matrix/">java-matrix</a>
 */
public class Matrix extends RectangularMatrix implements Cloneable, Serializable {
    /**
     * For the Serializeable interface.
     */
    private static final long serialVersionUID = 1;

    /**
     * Construct a matrix from a copy of a 2-D array.
     *
     * @param A
     *            Two-dimensional array of doubles.
     * @return Matrix with copied array
     * @throws IllegalArgumentExce
     * ption
     *                All rows must have the same length
     * @throws NullPointerException
     *                Iff A or any sub-array is {@code null}.
     */
    public static Matrix constructWithCopy(final double[][] A) {
        final var nRows = A.length;
        final var nColumns = A[0].length;
        final var C = new double[nRows][];
        for (int i = 0; i < C.length; i++) {
            if (A[i].length != nColumns) {
                throw new IllegalArgumentException("All rows must have the same length."); //$NON-NLS-1$
            } else {
                C[i] = Arrays.copyOf(A[i], nColumns);
            }
        }
        return new Matrix(nRows, nColumns, C);
    }

    /**
     * Creates a diagonal Matrix.
     *
     * @param d
     *            diagonal elements
     * @return diagonal Matrix
     * @throws NullPointerException
     *    Iff {@code d == null}
     */
    public static Matrix diag(final double... d) {
        final var n = d.length;
        final var diag = new double[n][n];
        for (var i = 0; i < n; i++) {
            diag[i][i] = d[i];
        }
        return new Matrix(n, diag);
    }

    /**
     * Generate identity matrix.
     *
     * @param m
     *            Number of rows and colums.
     * @return An m-by-n square matrix with ones on the diagonal and zeros elsewhere.
     */
    public static Matrix identity(final int m) {
        return Matrix.identity(m, m);
    }

    /**
     * Generate identity matrix.
     *
     * @param nRows
     *            Number of rows.
     * @param nColumns
     *            Number of colums.
     * @return An m-by-n matrix with ones on the diagonal and zeros elsewhere.
     */
    public static Matrix identity(final int nRows, final int nColumns) {
        final var I = new Matrix(nRows, nColumns);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                I.A[r][c] = r == c ? 1D : 0D;
            }
        }
        return I;
    }

    /**
     * Generate matrix with random elements.
     *
     * @param n
     *            Number of rows and colums.
     * @return An n-by-n matrix with uniformly distributed random elements.
     */
    public static Matrix random(final int n) {
        return Matrix.random(n, n);
    }

    /**
     * Generate matrix with random elements.
     *
     * @param nRows
     *            Number of rows.
     * @param nColumns
     *            Number of colums.
     * @return An m-by-n matrix with uniformly distributed random elements.
     */
    public static Matrix random(final int nRows, final int nColumns) {
        final var R = new Matrix(nRows, nColumns);
        final var rnd = new Random();
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                R.A[r][c] = rnd.nextDouble();
            }
        }
        return R;
    }

    /**
     * Generate matrix with random elements.
     *
     * @param n
     *            Number of rows and colums.
     * @return An n-by-n matrix with uniformly distributed random elements.
     */
    public static Matrix randomInt(final int n) {
        return Matrix.randomInt(n, n);
    }

    /**
     * Generate matrix with random elements.
     *
     * @param nRows
     *            Number of rows.
     * @param nColumns
     *            Number of colums.
     * @return An m-by-n matrix with uniformly distributed random elements.
     */
    public static Matrix randomInt(final int nRows, final int nColumns) {
        final var M = new Matrix(nRows, nColumns);
        final var rnd = new Random();
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                M.A[r][c] = rnd.nextInt();
            }
        }
        return M;
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
     *             if any IOException occurs while reading
     */
    public static Matrix read(final BufferedReader input) throws IOException {
        final var tokenizer = new StreamTokenizer(input);

        // Although StreamTokenizer will parse numbers, it doesn't recognize
        // scientific notation (E or D); however, Double.valueOf does.
        // The strategy here is to disable StreamTokenizer's number parsing.
        // We'll only get whitespace delimited words, EOL's and EOF's.
        // These words should all be numbers, for Double.valueOf to parse.

        tokenizer.resetSyntax();
        tokenizer.wordChars(0, 255);
        tokenizer.whitespaceChars(0, ' ');
        tokenizer.eolIsSignificant(true);

        while (tokenizer.nextToken() == StreamTokenizer.TT_EOL) {
            // Ignore initial empty lines
        }

        if (tokenizer.ttype == StreamTokenizer.TT_EOF) {
            throw new IOException("Unexpected EOF on matrix read."); //$NON-NLS-1$
        }

        final var vD = new Vector<Double>();
        do {
            // Read & store 1st row.
            vD.addElement(Double.valueOf(tokenizer.sval));
        } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);

        // Now we've got the number of columns!
        final var n = vD.size();
        var row = new double[n];
        var j = 0;
        for (final var d : vD) {
            // extract the elements of the 1st row.
            row[j++] = d.doubleValue();
        }

        // Start storing rows instead of columns.
        final var v = new Vector<double[]>();
        v.addElement(row);
        while (tokenizer.nextToken() == StreamTokenizer.TT_WORD) {
            // While non-empty lines
            v.addElement(row = new double[n]);
            j = 0;
            do {
                if (j >= n) {
                    throw new IOException("Row " + v.size() + " is too long."); //$NON-NLS-1$ //$NON-NLS-2$
                }
                row[j++] = Double.valueOf(tokenizer.sval).doubleValue();
            } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
            if (j < n) {
                throw new IOException("Row " + v.size() + " is too short."); //$NON-NLS-1$ //$NON-NLS-2$
            }
        }

        // Now we've got the number of rows.
        final var m = v.size();
        final var A = new double[m][];

        // copy the rows out of the vector
        v.copyInto(A);
        return new Matrix(A);
    }

    /**
     * Array for internal storage of elements.
     *
     * @serial internal array storage.
     */
    private final double[][] A;

    /**
     * Construct a matrix from a one-dimensional packed array.
     *
     * @param vals
     *            One-dimensional array of doubles, packed by columns (ala
     *            Fortran).
     * @param nRows
     *            Number of rows.
     * @throws IllegalArgumentException
     *             Array length must be a multiple of m.
     */
    public Matrix(final double[] vals, final int nRows) {
        this(nRows, nRows != 0 ? vals.length / nRows : 0);
        if (nRows * nColumns != vals.length) {
            throw new IllegalArgumentException("Array length must be a multiple of m."); //$NON-NLS-1$
        } else {
            for (var i = 0; i < nRows; i++) {
                for (var j = 0; j < nColumns; j++) {
                    this.A[i][j] = vals[i + j * nRows];
                }
            }
        }
    }

    /**
     * Construct a matrix from a 2-D array.
     *
     * @param A
     *            Two-dimensional array of doubles. This does <b>not</b> copy
     *            the given array, it just stores it's reference.
     *            You should avoid modifying that array afterwards,
     *            as this matrix is backed by that array
     * @throws IllegalArgumentException
     *                All rows must have the same length
     * @see #constructWithCopy 
     */
    public Matrix(final double[][] A) {
        // TODO: how to handle null here?
        this(A.length, A[0].length, A);

        // check if each row has the same length
        for (var r : this.A) {
            if (r.length != nColumns) {
                throw new IllegalArgumentException("All rows must have the same length."); //$NON-NLS-1$
            }
        }
    }

    /**
     * Construct a matrix quickly without checking arguments.
     *
     * @param nRows
     *            Number of rows.
     * @param nColumns
     *            Number of colums.
     * @param A
     *            Two-dimensional array of doubles. This does <b>not</b> copy
     *            the given array, it just stores it's reference.
     *            You should avoid modifying that array afterwards,
     *            as this matrix is backed by that array
     */
    public Matrix(final int nRows, final int nColumns, final double[][] A) {
        super(nRows, nColumns);
        this.A = A;
    }

    /**
     * Construct a matrix quickly without checking arguments.
     *
     * @param n
     *            Number of rows and colums.
     * @param A
     *            Two-dimensional array of doubles.
     */
    public Matrix(final int n, final double[][] A) {
        this(n, n, A);
    }

    /**
     * Construct an n-by-n matrix of zeros.
     *
     * @param n
     *            Number of rows and colums.
     */
    public Matrix(final int n) {
        this(n, n);
    }

    /**
     * Construct an m-by-n matrix of zeros.
     *
     * @param nRows
     *            Number of rows.
     * @param nColumns
     *            Number of colums
     */
    public Matrix(final int nRows, final int nColumns) {
        super(nRows, nColumns);
        this.A = new double[nRows][nColumns];
    }

    /**
     * Construct an m-by-n constant matrix.
     *
     * @param nRows
     *            Number of rows.
     * @param nColumns
     *            Number of colums.
     * @param s
     *            Fill the matrix with this scalar value.
     */
    public Matrix(final int nRows, final int nColumns, final double s) {
        this(nRows, nColumns);
        this.fill(s);
    }

    /**
     * Copies Matrix X.
     *
     * @param X
     *            Matrix to be copied
     */
    public Matrix(final Matrix X) {
        this(X.nRows, X.nColumns, X.getArrayCopy());
    }
    

   /**
    * This is used to build a Matrix with a functional interface.
    *
    * @param nRows
    *            Number of rows.
    * @param nColumns
    *            Number of colums.
    * @param matrix
    *            A lambda expression that describes the given Matrix.
    *            It must return valid numbers within the given number of rows and columns.
    *            Any {@link RuntimeException} in {@link FunctionalMatrix#get(int, int)} is thrown back to the caller.
    */
    public Matrix(final int nRows, final int nColumns, final FunctionalMatrix matrix) {
        this(nRows, nColumns);
        for (var i = 0; i < nRows; i++) {
            for (var j = 0; j < nColumns; j++) {
                this.A[i][j] = matrix.get(i, j);
            }
        }
    }

    /**
     * This is used to build a Matrix with a functional interface.
     *
     * @param nRows
     *            Number of rows.
     * @param nColumns
     *            Number of colums.
     * @param matrix
     *            A lambda expression that describes the given Matrix.
     *            It must return valid numbers within the given number of rows and columns.
     *            Any {@link RuntimeException} in {@link FunctionalMatrix#get(int, int)} is thrown back to the caller.
     * @param f This is applied after parameter matrix is evaluated.
     */
     public Matrix(final int nRows, final int nColumns,
             final FunctionalMatrix matrix, final DoubleUnaryOperator f) {
         this(nRows, nColumns);
         for (var i = 0; i < nRows; i++) {
             for (var j = 0; j < nColumns; j++) {
                 this.A[i][j] = f.applyAsDouble(matrix.get(i, j));
             }
         }
     }

    /**
     * Returns true if all elements of the Matrix match the provided predicate.
     * This uses lazy evaluation.
     *
     * @param predicate
     *            A predicate to test all elements of the Matrix
     * @return True if all elements match the given predicate, false otherwise.
     * @throws NullPointerException
     *             iff predicate == null
     * @see #anyMatch(DoublePredicate)
     */
    public boolean allMatch(final DoublePredicate predicate) {
        for (final var row : this.A) {
            for (final var d : row) {
                if (!predicate.test(d)) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Returns true if any element of the Matrix matches the provided
     * predicate. This uses lazy evaluation.
     *
     * @param predicate
     *            A predicate to test all elements of the Matrix
     * @return True if any element matches the given predicate, false otherwise.
     * @throws NullPointerException
     *             iff predicate == null
     * @see #allMatch(DoublePredicate)
     */
    public boolean anyMatch(final DoublePredicate predicate) {
        for (final var row : this.A) {
            for (final var d : row) {
                if (predicate.test(d)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Element-by-element left division, C = A.\B.
     *
     * @param B
     *            another matrix
     * @return A.\B
     */
    public Matrix arrayLeftDivide(final Matrix B) {
        this.checkMatrixDimensions(B);
        final var M = new Matrix(nRows, nColumns);
        for (var r = 0; r < nRows; r++) {
            for (int c = 0; c < nColumns; c++) {
                M.A[r][c] = B.A[r][c] / this.A[r][c];
            }
        }
        return M;
    }

    /**
     * Element-by-element left division in place, A = A.\B.
     *
     * @param B
     *            another matrix
     */
    public void arrayLeftDivideEquals(final Matrix B) {
        this.checkMatrixDimensions(B);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                this.A[r][c] = B.A[r][c] / this.A[r][c];
            }
        }
    }

    /**
     * Element-by-element right division, C = A./B.
     *
     * @param B
     *            another matrix
     * @return A./B
     */
    public Matrix arrayRightDivide(final Matrix B) {
        this.checkMatrixDimensions(B);
        final var M = new Matrix(nRows, nColumns);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                M.A[r][c] = this.A[r][c] / B.A[r][c];
            }
        }
        return M;
    }

    /**
     * Element-by-element right division in place, A = A./B.
     *
     * @param B
     *            another matrix
     */
    public void arrayRightDivideEquals(final Matrix B) {
        this.checkMatrixDimensions(B);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                this.A[r][c] /= B.A[r][c];
            }
        }
    }

    /**
     * Element-by-element multiplication, C = A.*B.
     *
     * @param B
     *            another matrix
     * @return A.*B
     * @see #arrayTimesEquals
     */
    public Matrix arrayTimes(final Matrix B) {
        this.checkMatrixDimensions(B);
        final var M = new Matrix(nRows, nColumns);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                M.A[r][c] = this.A[r][c] * B.A[r][c];
            }
        }
        return M;
    }

    /**
     * Element-by-element multiplication in place, A = A.*B.
     *
     * @param B
     *            another matrix
     * @see #arrayTimes
     */
    public void arrayTimesEquals(final Matrix B) {
        this.checkMatrixDimensions(B);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                this.A[r][c] *= B.A[r][c];
            }
        }
    }

    /**
     * Check if size(A) == size(B).
     *
     * @param B another Matrix
     * @throws IllegalArgumentException
     *             if they don't agree.
     */
    private void checkMatrixDimensions(final Matrix B) {
        if (!this.equalDimensions(B)) {
            throw new IllegalArgumentException("Matrix dimensions must agree."); //$NON-NLS-1$
        }
    }

    /**
     * Compares dimension of both matrices.
     *
     * @param B
     *            another Matrix
     * @return true if both have the same dimension, false otherwise, returns
     *         false if B==null.
     */
    public boolean equalDimensions(final Matrix B) {
        return B != null && B.nRows == nRows && B.nColumns == nColumns;
    }

    /**
     * Cholesky Decomposition.
     *
     * @return CholeskyDecomposition
     * @see CholeskyDecomposition
     */
    public CholeskyDecomposition chol() {
        return new CholeskyDecomposition(this);
    }

    /**
     * Matrix condition (2 norm).
     *
     * @return ratio of largest to smallest singular value.
     * @see SingularValueDecomposition#cond()
     */
    public double cond() {
        return this.svd().cond();
    }

    /**
     * Returns a copy of the matrix.
     *
     * @return copy of this Matrix
     * @see #Matrix(Matrix)
     */
    @Override
    public Matrix clone() {
        return new Matrix(this);
    }

    /**
     * Matrix determinant.
     *
     * @return determinant
     */
    public double det() {
        return this.lu().det();
    }

    /**
     * Eigenvalue Decomposition.
     *
     * @return EigenvalueDecomposition
     * @see EigenvalueDecomposition
     */
    public EigenvalueDecomposition eig() {
        return new EigenvalueDecomposition(this);
    }

    /**
     * Compares a Matrix to another Matrix.
     *
     * @param other
     *            another Matrix
     * @return true if other equals A
     */
    public boolean equals(final Matrix other) {
        return (other == this)
                || (other != null && nRows == other.nRows && nColumns == other.nColumns && Arrays.deepEquals(this.A, other.A));
    }

    /**
     * Compares a Matrix to another Matrix.
     *
     * @param obj
     *            another Object
     * @return true if obj is a Matrix and equals A
     */
    @Override
    public boolean equals(final Object obj) {
        return (obj instanceof Matrix) && this.equals((Matrix) obj);
    }

    /**
     * Returns index at index i, row-by-row (MATLAB like).
     *
     * @param i
     *            index
     * @return entry
     */
    public double get(final int i) {
        return this.A[i / nRows][i % nColumns];
    }

    /**
     * Get a single element.
     *
     * @param r
     *            Row index.
     * @param c
     *            Column index.
     * @return A(i,j)
     * @throws ArrayIndexOutOfBoundsException
     *             if indices are invalid
     */
    @Override
    public double get(final int r, final int c) {
        return this.A[r][c];
    }

    /**
     * Access the internal two-dimensional array.
     * You should avoid modifying the returned array,
     * as this matrix is backed by that array
     *
     * @return Pointer to the two-dimensional array of matrix elements.
     */
    public double[][] getArray() {
        return this.A;
    }

    /**
     * Copy the internal two-dimensional array.
     *
     * @return Two-dimensional array copy of matrix elements.
     */
    @Override
    public double[][] getArrayCopy() {
        // Copying the array itsself is faster then the default implementation.
        final var C = new double[nRows][];
        for (var r = 0; r < nRows; r++) {
            C[r] = Arrays.copyOf(this.A[r], nColumns);
        }
        return C;
    }

    /**
     * Make a one-dimensional column packed copy of the internal array.
     *
     * @return Matrix elements packed in a one-dimensional array by columns.
     */
    public double[] getColumnPackedCopy() {
        final var vals = new double[nRows * nColumns];
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                vals[r + c * nRows] = this.A[r][c];
            }
        }
        return vals;
    }

    /**
     * Make a one-dimensional row packed copy of the internal array.
     *
     * @return Matrix elements packed in a one-dimensional array by rows.
     */
    public double[] getRowPackedCopy() {
        final var vals = new double[nRows * nColumns];
        var i = 0;
        for (var row : this.A) {
            for (var d : row) {
                vals[i++] = d;
            }
        }
        return vals;
    }

    /**
     * Get a submatrix.
     *
     * @param r0
     *            Initial row index (inclusive)
     * @param r1
     *            Final row index (inclusive)
     * @param c
     *            Array of column indices.
     * @return A(r0:r1,c(:))
     * @throws ArrayIndexOutOfBoundsException
     *             Submatrix indices
     */
    public Matrix getMatrix(final int r0, final int r1, final int[] c) {
        final var m1 = r1 - r0 + 1;
        final var n1 = c.length;
        final var M = new Matrix(m1, n1);
        for (var i = r0; i <= r1; i++) {
            for (var j = 0; j < n1; j++) {
                M.A[i - r0][j] = this.A[i][c[j]];
            }
        }
        return M;
    }

    /**
     * Get a submatrix.
     *
     * @param r0
     *            Initial row index (inclusive)
     * @param r1
     *            Final row index (inclusive)
     * @param c0
     *            Initial column index (inclusive)
     * @param c1
     *            Final column index (inclusive)
     * @return A(r0:r1,c0:c1)
     * @exception ArrayIndexOutOfBoundsException
     *                Submatrix indices
     */
    public Matrix getMatrix(final int r0, final int r1, final int c0, final int c1) {
        final var m1 = r1 - r0 + 1;
        final var n1 = c1 - c0 + 1;
        final var M = new Matrix(m1, n1);
        for (var r = r0; r <= r1; r++) {
            for (var c = c0; c <= c1; c++) {
                M.A[r - r0][c - c0] = this.A[r][c];
            }
        }
        return M;
    }

    /**
     * Get a submatrix.
     *
     * @param r
     *            Array of row indices
     * @param c0
     *            Initial column index (inclusive)
     * @param c1
     *            Final column index (inclusive)
     * @return A(r(:),c0:c1)
     * @exception ArrayIndexOutOfBoundsException
     *                Submatrix indices
     */
    public Matrix getMatrix(final int r[], final int c0, final int c1) {
        final var m1 = r.length;
        final var n1 = c1 - c0 + 1;
        final var M = new Matrix(m1, n1);
        for (int i = 0; i < m1; i++) {
            for (int c = c0; c <= c1; c++) {
                M.A[i][c - c0] = this.A[r[i]][c];
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
    public Matrix getMatrix(final int r[], final int c[]) {
        final var m1 = r.length;
        final var n1 = c.length;
        final var A = new double[m1][n1];
        for (var i = 0; i < m1; i++) {
            for (var j = 0; j < n1; j++) {
                A[i][j] = this.A[r[i]][c[j]];
            }
        }
        return new Matrix(m1, n1, A);
    }

    /**
     * Returns a deep hash code for this Matrix.
     * @return hash code
     */
    @Override
    public int hashCode() {
        return Arrays.deepHashCode(this.A);
    }


    /**
     * Matrix inverse or pseudoinverse.
     *
     * @return inverse(A) if A is square, pseudoinverse otherwise.
     * @see #solve
     */
    public Matrix inverse() {
        return solve(identity(nRows));
    }

    /**
     * LU Decomposition.
     *
     * @return LUDecomposition
     * @see LUDecomposition
     */
    public LUDecomposition lu() {
        return new LUDecomposition(this);
    }

    /**
     * C = A - B.
     *
     * @param B
     *            another matrix
     * @return new Matrix A - B
     * @see #minusEquals
     */
    public Matrix minus(final Matrix B) {
        this.checkMatrixDimensions(B);
        final var M = new Matrix(nRows, nColumns);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                M.A[r][c] = this.A[r][c] - B.A[r][c];
            }
        }
        return M;
    }

    /**
     * A = A - B.
     *
     * @param B
     *            another matrix
     * @see #minus
     */
    public void minusEquals(final Matrix B) {
        this.checkMatrixDimensions(B);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                this.A[r][c] -= B.A[r][c];
            }
        }
    }

    /**
     * One norm.
     *
     * @return maximum column sum.
     */
    public double norm1() {
        var norm1 = 0D;
        for (var c = 0; c < nColumns; c++) {
            var colSum = 0D;
            for (var r = 0; r < nRows; r++) {
                colSum += abs(this.A[r][c]);
            }
            if (norm1 < colSum) {
                norm1 = colSum;
            }
        }
        return norm1;
    }

    /**
     * Two norm.
     *
     * @return maximum singular value.
     */
    public double norm2() {
        return this.svd().norm2();
    }

    /**
     * Frobenius norm.
     *
     * @return sqrt of sum of squares of all elements.
     */
    public double normF() {
        double normF = 0;
        for (final double[] row : this.A) {
            for (final double d : row) {
                normF = hypot(normF, d);
            }
        }
        return normF;
    }

    /**
     * Infinity norm.
     *
     * @return maximum row sum.
     */
    public double normInf() {
        var normInf = 0D;
        for (final var row : this.A) {
            var rowSum = 0D;
            for (final var d : row) {
                rowSum += abs(d);
            }
            if (normInf < rowSum) {
                normInf = rowSum;
            }
        }
        return normInf;
    }

    /**
     * C = A + B.
     *
     * @param B
     *            another matrix
     * @return new Matrix A + B
     * @throws IllegalArgumentException
     *     Iff matrix dimensions don't agree
     * @see #plusEquals
     */
    public Matrix plus(final Matrix B) {
        this.checkMatrixDimensions(B);
        // urgh, casting needed, to avoid a recursion
        // cast safe, since Matrix is subtype of FunctionalMatrix
        return plus((FunctionalMatrix)B);
    }

    /**
     * This adds all elements of the argument to this Matrix,
     * without doing any bounds checks.
     * 
     * @param m
     *   Another functional Matrix
     * @return
     *  Returns a new matrix of same size as this,
     *  with the result of the addition.
     */
    public Matrix plus(final FunctionalMatrix m) {
        final var M = new Matrix(nRows, nColumns);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                M.A[r][c] = this.A[r][c] + m.get(r, c);
            }
        }
        return M;
    }

    /**
     * A = A + B.
     *
     * @param B
     *            another matrix
     * @see #plus
     */
    public void plusEquals(final Matrix B) {
        this.checkMatrixDimensions(B);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                this.A[r][c] += B.A[r][c];
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

    public void print(final NumberFormat format, final int width) {
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
    public void print(final PrintWriter output, final NumberFormat format, final int width) {
        // start on new line.
        output.println();

        for (final var row : this.A) {
            for (final var d : row) {
                // format the number
                final var s = format.format(d);
                // At _least_ 1 space
                final var padding = width > s.length() ? width - s.length() : 1;
                for (var k = 0; k < padding; k++) {
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
    public void print(final PrintWriter output, final int w, final int d) {
        // DecimalFormat is a little disappointing coming from Fortran or C's printf.
        // Since it doesn't pad on the left, the elements will come out different
        // widths. Consequently, we'll pass the desired column width in as an
        // argument and do the extra padding ourselves.
        final var format = new DecimalFormat();
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
    public void print(final int w, final int d) {
        this.print(new PrintWriter(System.out, true), w, d);
    }

    /**
     * QR Decomposition.
     *
     * @return QRDecomposition
     * @see QRDecomposition
     */
    public QRDecomposition qr() {
        return new QRDecomposition(this);
    }


    /**
     * Matrix rank.
     *
     * @return effective numerical rank, obtained from SVD.
     * @see SingularValueDecomposition#rank
     */
    public int rank() {
        return svd().rank();
    }

    /**
     * Sets all values to s.
     *
     * @param s
     *            scalar
     */
    public void fill(final double s) {
        for (var row : A) {
            Arrays.fill(row, s);
        }
    }

    /**
     * Sets all values to s.
     *
     * @param i
     *            index (counting by columns)
     * @param s
     *            scalar
     */
    public void set(final int i, final double s) {
        this.A[i / nRows][i % nColumns] = s;
    }

    /**
     * Set a single element.
     *
     * @param r
     *            Row index.
     * @param c
     *            Column index.
     * @param s
     *            A(i,j).
     */
    public void set(final int r, final int c, final double s) {
        this.A[r][c] = s;
    }

    /**
     * Set a submatrix.
     *
     * @param r0
     *            Initial row index
     * @param r1
     *            Final row index
     * @param c0
     *            Initial column index
     * @param c1
     *            Final column index
     * @param X
     *            A(r0:r1,c0:c1)
     */
    public void setMatrix(final int r0, final int r1, final int c0, final int c1, final Matrix X) {
        for (var r = r0; r <= r1; r++) {
            for (var c = c0; c <= c1; c++) {
                this.A[r][c] = X.A[r - r0][c - c0];
            }
        }
    }

    /**
     * Set a submatrix.
     *
     * @param r0
     *            Initial row index
     * @param r1
     *            Final row index
     * @param c
     *            Array of column indices.
     * @param X
     *            A(r0:r1,c(:))
     */
    public void setMatrix(final int r0, final int r1, final int[] c, final Matrix X) {
        for (var i = r0; i <= r1; i++) {
            for (var j = 0; j < c.length; j++) {
                this.A[i][c[j]] = X.A[i - r0][j];
            }
        }
    }

    /**
     * Set a submatrix.
     *
     * @param r
     *            Array of row indices.
     * @param c0
     *            Initial column index
     * @param c1
     *            Final column index
     * @param X
     *            A(r(:),c0:c1)
     */
    public void setMatrix(final int[] r, final int c0, final int c1, final Matrix X) {
        for (var i = 0; i < r.length; i++) {
            for (var j = c0; j <= c1; j++) {
                this.A[r[i]][j] = X.A[i][j - c0];
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
    public void setMatrix(final int[] r, final int[] c, final Matrix X) {
        for (var i = 0; i < r.length; i++) {
            for (var j = 0; j < c.length; j++) {
                this.A[r[i]][c[j]] = X.get(i, j);
            }
        }
    }

    /**
     * Solve A*X = B.
     *
     * @param B
     *            right hand side
     * @return solution if A is square, least squares solution otherwise.
     *         Returns null if no such solution exists (matrix is singular or
     *         rank deficient).
     */
    public Matrix solve(final Matrix B) {
        return (nRows == nColumns ? lu().solve(B) : qr().solve(B));
    }


    /**
     * Solve X*A = B, which is also A'*X' = B'.
     *
     * @param B
     *            right hand side
     * @return solution if A is square, least squares solution otherwise.
     *         Returns null if no such solution exists (matrix is singular or
     *         rank deficient).
     */
    public Matrix solveTranspose(final Matrix B) {
        return this.transpose().solve(B.transpose());
    }

    /**
     * Singular Value Decomposition.
     *
     * @return SingularValueDecomposition
     * @see SingularValueDecomposition
     */
    public SingularValueDecomposition svd() {
        return new SingularValueDecomposition(this);
    }

    /**
     * Multiply a matrix by a scalar, C = s*A.
     *
     * @param s
     *            scalar
     * @return new Matrix s*A
     */
    public Matrix times(final double s) {
        final var M = new Matrix(nRows, nColumns);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                M.A[r][c] = s * this.A[r][c];
            }
        }
        return M;
    }

    /**
     * Linear algebraic matrix multiplication, A * B.
     *
     * @param B
     *            another matrix
     * @return Matrix product, A * B
     * @throws
     *   IndexOutOfBoundsException Iff matrix dimension don't agree
     */
    public Matrix times(final Matrix B) {
        if (B.nRows != nColumns) {
            throw new IndexOutOfBoundsException(B.nRows);
        } else {
            final var X = new Matrix(nRows, B.nColumns);
            for (var c1 = 0; c1 < B.nColumns; c1++) {
                for (var r = 0; r < nRows; r++) {
                    var s = 0D;
                    for (var c2 = 0; c2 < nColumns; c2++) {
                        s += this.A[r][c2] * B.A[c2][c1];
                    }
                    X.A[r][c1] = s;
                }
            }
            return X;
        }
    }

    /**
     * Multiply a matrix by a scalar in place, A = s*A.
     *
     * @param s
     *            scalar
     */
    public void timesEquals(final double s) {
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                this.A[r][c] *= s;
            }
        }
    }

    /**
     * Returns a human-readable representation of this Matrix,
     * with each value in a row separated by tabulator.
     * @return String representation
     */
    @Override
    public String toString() {
        final var st = new StringBuffer();
        for (final var row : this.A) {
            for (final var d : row) {
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
    public double trace() {
        final var dim = nRows < nColumns ? nRows : nColumns;
        var trace = 0D;
        for (var i = 0; i < dim; i++) {
            trace += this.A[i][i];
        }
        return trace;
    }

    /**
     * Applies the given operator to all elements, returning a new Matrix.
     *
     * @param operator
     *            Operator to be applied to this Matrix and B
     * @return new Matrix with the result
     * @throws NullPointerException
     *             iff operator == null or B == null
     * @see #transformEquals(DoubleUnaryOperator)
     */
    public Matrix transform(final DoubleUnaryOperator operator) {
        final var M = new Matrix(nRows, nColumns);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                M.A[r][c] = operator.applyAsDouble(this.A[r][c]);
            }
        }
        return M;
    }

    /**
     * Applies the given operator to all elements, modifying this Matrix.
     *
     * @param operator
     *            Operator to be applied
     * @throws NullPointerException
     *             iff operator == null
     * @see #transform(DoubleUnaryOperator)
     */
    public void transformEquals(final DoubleUnaryOperator operator) {
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                this.A[r][c] = operator.applyAsDouble(this.A[r][c]);
            }
        }
    }

    /**
     * Applies the given operator to all elements, returning a new Matrix.
     *
     * @param B
     *            another Matrix
     * @param operator
     *            Operator to be applied
     * @return new Matrix with the result
     * @throws NullPointerException
     *             iff operator == null
     * @see #transformEquals(Matrix, DoubleBinaryOperator)
     */
    public Matrix transform(final FunctionalMatrix B, final DoubleBinaryOperator operator) {
        final var M = new Matrix(nRows, nColumns);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                M.A[r][c] = operator.applyAsDouble(this.get(r, c), B.get(r, c));
            }
        }
        return M;
    }

    /**
     * Applies the given operator to all elements, modifying this Matrix.
     *
     * @param B
     *            another Matrix
     * @param operator
     *            Operator to be applied to this Matrix and B
     * @throws NullPointerException
     *             iff operator == null or B == null
     * @see #transform(Matrix, DoubleBinaryOperator)
     */
    public void transformEquals(final FunctionalMatrix B, final DoubleBinaryOperator operator) {
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                this.A[r][c] = operator.applyAsDouble(this.get(r, c), B.get(r, c));
            }
        }
    }

    /**
     * Matrix transpose.
     *
     * @return new Matrix A'
     */
    @Override
    public Matrix transpose() {
        final var M = new Matrix(nColumns, nRows);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                M.A[c][r] = this.A[r][c];
            }
        }
        return M;
    }

    /**
     * Transposes this matrix.
     */
    public void transposeThis() {
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < r; c++) {
                final var t = this.A[r][c];
                this.A[r][c] = this.A[c][r];
                this.A[c][r] = t;
            }
        }
    }

    /**
     * Unary minus.
     *
     * @return new Matrix -A
     */
    @Override
    public Matrix uminus() {
        final var M = new Matrix(nRows, nColumns);
        for (var r = 0; r < nRows; r++) {
            for (var c = 0; c < nColumns; c++) {
                M.A[r][c] = -this.A[r][c];
            }
        }
        return M;
    }
}
