package jama2;

import java.io.Serializable;
import java.util.Arrays;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;

/**
 * This class represents a diagonal matrix of a fixed square size.
 *
 * @author Tobias Weimer
 *
 */
public class DiagonalMatrix implements Serializable, Cloneable, IMatrix {
    private static final long serialVersionUID = 1L;

    /** array of the diagonal elements */
    final double[] diag;

    /**
     * Constructs a diagonal matrix
     *
     * @param size Size of the matrix
     */
    public DiagonalMatrix(final int size) {
        diag = new double[size];
    }

    /**
     * Makes a deep copy of the underlying array
     *
     * @param d Another diagonal Matrix
     */
    public DiagonalMatrix(final DiagonalMatrix D) {
        diag = Arrays.copyOf(D.diag, D.diag.length);
    }
    
    /**
     * Creates a diagonal matrix without copying the diagonal elements
     * @param diag array of diagonal elements
     */
    DiagonalMatrix(final double[] diag) {
        this.diag = diag;
    }

    /**
     * Creates a diagonal Matrix.
     *
     * @param diag array diagonal elements
     * @return diagonal Matrix
     * @throws NullPointerException Iff {@code diag == null}
     */
    public static DiagonalMatrix diag(final double... diag) {
        final var diag2 = Arrays.copyOf(diag, diag.length);
        return new DiagonalMatrix(diag2);
    }

    /**
     * Get the element at position (i,j).
     *
     * @param i row index
     * @param j column index
     * @return Iff {@code i == j}, the corresponding diagonal element is returned,
     *         otherwise 0
     * @throws ArrayIndexOutOfBoundsException Iff any index is invalid
     * @see #get(int)
     */
    @Override
    public double get(final int i, final int j) {
        checkRange(i);
        checkRange(j);
        return i == j ? diag[i] : 0D;
    }

    /**
     * Get the value of the i-th diagonal entry
     *
     * @param i index
     * @return i-th diagonal value
     * @see #get(int, int)
     * @throws ArrayIndexOutOfBoundsException Iff the index is invalid
     */
    public double get(final int i) {
        checkRange(i);
        return diag[i];
    }

    /**
     * Get the dimension of the Matrix
     *
     * @return
     */
    public int getSize() {
        return diag.length;
    }

    /**
     * Set the i-th diagonal value
     *
     * @param i diagonal index
     * @param d value to be set
     * @throws IndexOutOfBoundsException Iff i is an invalid index
     */
    public void set(final int i, final double d) {
        diag[i] = d;
    }

    /*
     * (non-Javadoc)
     *
     * @see java.lang.Object#hashCode()
     */
    @Override
    public int hashCode() {
        return Arrays.hashCode(diag);
    }

    /*
     * (non-Javadoc)
     *
     * @see java.lang.Object#equals(java.lang.Object)
     */
    @Override
    public boolean equals(final Object obj) {
        if (this == obj)
            return true;
        if (!(obj instanceof DiagonalMatrix))
            return false;
        final DiagonalMatrix other = (DiagonalMatrix) obj;
        return Arrays.equals(diag, other.diag);
    }

    /**
     * Multiply a matrix by a scalar, C = s*A.
     *
     * @param s scalar
     * @return new Matrix s*A
     */
    public DiagonalMatrix times(final double s) {
        final var D2 = new DiagonalMatrix(diag.length);
        // Just multiply the diagonal values
        for (var i = 0; i < diag.length; i++) {
            D2.set(i, s * get(i));
        }
        return D2;
    }

    /**
     * Multiply a matrix by a scalar in place, A = s*A.
     *
     * @param s scalar
     */
    public void timesEquals(final double s) {
        for (var i = 0; i < diag.length; i++) {
            diag[i] *= s;
        }
    }

    /**
     * C = A + B.
     *
     * @param D another matrix
     * @return new DiagonalMatrix A + B
     * @see #plusEquals
     */
    public DiagonalMatrix plus(final DiagonalMatrix D) {
        checkMatrixDimensions(D);
        final var D2 = new DiagonalMatrix(diag.length);
        for (var i = 0; i < diag.length; i++) {
            D2.diag[i] = this.diag[i] + D.diag[i];
        }
        return D2;
    }

    /**
     * A = A + B.
     *
     * @param D another diagonal matrix
     * @see #plus
     */
    public void plusEquals(final DiagonalMatrix D) {
        checkMatrixDimensions(D);
        for (var i = 0; i < diag.length; i++) {
            diag[i] += D.diag[i];
        }
    }

    /**
     * C = A - B.
     *
     * @param B another matrix
     * @return new Matrix A - B
     * @see #minusEquals
     */
    public DiagonalMatrix minus(final DiagonalMatrix D) {
        checkMatrixDimensions(D);
        final var D2 = new DiagonalMatrix(diag.length);
        for (var i = 0; i < diag.length; i++) {
            D2.diag[i] = this.diag[i] - D.diag[i];
        }
        return D2;
    }

    /**
     * A = A - B.
     *
     * @param B another matrix
     * @see #minus
     */
    public void minusEquals(final DiagonalMatrix D) {
        checkMatrixDimensions(D);
        for (var i = 0; i < diag.length; i++) {
            diag[i] -= D.diag[i];
        }
    }

    /**
     * Applies the given operator to all elements, returning a new Matrix.
     *
     * @param operator Operator to be applied to this Matrix and B
     * @return new Matrix with the result
     * @throws NullPointerException iff operator == null or B == null
     * @see #transformEquals(DoubleUnaryOperator)
     */
    public DiagonalMatrix transform(final DoubleUnaryOperator operator) {
        final var D = new DiagonalMatrix(diag.length);
        for (var i = 0; i < diag.length; i++) {
            D.diag[i] = operator.applyAsDouble(diag[i]);
        }
        return D;
    }

    /**
     * Applies the given operator to all elements, modifying this Matrix.
     *
     * @param operator Operator to be applied
     * @throws NullPointerException iff operator == null
     * @see #transform(DoubleUnaryOperator)
     */
    public void transformEquals(final DoubleUnaryOperator operator) {
        for (var i = 0; i < diag.length; i++) {
            diag[i] = operator.applyAsDouble(diag[i]);
        }
    }

    /**
     * Applies the given operator to all elements, returning a new Matrix.
     *
     * @param D  another diagonal Matrix
     * @param operator Operator to be applied
     * @return new Matrix with the result
     * @throws NullPointerException iff operator == null
     * @see #transformEquals(Matrix, DoubleBinaryOperator)
     */
    public DiagonalMatrix transform(final DiagonalMatrix D, final DoubleBinaryOperator operator) {
        checkMatrixDimensions(D);
        final var D2 = new DiagonalMatrix(diag.length);
        for (var i = 0; i < diag.length; i++) {
            D2.diag[i] = operator.applyAsDouble(diag[i], D.diag[i]);
        }
        return D2;
    }

    /**
     * Applies the given operator to all elements, modifying this Matrix.
     *
     * @param D  another diagonal Matrix
     * @param operator Operator to be applied to this Matrix and B
     * @throws NullPointerException iff operator == null or B == null
     * @see #transform(Matrix, DoubleBinaryOperator)
     */
    public void transformEquals(final DiagonalMatrix D, final DoubleBinaryOperator operator) {
        checkMatrixDimensions(D);
        for (var i = 0; i < diag.length; i++) {
            diag[i] = operator.applyAsDouble(diag[i], D.diag[i]);
        }
    }

    /**
     * Linear algebraic matrix multiplication, A * B.
     *
     *
     * @param D another diagonal matrix
     * @return Matrix product, A * B, null if matrix dimensions don't agree. Note
     *         that the product of two diagonal matrixes is a diagonal matrix.
     */
    public DiagonalMatrix times(final DiagonalMatrix D) {
        checkMatrixDimensions(D);
        final var D2 = new DiagonalMatrix(diag.length);
        // We just have to multiply the diagonal entries
        for (var i = 0; i < diag.length; i++) {
            D2.diag[i] = diag[i] * D.diag[i];
        }
        return D2;
    }

    /**
     * The inverse of the diagonal matrix.
     *
     * @return the inverse, which is itsself a diagonal Matrix
     * @throws ArithmeticException Iff any diagonal entry is zero, hence there is no
     *                             inverse
     */
    public DiagonalMatrix inverse() {
        final var D2 = new DiagonalMatrix(diag.length);
        // We just have to inverse the diagonal entries
        for (var i = 0; i < diag.length; i++) {
            // Throws ArithmeticException iff diag[i] == 0D
            D2.diag[i] = 1d / diag[i];
        }
        return D2;
    }

    /**
     * Matrix trace.
     *
     * @return sum of the diagonal elements.
     */
    public double trace() {
        var trace = 0D;
        for (final var d : diag) {
            trace += d;
        }
        return trace;
    }

    /**
     * Returns the determinant. The determinant can be computed by multiplying all
     * diagonal items.
     *
     * @return determinant
     */
    public double det() {
        var det = 1D;
        for (final var d : diag) {
            det *= d;
        }
        return det;
    }
    

    /**
     * Check if index is in range of the Matrix
     *
     * @param i index
     * @return true iff in range
     */
    private boolean isInRange(final int i) {
        return i >= 0 && i < diag.length;
    }

    /**
     * Check range
     *
     * @param i index
     * @throws ArrayIndexOutOfBoundsException Iff the index is not valid
     */
    private void checkRange(final int i) {
        if (!isInRange(i))
            throw new ArrayIndexOutOfBoundsException(i);
    }
    
    /**
     * Check if size(A) == size(D).
    *
    * @param D another Matrix
    * @throws IllegalArgumentException
    *             if they don't agree.
    */
   private void checkMatrixDimensions(final DiagonalMatrix D) {
       if (this.diag.length != D.diag.length) {
           throw new IllegalArgumentException("Matrix dimensions must agree."); //$NON-NLS-1$
       }
   }

}
