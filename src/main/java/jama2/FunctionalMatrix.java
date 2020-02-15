package jama2;

/**
 * This functional interface represent a {@link Matrix}.
 * This is the base interface for all implementations.
 * 
 * 
 * @author The MathWorks, Inc. and the National Institute of Standards and
 *         Technology.
 * @version 2.0
 * @see jama2.Matrix#Matrix(int, int, IMatrix)
 * @see <a href="http://tweimer.github.io/java-matrix/">java-matrix</a>
 */
@FunctionalInterface
public interface FunctionalMatrix {
    /**
     * Get a single element.
     * @param r row index
     * @param c column index
     * @return value
     * @see Matrix#Matrix(int, int, FunctionalMatrix)
     */
    double get(int r, int c);
    
    /**
     * Matrix transpose
     * @return
     *    transpose
     */
    default FunctionalMatrix transpose() {
        return (r, c) -> get(c, r);
    }
    
    /**
     * C = A + B
     * @param m
     *   Another functional Matrix
     * @return
     *   {@code get(r, c) + m.get(r, c)}
     */
    default FunctionalMatrix plus(FunctionalMatrix m) {
        return (r, c) -> this.get(r, c) + m.get(r, c);
    }

    /**
     * C = A - B
     * @param m
     *   Another functional Matrix
     * @return
     *   {@code get(r, c) - m.get(r, c)}
     */
    default FunctionalMatrix minus(FunctionalMatrix m) {
        return (r, c) -> this.get(r, c) - m.get(r, c);
    }

    /**
     * Multiply a matrix by a scalar, C = s*A.
     *
     * @param s scalar
     * @return Returns a new FunctionalMatrix s*A.
     */
    default FunctionalMatrix times(double s) {
        return (r, c) -> s * get(r, c);
    }
    
    /**
     * Unary minus.
     *
     * @return new FunctionalMatrix -A
     */
    default FunctionalMatrix uminus() {
        return (r, c) -> -get(r, c);
    }
}
