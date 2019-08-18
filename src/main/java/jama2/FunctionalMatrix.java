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
     * @param i row index
     * @param j column index
     * @return value
     * @see Matrix#Matrix(int, int, FunctionalMatrix)
     */
    double get(int i, int j);
}
