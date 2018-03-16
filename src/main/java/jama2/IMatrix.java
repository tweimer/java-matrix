/**
 * 
 */
package jama2;

/**
 * This functional interface represent a {@link Matrix} and can be used to build a Matrix.
 * @author The MathWorks, Inc. and the National Institute of Standards and
 *         Technology.
 * @version 2.0
 * @see Matrix#Matrix(int, int, IMatrix)
 * @see <a href="http://tweimer.github.io/java-matrix/">java-matrix</a>
 */
@FunctionalInterface
public interface IMatrix {
    /**
     * Get a single element.
     * @param i row index
     * @param j column index
     * @return value
     * @see Matrix#Matrix(int, int, IMatrix)
     */
    double get(int i, int j);
}
