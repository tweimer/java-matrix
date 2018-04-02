package jama2;

/**
 * @author tobias
 *
 */
@FunctionalInterface
public interface ISolver {
    /**
     * Solves an equation system <code>A*X=B</code> for a fixed Matrix A and a given Matrix B.
     * @param B Given matrix
     * @return Solution of the equation system <code>A*X=B</code>.
     */
    Matrix solve(Matrix B);
}
