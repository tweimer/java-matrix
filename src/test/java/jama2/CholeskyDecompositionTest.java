package jama2;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class CholeskyDecompositionTest {
    @Test
    void testCholesky() {
        final double[][] pvals = { { 4.0, 1.0, 1.0 },
                                   { 1.0, 2.0, 3.0 },
                                   { 1.0, 3.0, 6.0 } };
        final var A = new Matrix(pvals);
        // Cholesky Decomposition
        final var chol = A.chol();
        final var L = chol.getL();
        // Check A==L*L^T
        assertEquals(A, L.times(L.transpose()), "CholeskyDecomposition: incorrect Cholesky decomposition calculation");

        // Solve A * X = I
        final var I = Matrix.identity(3);
        final var X = chol.solve(I);
        I.minusEquals(A.times(X));
        final var normInf = I.normInf();
        assertEquals(0.0, normInf, 1E-14, "CholeskyDecomposition solve(): incorrect Choleskydecomposition solve calculation"); //$NON-NLS-1$

    }

}
