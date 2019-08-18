package jama2;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class CholeskyDecompositionTest {
    @Test
    void testCholesky() {
        final double[][] pvals = { { 4D, 1D, 1D }, { 1D, 2D, 3D }, { 1D, 3D, 6D } };
        // Cholesky Decomposition
        final var A = new Matrix(pvals);
        final var chol = A.chol();
        final var L = chol.getL();
        // Check A==L*L^T
        assertEquals(A, L.times(L.transpose()), "CholeskyDecomposition: incorrect Cholesky decomposition calculation");

        // Solve A * X = I
        final var I = Matrix.identity(3);
        final var X = chol.solve(I);
        assertEquals(I, A.times(X), "CholeskyDecomposition solve(): incorrect Choleskydecomposition solve calculation"); //$NON-NLS-1$

    }

}
