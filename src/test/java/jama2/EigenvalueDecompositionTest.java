package jama2;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import jama2.test.TestMatrix;

class EigenvalueDecompositionTest {

    @Test
    void test() {
        /*
         * EigenvalueDecomposition Eig = A.eig(); Matrix D = Eig.getD(), V = Eig.getV();
         * try { // Check A*V==V*D TestMatrix.check(A.times(V), V.times(D));
         * TestMatrix.try_success("EigenvalueDecomposition (symmetric)...");
         * //$NON-NLS-1$ } catch (final RuntimeException e) {
         * TestMatrix.try_failure("EigenvalueDecomposition (symmetric)...",
         * //$NON-NLS-1$ "incorrect symmetric Eigenvalue decomposition calculation");
         * //$NON-NLS-1$ }
         *
         * try { // Eigenvalues A = new Matrix(evals); Eig = A.eig(); D = Eig.getD(); V
         * = Eig.getV(); TestMatrix.check(A.times(V), V.times(D));
         * TestMatrix.try_success("EigenvalueDecomposition (nonsymmetric)...");
         * //$NON-NLS-1$ } catch (final RuntimeException e) {
         * TestMatrix.try_failure("EigenvalueDecomposition (nonsymmetric)...",
         * //$NON-NLS-1$ "incorrect nonsymmetric Eigenvalue decomposition calculation");
         * //$NON-NLS-1$ }
         *
         */
    }

  /*  @Test(timeout = 1000)
    @DisplayName("Test bad Eigenvalues")
    void testBadEigs() {
        final double[][] badeigs = { { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 1 }, { 0, 0, 0, 1, 0 }, { 1, 1, 0, 0, 1 },
                { 1, 0, 1, 0, 1 } };
        try {
            final Matrix bA = new Matrix(badeigs);
            System.out.println();
            System.out.println("Testing Eigenvalue; If this hangs, we've failed"); //$NON-NLS-1$
            bA.eig();
            TestMatrix.try_success("EigenvalueDecomposition (hang)..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("EigenvalueDecomposition (hang)...", //$NON-NLS-1$
                    "incorrect termination"); //$NON-NLS-1$
        }
    }*/

}
