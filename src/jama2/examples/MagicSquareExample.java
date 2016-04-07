package jama2.examples;

import jama2.EigenvalueDecomposition;
import jama2.LUDecomposition;
import jama2.Matrix;
import jama2.QRDecomposition;
import static jama2.util.Maths.eps;

import java.text.DecimalFormat;

/**
 * Example of use of Matrix Class, featuring magic squares.
 * 
 * @author The MathWorks, Inc. and the National Institute of Standards and
 *         Technology.
 * @version 8. Februar 2016
 * @see http://math.nist.gov/javanumerics/jama/
 **/
public class MagicSquareExample
{
    /**
     * Format double with Fw.d.
     * 
     * @param x
     *            integer
     * @param w
     *            width
     * @param d
     *            fraction digits
     * @return formatted String
     * **/
    public static String fixedWidthDoubletoString(final double x, final int w, final int d)
    {
        final DecimalFormat fmt = new DecimalFormat();
        fmt.setMaximumFractionDigits(d);
        fmt.setMinimumFractionDigits(d);
        fmt.setGroupingUsed(false);
        String s = fmt.format(x);
        while (s.length() < w)
        {
            s = " " + s; //$NON-NLS-1$
        }
        return s;
    }

    /**
     * Format integer with Iw.
     * 
     * @param n
     *            integer
     * @param w
     *            width
     * @return formatted String
     **/
    public static String fixedWidthIntegertoString(final int n, final int w)
    {
        String s = Integer.toString(n);
        while (s.length() < w)
        {
            s = " " + s; //$NON-NLS-1$
        }
        return s;
    }

    /**
     * Generate magic square test matrix.
     * 
     * @param n
     * @return Matrix
     **/
    public static Matrix magic(final int n)
    {

        final double[][] M = new double[n][n];

        // Odd order
        if ((n % 2) != 0)
        {
            final int b = (n + 1), a = b / 2;
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    M[i][j] = (n * ((i + j + a) % n)) + ((i + (2 * j) + b) % n) + 1;
                }
            }
        }
        // Doubly Even Order
        else if ((n % 4) == 0)
        {
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    M[i][j] = ((((i + 1) / 2) % 2) == (((j + 1) / 2) % 2)) ? ((n * n) - (n * i) - j) : ((n * i) + j + 1);
                }
            }
        }
        // Sing|ly Even Order
        else
        {
            final int p = n / 2, k = (n - 2) / 4;
            final Matrix A = MagicSquareExample.magic(p);
            for (int j = 0; j < p; j++)
            {
                for (int i = 0; i < p; i++)
                {
                    final double aij = M[i][j] = A.get(i, j);
                    M[i][j + p] = aij + (2 * p * p);
                    M[i + p][j] = aij + (3 * p * p);
                    M[i + p][j + p] = aij + (p * p);
                }
            }
            for (int i = 0; i < p; i++)
            {
                for (int j = 0; j < k; j++)
                {
                    final double t = M[i][j];
                    M[i][j] = M[i + p][j];
                    M[i + p][j] = t;
                }
                for (int j = (n - k) + 1; j < n; j++)
                {
                    final double t = M[i][j];
                    M[i][j] = M[i + p][j];
                    M[i + p][j] = t;
                }
            }
            double t = M[k][0];
            M[k][0] = M[k + p][0];
            M[k + p][0] = t;

            t = M[k][k];
            M[k][k] = M[k + p][k];
            M[k + p][k] = t;
        }
        return new Matrix(M);
    }

    /**
     * Tests LU, QR, SVD and symmetric Eig decompositions.
     * 
     * n = order of magic square.<br />
     * trace = diagonal sum, should be the magic sum, (n^3 + n)/2.<br />
     * max_eig = maximum eigenvalue of (A + A')/2, should equal trace.<br />
     * rank = linear algebraic rank, should equal n if n is odd, be less than n
     * if n is even.<br />
     * cond = L_2 condition number, ratio of singular values.<br />
     * lu_res = test of LU factorization, norm1(L*U-A(p,:))/(n*eps).<br />
     * qr_res = test of QR factorization, norm1(Q*R-A)/(n*eps).
     * 
     * @param argv
     *            unused
     */
    public static void main(final String argv[])
    {
        System.out.println("\n    Test of Matrix Class, using magic squares."); //$NON-NLS-1$
        System.out.println("    See MagicSquareExample.main() for an explanation."); //$NON-NLS-1$
        System.out.println("\n      n     trace       max_eig   rank        cond      lu_res      qr_res\n"); //$NON-NLS-1$

        final long start_time = System.currentTimeMillis();
        for (int n = 3; n <= 32; n++)
        {
            System.out.print(MagicSquareExample.fixedWidthIntegertoString(n, 7));

            final Matrix M = MagicSquareExample.magic(n);

            final int t = (int) M.trace();
            System.out.print(MagicSquareExample.fixedWidthIntegertoString(t, 10));

            final EigenvalueDecomposition E = M.plus(M.transpose()).times(0.5).eig();
            final double[] d = E.getRealEigenvalues();
            System.out.print(MagicSquareExample.fixedWidthDoubletoString(d[n - 1], 14, 3));

            final int r = M.rank();
            System.out.print(MagicSquareExample.fixedWidthIntegertoString(r, 7));

            final double c = M.cond();
            System.out.print(c < (1 / eps) ? MagicSquareExample.fixedWidthDoubletoString(c, 12, 3) : "         Inf"); //$NON-NLS-1$

            final LUDecomposition LU = M.lu();
            final Matrix L = LU.getL();
            final Matrix U = LU.getU();
            final int[] p = LU.getPivot();
            Matrix R = L.times(U).minus(M.getMatrix(p, 0, n - 1));
            double res = R.norm1() / (n * eps);
            System.out.print(MagicSquareExample.fixedWidthDoubletoString(res, 12, 3));

            final QRDecomposition QR = M.qr();
            final Matrix Q = QR.getQ();
            R = QR.getR();
            R = Q.times(R).minus(M);
            res = R.norm1() / (n * eps);
            System.out.println(MagicSquareExample.fixedWidthDoubletoString(res, 12, 3));
        }
        final long stop_time = System.currentTimeMillis();
        final double etime = (stop_time - start_time) / 1000.;
        System.out.println("\nElapsed Time = " + MagicSquareExample.fixedWidthDoubletoString(etime, 12, 3) + " seconds"); //$NON-NLS-1$ //$NON-NLS-2$
        System.out.println("Adios"); //$NON-NLS-1$
    }

}
