package jama2.test;

import jama2.CholeskyDecomposition;
import jama2.EigenvalueDecomposition;
import jama2.LUDecomposition;
import jama2.Matrix;
import jama2.QRDecomposition;
import jama2.SingularValueDecomposition;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

import static jama2.util.Maths.eps;

/**
 * TestMatrix tests the functionality of the Jama Matrix class and associated
 * decompositions.
 * <P>
 * Run the test from the command line using <BLOCKQUOTE>
 * 
 * <PRE>
 * <CODE>
 *  java jama2.test.TestMatrix 
 * </CODE>
 * </PRE>
 * 
 * </BLOCKQUOTE>
 * 
 * Detailed output is provided indicating the functionality being tested and
 * whether the functionality is correctly implemented. Exception handling is
 * also tested.
 * <P>
 * The test is designed to run to completion and give a summary of any
 * implementation errors encountered. The final output should be: <BLOCKQUOTE>
 * 
 * <PRE>
 * <CODE>
 *       TestMatrix completed.
 *       Total errors reported: n1
 *       Total warning reported: n2
 * </CODE>
 * </PRE>
 * 
 * </BLOCKQUOTE>
 * 
 * If the test does not run to completion, this indicates that there is a
 * substantial problem within the implementation that was not anticipated in the
 * test design. The stopping point should give an indication of where the
 * problem exists.
 * 
 * @author The MathWorks, Inc. and the National Institute of Standards and
 *         Technology.
 * @version 2.0
 * @see <a href="http://tweimer.github.io/java-matrix/">java-matrix</a>
 **/
public class TestMatrix {
    static int errorCount = 0, warningCount = 0;
    /**
     * 
     */
    final static String tmpname = "TMPMATRIX.serial"; //$NON-NLS-1$

    /**
     * Check magnitude of difference of scalars.
     * 
     * @param x
     * @param y
     * @throws RuntimeException
     *             on failure
     */
    private static void check(final double x, final double y) {
        if ((x == 0) && (Math.abs(y) < (10 * eps))) {
            return;
        }
        if ((y == 0) && (Math.abs(x) < (10 * eps))) {
            return;
        }
        if (Math.abs(x - y) > (10 * eps * Math.max(Math.abs(x), Math.abs(y)))) {
            throw new RuntimeException("The difference x-y is too large: x = " + x + "  y = " + y); //$NON-NLS-1$ //$NON-NLS-2$
        }
    }

    /**
     * Check norm of difference of "vectors".
     * 
     * @param x
     * @param y
     * @throws RuntimeException
     *             on failure
     */
    private static void check(final double[] x, final double[] y) {
        if (x.length == y.length) {
            for (int i = 0; i < x.length; i++) {
                TestMatrix.check(x[i], y[i]);
            }
        } else {
            throw new RuntimeException("Attempt to compare vectors of different lengths"); //$NON-NLS-1$
        }
    }

    /**
     * Check norm of difference of arrays.
     * 
     * @param x
     *            Matrix
     * @param y
     *            Matrix
     * @throws RuntimeException
     *             on failure
     */
    private static void check(final double[][] x, final double[][] y) {
        final Matrix A = new Matrix(x), B = new Matrix(y);
        TestMatrix.check(A, B);
    }

    /**
     * Check norm of difference of Matrices.
     * 
     * @param X
     * @param Y
     * @throws RuntimeException
     *             on failure
     */
    private static void check(final Matrix X, final Matrix Y) {
        final double Xnorm1 = X.norm1(), Ynorm1 = Y.norm1();
        if ((Xnorm1 == 0D) && (Ynorm1 < (10D * eps))) {
            return;
        }
        if ((Ynorm1 == 0D) && (Xnorm1 < (10D * eps))) {
            return;
        }
        if (X.minus(Y).norm1() > (1000 * eps * Math.max(Xnorm1, Ynorm1))) {
            throw new RuntimeException("The norm of (X-Y) is too large: " + X.minus(Y).norm1()); //$NON-NLS-1$
        }
    }

    /**
     * Print appropriate messages for unsuccessful outcome try
     * 
     * @param s
     *            failure string
     * @param e
     *            error message
     */
    private static void try_failure(final String s, final String e) {
        System.out.println(">\t" + s + "*** failure ***\n>\t\tMessage: " + e); //$NON-NLS-1$ //$NON-NLS-2$
        errorCount++;
    }

    /**
     * Print appropriate messages for successful outcome try
     * 
     * @param s
     *            failure string
     */
    private static void try_success(final String s) {
        System.out.println(">\t" + s + "success"); //$NON-NLS-1$ //$NON-NLS-2$
    }

    /**
     * Print appropriate messages for successful outcome try
     * 
     * @param s
     *            failure string
     * @param e
     *            error message
     */
    private static void try_success(final String s, final String e) {
        TestMatrix.try_success(s);
        System.out.println(">\t\tMessage: " + e); //$NON-NLS-1$
    }

    /**
     * Print appropriate messages for unsuccessful outcome try
     * 
     * @param s
     *            failure string
     * @param e
     *            failure message
     * 
     * @return number of errors
     */
    private static void try_warning(final String s, final String e) {
        System.out.println(">\t" + s + "*** warning ***\n>\tMessage: " + e); //$NON-NLS-1$ //$NON-NLS-2$
        warningCount++;
    }

    /**
     * @param argv
     *            unused
     */
    public static void main(final String argv[]) {
        Matrix A, blablabla, M;
        // Uncomment this to test IO in a different locale.
        // Locale.setDefault(Locale.GERMAN);
        double tmp;
        final double[] columnwise = { 1D, 2D, 3D, 4D, 5D, 6D, 7D, 8D, 9D, 10D, 11D, 12D };
        final double[] rowwise = { 1., 4., 7., 10., 2., 5., 8., 11., 3., 6., 9., 12. };
        final double[][] avals = { { 1., 4., 7., 10. }, { 2., 5., 8., 11. }, { 3., 6., 9., 12. } };
        final double[][] rankdef = avals;
        final double[][] tvals = { { 1D, 2D, 3D }, { 4D, 5D, 6D }, { 7D, 8D, 9D }, { 10D, 11D, 12D } };
        final double[][] subavals = { { 5D, 8D, 11D }, { 6D, 9D, 12D } };
        final double[][] rvals = { { 1D, 4D, 7D }, { 2D, 5D, 8D, 11D }, { 3D, 6D, 9D, 12D } };
        final double[][] pvals = { { 4D, 1D, 1D }, { 1D, 2D, 3D }, { 1D, 3D, 6D } };
        final double[][] ivals = { { 1D, 0D, 0D, 0D }, { 0D, 1D, 0D, 0D }, { 0D, 0D, 1D, 0D } };
        final double[][] evals = { { 0D, 1D, 0D, 0D }, { 1D, 0D, 2.e-7D, 0D }, { 0D, -2.e-7D, 0D, 1D },
                { 0D, 0D, 1D, 0D } };
        final double[][] square = { { 166D, 188D, 210D }, { 188D, 214D, 240D }, { 210D, 240D, 270D } };
        final double[][] sqSolution = { { 13D }, { 15D } };
        final double[][] condmat = { { 1D, 3D }, { 7D, 9D } };
        final double[][] badeigs = { { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 1 }, { 0, 0, 0, 1, 0 }, { 1, 1, 0, 0, 1 },
                { 1, 0, 1, 0, 1 } };
        final int rows = 3, cols = 4;

        // should trigger bad shape for construction with val
        final int invalidld = 5;

        // (raggedr,raggedc) should be out of bounds in ragged array
        final int raggedr = 0, raggedc = 4;

        // leading dimension of intended test Matrices
        final int validld = 3;

        // leading dimension which is valid, but nonconforming
        final int nonconformld = 4;

        // index ranges for sub Matrix
        final int ib = 1, ie = 2, jb = 1, je = 3;
        final int[] rowindexset = { 1, 2 }, badrowindexset = { 1, 3 }, columnindexset = { 1, 2, 3 },
                badcolumnindexset = { 1, 2, 4 };
        final double columnsummax = 33D, rowsummax = 30D, sumofdiagonals = 15D, sumofsquares = 650D;

        // Constructors and constructor-like methods:
        // double[], int
        // double[][]
        // int, int
        // int, int, double
        // int, int, double[][]
        // constructWithCopy(double[][])
        // random(int,int)
        // identity(int)

        System.out.println();
        System.out.println("Testing constructors and constructor-like methods..."); //$NON-NLS-1$
        try {
            // check that exception is thrown in packed constructor with invalid
            // length
            A = new Matrix(columnwise, invalidld);
            TestMatrix.try_failure("Catch invalid length in packed constructor... ", //$NON-NLS-1$
                    "exception not thrown for invalid input"); //$NON-NLS-1$
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_success("Catch invalid length in packed constructor... ", //$NON-NLS-1$
                    e.getMessage());
        }

        try {
            // check that exception is thrown in default constructor
            // if input array is 'ragged'
            A = new Matrix(rvals);
            tmp = A.get(raggedr, raggedc);
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_success("Catch ragged input to default constructor... ", //$NON-NLS-1$
                    e.getMessage());
        } catch (final ArrayIndexOutOfBoundsException e) {
            TestMatrix.try_failure("Catch ragged input to constructor... ", //$NON-NLS-1$
                    "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later"); //$NON-NLS-1$
        }

        try {
            // check that exception is thrown in constructWithCopy
            // if input array is 'ragged'
            A = Matrix.constructWithCopy(rvals);
            tmp = A.get(raggedr, raggedc);
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_success("Catch ragged input to constructWithCopy... ", //$NON-NLS-1$
                    e.getMessage());
        } catch (final ArrayIndexOutOfBoundsException e) {
            TestMatrix.try_failure("Catch ragged input to constructWithCopy... ", //$NON-NLS-1$
                    "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later"); //$NON-NLS-1$
        }

        A = new Matrix(columnwise, validld);
        Matrix B = new Matrix(avals);
        tmp = B.get(0, 0);
        avals[0][0] = 0D;
        Matrix C = B.minus(A);
        avals[0][0] = tmp;
        B = Matrix.constructWithCopy(avals);
        tmp = B.get(0, 0);
        avals[0][0] = 0D;
        if (tmp != B.get(0, 0)) {
            // check that constructWithCopy behaves properly
            TestMatrix.try_failure("constructWithCopy... ", //$NON-NLS-1$
                    "copy not effected... data visible outside"); //$NON-NLS-1$
        } else {
            TestMatrix.try_success("constructWithCopy... "); //$NON-NLS-1$
        }

        avals[0][0] = columnwise[0];
        final Matrix I = new Matrix(ivals);
        try {
            TestMatrix.check(I, Matrix.identity(3, 4));
            TestMatrix.try_success("identity... "); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("identity... ", //$NON-NLS-1$
                    "identity Matrix not successfully created"); //$NON-NLS-1$
        }

        // Access Methods:
        // getColumnDimension()
        // getRowDimension()
        // getArray()
        // getArrayCopy()
        // getColumnPackedCopy()
        // getRowPackedCopy()
        // get(int,int)
        // getMatrix(int,int,int,int)
        // getMatrix(int,int,int[])
        // getMatrix(int[],int,int)
        // getMatrix(int[],int[])
        // set(int,int,double)
        // setMatrix(int,int,int,int,Matrix)
        // setMatrix(int,int,int[],Matrix)
        // setMatrix(int[],int,int,Matrix)
        // setMatrix(int[],int[],Matrix)

        System.out.println();
        System.out.println("Testing access methods..."); //$NON-NLS-1$

        // Various get methods:
        B = new Matrix(avals);
        if (B.getRowDimension() != rows) {
            TestMatrix.try_failure("getRowDimension... ", ""); //$NON-NLS-1$ //$NON-NLS-2$
        } else {
            TestMatrix.try_success("getRowDimension... "); //$NON-NLS-1$
        }

        if (B.getColumnDimension() != cols) {
            TestMatrix.try_failure("getColumnDimension... ", ""); //$NON-NLS-1$ //$NON-NLS-2$
        } else {
            TestMatrix.try_success("getColumnDimension... "); //$NON-NLS-1$
        }

        // Check if getArray returns the same array
        B = new Matrix(avals);
        double[][] barray = B.getArray();
        if (barray != avals) {
            TestMatrix.try_failure("getArray... ", ""); //$NON-NLS-1$ //$NON-NLS-2$
        } else {
            TestMatrix.try_success("getArray... "); //$NON-NLS-1$
        }

        // Check if getArrayCopy returns a copy of the attay
        barray = B.getArrayCopy();
        if (barray == avals) {
            TestMatrix.try_failure("getArrayCopy... ", //$NON-NLS-1$
                    "data not (deep) copied"); //$NON-NLS-1$
        }

        try {
            TestMatrix.check(barray, avals);
            TestMatrix.try_success("getArrayCopy... "); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("getArrayCopy... ", //$NON-NLS-1$
                    "data not successfully (deep) copied"); //$NON-NLS-1$
        }

        double bpacked[] = B.getColumnPackedCopy();
        try {
            TestMatrix.check(bpacked, columnwise);
            TestMatrix.try_success("getColumnPackedCopy... "); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("getColumnPackedCopy... ", //$NON-NLS-1$
                    "data not successfully (deep) copied by columns"); //$NON-NLS-1$
        }

        bpacked = B.getRowPackedCopy();
        try {
            TestMatrix.check(bpacked, rowwise);
            TestMatrix.try_success("getRowPackedCopy... "); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("getRowPackedCopy... ", //$NON-NLS-1$
                    "data not successfully (deep) copied by rows"); //$NON-NLS-1$
        }
        try {
            tmp = B.get(B.getRowDimension(), B.getColumnDimension() - 1);
            TestMatrix.try_failure("get(int,int)... ", //$NON-NLS-1$
                    "OutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        } catch (final ArrayIndexOutOfBoundsException e) {
            try {
                tmp = B.get(B.getRowDimension() - 1, B.getColumnDimension());
                TestMatrix.try_failure("get(int,int)... ", //$NON-NLS-1$
                        "OutOfBoundsException expected but not thrown"); //$NON-NLS-1$
            } catch (final ArrayIndexOutOfBoundsException e1) {
                TestMatrix.try_success("get(int,int)... OutofBoundsException... "); //$NON-NLS-1$
            }
        } catch (final IllegalArgumentException e1) {
            TestMatrix.try_failure("get(int,int)... ", //$NON-NLS-1$
                    "OutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        }

        try {
            if (B.get(B.getRowDimension() - 1,
                    B.getColumnDimension() - 1) != avals[B.getRowDimension() - 1][B.getColumnDimension() - 1]) {
                TestMatrix.try_failure("get(int,int)... ", //$NON-NLS-1$
                        "Matrix entry (i,j) not successfully retreived"); //$NON-NLS-1$
            } else {
                TestMatrix.try_success("get(int,int)... "); //$NON-NLS-1$
            }
        } catch (final ArrayIndexOutOfBoundsException e) {
            TestMatrix.try_failure("get(int,int)... ", //$NON-NLS-1$
                    "Unexpected ArrayIndexOutOfBoundsException"); //$NON-NLS-1$
        }

        final Matrix SUB = new Matrix(subavals);
        try {
            M = B.getMatrix(ib, ie + B.getRowDimension() + 1, jb, je);
            TestMatrix.try_failure("getMatrix(int,int,int,int)... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        } catch (final ArrayIndexOutOfBoundsException e) {
            try {
                M = B.getMatrix(ib, ie, jb, je + B.getColumnDimension() + 1);
                TestMatrix.try_failure("getMatrix(int,int,int,int)... ", //$NON-NLS-1$
                        "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
            } catch (final ArrayIndexOutOfBoundsException e1) {
                TestMatrix.try_success("getMatrix(int,int,int,int)... ArrayIndexOutOfBoundsException... "); //$NON-NLS-1$
            }
        } catch (final IllegalArgumentException e1) {
            TestMatrix.try_failure("getMatrix(int,int,int,int)... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        }

        try {
            M = B.getMatrix(ib, ie, jb, je);
            try {
                TestMatrix.check(SUB, M);
                TestMatrix.try_success("getMatrix(int,int,int,int)... "); //$NON-NLS-1$
            } catch (final RuntimeException e) {
                TestMatrix.try_failure("getMatrix(int,int,int,int)... ", //$NON-NLS-1$
                        "submatrix not successfully retreived"); //$NON-NLS-1$
            }
        } catch (final ArrayIndexOutOfBoundsException e) {
            TestMatrix.try_failure("getMatrix(int,int,int,int)... ", //$NON-NLS-1$
                    "Unexpected ArrayIndexOutOfBoundsException"); //$NON-NLS-1$
        }

        try {
            M = B.getMatrix(ib, ie, badcolumnindexset);
            TestMatrix.try_failure("getMatrix(int,int,int[])... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        } catch (final ArrayIndexOutOfBoundsException e) {
            try {
                M = B.getMatrix(ib, ie + B.getRowDimension() + 1, columnindexset);
                TestMatrix.try_failure("getMatrix(int,int,int[])... ", //$NON-NLS-1$
                        "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
            } catch (final ArrayIndexOutOfBoundsException e1) {
                TestMatrix.try_success("getMatrix(int,int,int[])... ArrayIndexOutOfBoundsException... "); //$NON-NLS-1$
            }
        } catch (final IllegalArgumentException e1) {
            TestMatrix.try_failure("getMatrix(int,int,int[])... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        }

        try {
            M = B.getMatrix(ib, ie, columnindexset);
            try {
                TestMatrix.check(SUB, M);
                TestMatrix.try_success("getMatrix(int,int,int[])... "); //$NON-NLS-1$
            } catch (final RuntimeException e) {
                TestMatrix.try_failure("getMatrix(int,int,int[])... ", //$NON-NLS-1$
                        "submatrix not successfully retreived"); //$NON-NLS-1$
            }
        } catch (final ArrayIndexOutOfBoundsException e) {
            TestMatrix.try_failure("getMatrix(int,int,int[])... ", //$NON-NLS-1$
                    "Unexpected ArrayIndexOutOfBoundsException"); //$NON-NLS-1$
        }

        try {
            M = B.getMatrix(badrowindexset, jb, je);
            TestMatrix.try_failure("getMatrix(int[],int,int)... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        } catch (final ArrayIndexOutOfBoundsException e) {
            try {
                M = B.getMatrix(rowindexset, jb, je + B.getColumnDimension() + 1);
                TestMatrix.try_failure("getMatrix(int[],int,int)... ", //$NON-NLS-1$
                        "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
            } catch (final ArrayIndexOutOfBoundsException e1) {
                TestMatrix.try_success("getMatrix(int[],int,int)... ArrayIndexOutOfBoundsException... "); //$NON-NLS-1$
            }
        } catch (final IllegalArgumentException e1) {
            TestMatrix.try_failure("getMatrix(int[],int,int)... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        }

        try {
            M = B.getMatrix(rowindexset, jb, je);
            try {
                TestMatrix.check(SUB, M);
                TestMatrix.try_success("getMatrix(int[],int,int)... "); //$NON-NLS-1$
            } catch (final RuntimeException e) {
                TestMatrix.try_failure("getMatrix(int[],int,int)... ", //$NON-NLS-1$
                        "submatrix not successfully retreived"); //$NON-NLS-1$
            }
        } catch (final ArrayIndexOutOfBoundsException e) {
            TestMatrix.try_failure("getMatrix(int[],int,int)... ", //$NON-NLS-1$
                    "Unexpected ArrayIndexOutOfBoundsException"); //$NON-NLS-1$
        }

        try {
            M = B.getMatrix(badrowindexset, columnindexset);
            TestMatrix.try_failure("getMatrix(int[],int[])... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        } catch (final ArrayIndexOutOfBoundsException e) {
            try {
                M = B.getMatrix(rowindexset, badcolumnindexset);
                TestMatrix.try_failure("getMatrix(int[],int[])... ", //$NON-NLS-1$
                        "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
            } catch (final ArrayIndexOutOfBoundsException e1) {
                TestMatrix.try_success("getMatrix(int[],int[])... ArrayIndexOutOfBoundsException... "); //$NON-NLS-1$
            }
        } catch (final IllegalArgumentException e1) {
            TestMatrix.try_failure("getMatrix(int[],int[])... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        }

        try {
            M = B.getMatrix(rowindexset, columnindexset);
            try {
                TestMatrix.check(SUB, M);
                TestMatrix.try_success("getMatrix(int[],int[])... "); //$NON-NLS-1$
            } catch (final RuntimeException e) {
                TestMatrix.try_failure("getMatrix(int[],int[])... ", //$NON-NLS-1$
                        "submatrix not successfully retreived"); //$NON-NLS-1$
            }
        } catch (final ArrayIndexOutOfBoundsException e) {
            TestMatrix.try_failure("getMatrix(int[],int[])... ", //$NON-NLS-1$
                    "Unexpected ArrayIndexOutOfBoundsException"); //$NON-NLS-1$
        }

        // Various set methods:
        try {
            B.set(B.getRowDimension(), B.getColumnDimension() - 1, 0.);
            TestMatrix.try_failure("set(int,int,double)... ", //$NON-NLS-1$
                    "OutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        } catch (final ArrayIndexOutOfBoundsException e) {
            try {
                B.set(B.getRowDimension() - 1, B.getColumnDimension(), 0.);
                TestMatrix.try_failure("set(int,int,double)... ", //$NON-NLS-1$
                        "OutOfBoundsException expected but not thrown"); //$NON-NLS-1$
            } catch (final ArrayIndexOutOfBoundsException e1) {
                TestMatrix.try_success("set(int,int,double)... OutofBoundsException... "); //$NON-NLS-1$
            }
        } catch (final IllegalArgumentException e1) {
            TestMatrix.try_failure("set(int,int,double)... ", //$NON-NLS-1$
                    "OutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        }

        try {
            // test set and get
            B.set(ib, jb, 0D);
            tmp = B.get(ib, jb);
        } catch (final ArrayIndexOutOfBoundsException e1) {
            TestMatrix.try_failure("set(int,int,double)... ", //$NON-NLS-1$
                    "Unexpected ArrayIndexOutOfBoundsException"); //$NON-NLS-1$
        }

        try {
            TestMatrix.check(tmp, 0D);
            TestMatrix.try_success("set(int,int,double)... "); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("set(int,int,double)... ", //$NON-NLS-1$
                    "Matrix element not successfully set"); //$NON-NLS-1$
        }

        M = new Matrix(2, 3, 0D);
        try {
            // This must throw ArrayIndexOutOfBoundsException
            B.setMatrix(ib, ie + B.getRowDimension() + 1, jb, je, M);
            TestMatrix.try_failure("setMatrix(int,int,int,int,Matrix)... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        } catch (final ArrayIndexOutOfBoundsException e) {
            // This is expected
            try {
                B.setMatrix(ib, ie, jb, je + B.getColumnDimension() + 1, M);
                TestMatrix.try_failure("setMatrix(int,int,int,int,Matrix)... ", //$NON-NLS-1$
                        "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
            } catch (final ArrayIndexOutOfBoundsException e1) {
                TestMatrix.try_success("setMatrix(int,int,int,int,Matrix)... ArrayIndexOutOfBoundsException... "); //$NON-NLS-1$
            }
        } catch (final IllegalArgumentException e1) {
            TestMatrix.try_failure("setMatrix(int,int,int,int,Matrix)... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        }

        try {
            B.setMatrix(ib, ie, jb, je, M);
            try {
                TestMatrix.check(M.minus(B.getMatrix(ib, ie, jb, je)), M);
                TestMatrix.try_success("setMatrix(int,int,int,int,Matrix)... "); //$NON-NLS-1$
            } catch (final RuntimeException e) {
                TestMatrix.try_failure("setMatrix(int,int,int,int,Matrix)... ", //$NON-NLS-1$
                        "submatrix not successfully set"); //$NON-NLS-1$
            }
            B.setMatrix(ib, ie, jb, je, SUB);
        } catch (final ArrayIndexOutOfBoundsException e1) {
            TestMatrix.try_failure("setMatrix(int,int,int,int,Matrix)... ", //$NON-NLS-1$
                    "Unexpected ArrayIndexOutOfBoundsException"); //$NON-NLS-1$
        }

        try {
            B.setMatrix(ib, ie + B.getRowDimension() + 1, columnindexset, M);
            TestMatrix.try_failure("setMatrix(int,int,int[],Matrix)... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        } catch (final ArrayIndexOutOfBoundsException e) {
            try {
                B.setMatrix(ib, ie, badcolumnindexset, M);
                TestMatrix.try_failure("setMatrix(int,int,int[],Matrix)... ", //$NON-NLS-1$
                        "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
            } catch (final ArrayIndexOutOfBoundsException e1) {
                TestMatrix.try_success("setMatrix(int,int,int[],Matrix)... ArrayIndexOutOfBoundsException... "); //$NON-NLS-1$
            }
        } catch (final IllegalArgumentException e1) {
            TestMatrix.try_failure("setMatrix(int,int,int[],Matrix)... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        }

        try {
            B.setMatrix(ib, ie, columnindexset, M);
            try {
                TestMatrix.check(M.minus(B.getMatrix(ib, ie, columnindexset)), M);
                TestMatrix.try_success("setMatrix(int,int,int[],Matrix)... "); //$NON-NLS-1$
            } catch (final RuntimeException e) {
                TestMatrix.try_failure("setMatrix(int,int,int[],Matrix)... ", //$NON-NLS-1$
                        "submatrix not successfully set"); //$NON-NLS-1$
            }
            B.setMatrix(ib, ie, jb, je, SUB);
        } catch (final ArrayIndexOutOfBoundsException e1) {
            TestMatrix.try_failure("setMatrix(int,int,int[],Matrix)... ", //$NON-NLS-1$
                    "Unexpected ArrayIndexOutOfBoundsException"); //$NON-NLS-1$
        }

        try {
            B.setMatrix(rowindexset, jb, je + B.getColumnDimension() + 1, M);
            TestMatrix.try_failure("setMatrix(int[],int,int,Matrix)... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        } catch (final ArrayIndexOutOfBoundsException e) {
            try {
                B.setMatrix(badrowindexset, jb, je, M);
                TestMatrix.try_failure("setMatrix(int[],int,int,Matrix)... ", //$NON-NLS-1$
                        "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
            } catch (final ArrayIndexOutOfBoundsException e1) {
                TestMatrix.try_success("setMatrix(int[],int,int,Matrix)... ArrayIndexOutOfBoundsException... "); //$NON-NLS-1$
            }
        } catch (final IllegalArgumentException e1) {
            TestMatrix.try_failure("setMatrix(int[],int,int,Matrix)... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        }

        try {
            B.setMatrix(rowindexset, jb, je, M);
            try {
                TestMatrix.check(M.minus(B.getMatrix(rowindexset, jb, je)), M);
                TestMatrix.try_success("setMatrix(int[],int,int,Matrix)... "); //$NON-NLS-1$
            } catch (final RuntimeException e) {
                TestMatrix.try_failure("setMatrix(int[],int,int,Matrix)... ", //$NON-NLS-1$
                        "submatrix not successfully set"); //$NON-NLS-1$
            }
            B.setMatrix(ib, ie, jb, je, SUB);
        } catch (final ArrayIndexOutOfBoundsException e1) {
            TestMatrix.try_failure("setMatrix(int[],int,int,Matrix)... ", //$NON-NLS-1$
                    "Unexpected ArrayIndexOutOfBoundsException"); //$NON-NLS-1$
        }

        try {
            B.setMatrix(rowindexset, badcolumnindexset, M);
            TestMatrix.try_failure("setMatrix(int[],int[],Matrix)... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        } catch (final ArrayIndexOutOfBoundsException e) {
            try {
                B.setMatrix(badrowindexset, columnindexset, M);
                TestMatrix.try_failure("setMatrix(int[],int[],Matrix)... ", //$NON-NLS-1$
                        "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
            } catch (final ArrayIndexOutOfBoundsException e1) {
                TestMatrix.try_success("setMatrix(int[],int[],Matrix)... ArrayIndexOutOfBoundsException... "); //$NON-NLS-1$
            }
        } catch (final IllegalArgumentException e1) {
            TestMatrix.try_failure("setMatrix(int[],int[],Matrix)... ", //$NON-NLS-1$
                    "ArrayIndexOutOfBoundsException expected but not thrown"); //$NON-NLS-1$
        }

        try {
            B.setMatrix(rowindexset, columnindexset, M);
            try {
                TestMatrix.check(M.minus(B.getMatrix(rowindexset, columnindexset)), M);
                TestMatrix.try_success("setMatrix(int[],int[],Matrix)... "); //$NON-NLS-1$
            } catch (final RuntimeException e) {
                TestMatrix.try_failure("setMatrix(int[],int[],Matrix)... ", //$NON-NLS-1$
                        "submatrix not successfully set"); //$NON-NLS-1$
            }
        } catch (final ArrayIndexOutOfBoundsException e1) {
            TestMatrix.try_failure("setMatrix(int[],int[],Matrix)... ", //$NON-NLS-1$
                    "Unexpected ArrayIndexOutOfBoundsException"); //$NON-NLS-1$
        }

        // Array-like methods:
        // minus
        // minusEquals
        // plus
        // plusEquals
        // arrayLeftDivide
        // arrayLeftDivideEquals
        // arrayRightDivide
        // arrayRightDivideEquals
        // arrayTimes
        // arrayTimesEquals
        // uminus

        System.out.println();
        System.out.println("Testing array-like methods..."); //$NON-NLS-1$
        Matrix S = new Matrix(columnwise, nonconformld);
        Matrix R = Matrix.random(A.getRowDimension(), A.getColumnDimension());
        A = R;
        try {
            S = A.minus(S);
            TestMatrix.try_failure("minus conformance check... ", //$NON-NLS-1$
                    "nonconformance not raised"); //$NON-NLS-1$
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_success("minus conformance check... "); //$NON-NLS-1$
        }

        if (A.minus(R).norm1() != 0D) {
            TestMatrix.try_failure("minus... ", //$NON-NLS-1$
                    "(difference of identical Matrices is nonzero,\nSubsequent use of minus should be suspect)"); //$NON-NLS-1$
        } else {
            TestMatrix.try_success("minus... "); //$NON-NLS-1$
        }

        // copy Matrix
        A = new Matrix(R);
        // Test equals
        if (A.equals(R)) {
            TestMatrix.try_success("copy constructor... "); //$NON-NLS-1$
        } else {
            TestMatrix.try_success("copy constructor", //$NON-NLS-1$
                    "copy must be equal to the copied matrix... "); // }
        }

        A.minusEquals(R);
        blablabla = new Matrix(A.getRowDimension(), A.getColumnDimension());
        try {
            A.minusEquals(S);
            TestMatrix.try_failure("minusEquals conformance check... ", //$NON-NLS-1$
                    "nonconformance not raised"); //$NON-NLS-1$
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_success("minusEquals conformance check... "); //$NON-NLS-1$
        }

        if (A.minus(blablabla).norm1() != 0D) {
            TestMatrix.try_failure("minusEquals... ", //$NON-NLS-1$
                    "(difference of identical Matrices is nonzero,\nSubsequent use of minus should be suspect)"); //$NON-NLS-1$
        } else {
            TestMatrix.try_success("minusEquals... "); //$NON-NLS-1$
        }

        // copy Matrix
        A = new Matrix(R);
        if (A.equals(R)) {
            TestMatrix.try_success("copy constructor... "); //$NON-NLS-1$
        } else {
            TestMatrix.try_success("copy constructor", //$NON-NLS-1$
                    "copy must be equal to the copied matrix... "); // }
        }

        B = Matrix.random(A.getRowDimension(), A.getColumnDimension());
        C = A.minus(B);
        try {
            S = A.plus(S);
            TestMatrix.try_failure("plus conformance check... ", //$NON-NLS-1$
                    "nonconformance not raised"); //$NON-NLS-1$
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_success("plus conformance check... ", //$NON-NLS-1$
                    e.getMessage());
        }

        try {
            TestMatrix.check(C.plus(B), A);
            TestMatrix.try_success("plus... "); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("plus... ", "(C = A - B, but C + B != A)"); //$NON-NLS-1$ //$NON-NLS-2$
        }

        // test plusEquals
        C = A.minus(B);
        C.plusEquals(B);
        try {
            A.plusEquals(S);
            TestMatrix.try_failure("plusEquals conformance check... ", //$NON-NLS-1$
                    "nonconformance not raised"); //$NON-NLS-1$
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_success("plusEquals conformance check... "); //$NON-NLS-1$
        }

        try {
            TestMatrix.check(C, A);
            TestMatrix.try_success("plusEquals... "); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("plusEquals... ", //$NON-NLS-1$
                    "(C = A - B, but C = C + B != A)"); //$NON-NLS-1$
        }

        A = R.uminus();
        try {
            TestMatrix.check(A.plus(R), blablabla);
            TestMatrix.try_success("uminus... "); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("uminus... ", "(-A + A != zeros)"); //$NON-NLS-1$ //$NON-NLS-2$
        }

        A = new Matrix(R);
        if (A.equals(R)) {
            TestMatrix.try_success("copy constructor... "); //$NON-NLS-1$
        } else {
            TestMatrix.try_success("copy constructor", //$NON-NLS-1$
                    "copy must be equal to the copied matrix... "); // }
        }

        // One-Matrix
        Matrix O = new Matrix(A.getRowDimension(), A.getColumnDimension(), 1.0);
        C = A.arrayLeftDivide(R);
        try {
            S = A.arrayLeftDivide(S);
            TestMatrix.try_failure("arrayLeftDivide conformance check... ", //$NON-NLS-1$
                    "nonconformance not raised"); //$NON-NLS-1$
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_success("arrayLeftDivide conformance check... "); //$NON-NLS-1$
        }

        try {
            TestMatrix.check(C, O);
            TestMatrix.try_success("arrayLeftDivide... "); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("arrayLeftDivide... ", "(M.\\M != ones)"); //$NON-NLS-1$ //$NON-NLS-2$
        }

        try {
            A.arrayLeftDivideEquals(S);
            TestMatrix.try_failure("arrayLeftDivideEquals conformance check... ", //$NON-NLS-1$
                    "nonconformance not raised"); //$NON-NLS-1$
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_success("arrayLeftDivideEquals conformance check... "); //$NON-NLS-1$
        }

        A.arrayLeftDivideEquals(R);
        try {
            TestMatrix.check(A, O);
            TestMatrix.try_success("arrayLeftDivideEquals... "); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("arrayLeftDivideEquals... ", //$NON-NLS-1$
                    "(M.\\M != ones)"); //$NON-NLS-1$
        }

        A = new Matrix(R);
        if (A.equals(R)) {
            TestMatrix.try_success("copy constructor... "); //$NON-NLS-1$
        } else {
            TestMatrix.try_success("copy constructor", //$NON-NLS-1$
                    "copy must be equal to the copied matrix... "); // }
        }

        try {
            A.arrayRightDivide(S);
            TestMatrix.try_failure("arrayRightDivide conformance check... ", //$NON-NLS-1$
                    "nonconformance not raised"); //$NON-NLS-1$
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_success("arrayRightDivide conformance check... "); //$NON-NLS-1$ $
        }

        try {
            C = A.arrayRightDivide(R);
            TestMatrix.check(C, O);
            TestMatrix.try_success("arrayRightDivide... "); //$NON-NLS-1$
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_failure("arrayRightDivide... ", e.getMessage()); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("arrayRightDivide... ", "(M./M != ones)"); //$NON-NLS-1$ //$NON-NLS-2$
        }

        try {
            A.arrayRightDivideEquals(S);
            TestMatrix.try_failure("arrayRightDivideEquals conformance check... ", //$NON-NLS-1$
                    "nonconformance not raised"); //$NON-NLS-1$
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_success("arrayRightDivideEquals conformance check... "); //$NON-NLS-1$
        }

        try {
            A.arrayRightDivideEquals(R);
            TestMatrix.check(A, O);
            TestMatrix.try_success("arrayRightDivideEquals... "); //$NON-NLS-1$
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_failure("arrayRightDivide... ", e.getMessage()); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("arrayRightDivideEquals... ", //$NON-NLS-1$
                    "(M./M != ones)"); //$NON-NLS-1$
        }

        A = new Matrix(R);
        B = Matrix.random(A.getRowDimension(), A.getColumnDimension());
        try {
            S = A.arrayTimes(S);
            TestMatrix.try_failure("arrayTimes conformance check... ", //$NON-NLS-1$
                    "nonconformance not raised"); //$NON-NLS-1$
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_success("arrayTimes conformance check... "); //$NON-NLS-1$
        }

        C = A.arrayTimes(B);
        try {
            C.arrayRightDivideEquals(B);
            TestMatrix.check(C, A);
            TestMatrix.try_success("arrayTimes... "); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("arrayTimes... ", //$NON-NLS-1$
                    "(A = R, C = A.*B, but C./B != A)"); //$NON-NLS-1$
        }

        try {
            A.arrayTimesEquals(S);
            TestMatrix.try_failure("arrayTimesEquals conformance check... ", //$NON-NLS-1$
                    "nonconformance not raised"); //$NON-NLS-1$
        } catch (final IllegalArgumentException e) {
            TestMatrix.try_success("arrayTimesEquals conformance check... "); //$NON-NLS-1$
        }

        A.arrayTimesEquals(B);
        try {
            A.arrayRightDivideEquals(B);
            TestMatrix.check(A, R);
            TestMatrix.try_success("arrayTimesEquals... "); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("arrayTimesEquals... ", //$NON-NLS-1$
                    "(A = R, A = A.*B, but A./B != R)"); //$NON-NLS-1$
        }

        // I/O methods:
        // read
        // printjava.io.
        // serializable:
        // writeObject
        // readObject

        System.out.println();
        System.out.println("Testing I/O methods..."); //$NON-NLS-1$
        final DecimalFormat fmt = new DecimalFormat("0.0000E00"); //$NON-NLS-1$
        fmt.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));

        try (final PrintWriter file = new PrintWriter(new FileOutputStream("JamaTestMatrix.out"))) //$NON-NLS-1$
        {
            A.print(file, fmt, 10);
            TestMatrix.try_success("print()..."); //$NON-NLS-1$
        } catch (final IOException ioe) {
            TestMatrix.try_warning("print()...", //$NON-NLS-1$
                    "unexpected I/O error, unable to run print/read test;  check write permission in current directory and retry"); //$NON-NLS-1$
        }

        try (final BufferedReader br = new BufferedReader(new FileReader("JamaTestMatrix.out"))) //$NON-NLS-1$
        {
            R = Matrix.read(br);
            if (A.minus(R).norm1() < 0.001) {
                TestMatrix.try_success("read()..."); //$NON-NLS-1$
            } else {
                TestMatrix.try_failure("read()...", //$NON-NLS-1$
                        "Matrix read from file does not match Matrix printed to file"); //$NON-NLS-1$
            }
        } catch (final IOException ioe) {
            TestMatrix.try_warning("read()...", //$NON-NLS-1$
                    "unexpected I/O error, unable to run print/read test;  check write permission in current directory and retry"); //$NON-NLS-1$
        }

        R = Matrix.random(A.getRowDimension(), A.getColumnDimension());
        try (final ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(tmpname))) {
            out.writeObject(R);

            try (final ObjectInputStream sin = new ObjectInputStream(new FileInputStream(tmpname))) {
                A = (Matrix) sin.readObject();
                TestMatrix.check(A, R);
                TestMatrix.try_success("readObject(Matrix)..."); //$NON-NLS-1$
            } catch (final RuntimeException e) {
                TestMatrix.try_failure("readObject(Matrix)...", //$NON-NLS-1$
                        "Matrix not serialized correctly"); //$NON-NLS-1$
            }
        } catch (final IOException ioe) {
            TestMatrix.try_warning("writeObject()...", //$NON-NLS-1$
                    "unexpected I/O error, unable to run serialization test;  check write permission in current directory and retry"); //$NON-NLS-1$
        } catch (final Exception e) {
            TestMatrix.try_failure("writeObject(Matrix)...", //$NON-NLS-1$
                    "unexpected error in serialization test"); //$NON-NLS-1$
        }

        // LA methods:
        // transpose
        // times
        // cond
        // rank
        // det
        // trace
        // norm1
        // norm2
        // normF
        // normInf
        // solve
        // solveTranspose
        // inverse
        // chol
        // eig
        // lu
        // qr
        // svd

        System.out.println();
        System.out.println("Testing linear algebra methods..."); //$NON-NLS-1$
        A = new Matrix(columnwise, 3);
        Matrix T = new Matrix(tvals);
        T = A.transpose();
        try {
            TestMatrix.check(A.transpose(), T);
            TestMatrix.try_success("transpose..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("transpose()...", "transpose unsuccessful"); //$NON-NLS-1$ //$NON-NLS-2$
        }
        A.transpose();
        try {
            TestMatrix.check(A.norm1(), columnsummax);
            TestMatrix.try_success("norm1..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("norm1()...", "incorrect norm calculation"); //$NON-NLS-1$ //$NON-NLS-2$
        }
        try {
            TestMatrix.check(A.normInf(), rowsummax);
            TestMatrix.try_success("normInf()..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("normInf()...", //$NON-NLS-1$
                    "incorrect norm calculation"); //$NON-NLS-1$
        }
        try {
            TestMatrix.check(A.normF(), Math.sqrt(sumofsquares));
            TestMatrix.try_success("normF..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("normF()...", "incorrect norm calculation"); //$NON-NLS-1$ //$NON-NLS-2$
        }
        try {
            TestMatrix.check(A.trace(), sumofdiagonals);
            TestMatrix.try_success("trace()..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("trace()...", "incorrect trace calculation"); //$NON-NLS-1$ //$NON-NLS-2$
        }
        try {
            TestMatrix.check(A.getMatrix(0, A.getRowDimension() - 1, 0, A.getRowDimension() - 1).det(), 0.);
            TestMatrix.try_success("det()..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("det()...", //$NON-NLS-1$
                    "incorrect determinant calculation"); //$NON-NLS-1$
        }

        Matrix SQ = new Matrix(square);
        try {
            TestMatrix.check(A.times(A.transpose()), SQ);
            TestMatrix.try_success("times(Matrix)..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("times(Matrix)...", //$NON-NLS-1$
                    "incorrect Matrix-Matrix product calculation"); //$NON-NLS-1$
        }
        try {
            TestMatrix.check(A.times(0.), blablabla);
            TestMatrix.try_success("times(double)..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("times(double)...", //$NON-NLS-1$
                    "incorrect Matrix-scalar product calculation"); //$NON-NLS-1$
        }

        A = new Matrix(columnwise, 4);
        final QRDecomposition QR = A.qr();
        R = QR.getR();
        try {
            TestMatrix.check(A, QR.getQ().times(R));
            TestMatrix.try_success("QRDecomposition..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("QRDecomposition...", //$NON-NLS-1$
                    "incorrect QR decomposition calculation"); //$NON-NLS-1$
        }

        SingularValueDecomposition SVD = A.svd();
        try {
            TestMatrix.check(A, SVD.getU().times(SVD.getS().times(SVD.getV().transpose())));
            TestMatrix.try_success("SingularValueDecomposition..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("SingularValueDecomposition...", //$NON-NLS-1$
                    "incorrect singular value decomposition calculation"); //$NON-NLS-1$
        }

        final Matrix DEF = new Matrix(rankdef);
        try {
            TestMatrix.check(DEF.rank(), Math.min(DEF.getRowDimension(), DEF.getColumnDimension()) - 1);
            TestMatrix.try_success("rank()..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("rank()...", "incorrect rank calculation"); //$NON-NLS-1$ //$NON-NLS-2$
        }
        B = new Matrix(condmat);
        SVD = B.svd();
        final double[] singularvalues = SVD.getSingularValues();
        try {
            TestMatrix.check(B.cond(),
                    singularvalues[0] / singularvalues[Math.min(B.getRowDimension(), B.getColumnDimension()) - 1]);
            TestMatrix.try_success("cond()..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("cond()...", //$NON-NLS-1$
                    "incorrect condition number calculation"); //$NON-NLS-1$
        }

        final int n = A.getColumnDimension();
        A = A.getMatrix(0, n - 1, 0, n - 1);
        A.set(0, 0, 0D);

        final LUDecomposition LU = A.lu();
        try {
            TestMatrix.check(A.getMatrix(LU.getPivot(), 0, n - 1), LU.getL().times(LU.getU()));
            TestMatrix.try_success("LUDecomposition..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("LUDecomposition...", //$NON-NLS-1$
                    "incorrect LU decomposition calculation"); //$NON-NLS-1$
        }

        Matrix X = A.inverse();
        try {
            TestMatrix.check(A.times(X), Matrix.identity(3));
            TestMatrix.try_success("inverse()..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("inverse()...", //$NON-NLS-1$
                    "incorrect inverse calculation"); //$NON-NLS-1$
        }

        O = new Matrix(SUB.getRowDimension(), 1, 1.0);
        final Matrix SOL = new Matrix(sqSolution);
        SQ = SUB.getMatrix(0, SUB.getRowDimension() - 1, 0, SUB.getRowDimension() - 1);
        try {
            TestMatrix.check(SQ.solve(SOL), O);
            TestMatrix.try_success("solve()..."); //$NON-NLS-1$
        } catch (final IllegalArgumentException e1) {
            TestMatrix.try_failure("solve()...", e1.getMessage()); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("solve()...", e.getMessage()); //$NON-NLS-1$
        }

        // Cholesky Decomposition
        A = new Matrix(pvals);
        final CholeskyDecomposition Chol = A.chol();
        final Matrix L = Chol.getL();
        try {
            // Check A==L*L^T
            TestMatrix.check(A, L.times(L.transpose()));
            TestMatrix.try_success("CholeskyDecomposition..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("CholeskyDecomposition...", //$NON-NLS-1$
                    "incorrect Cholesky decomposition calculation"); //$NON-NLS-1$
        }

        try {
            X = Chol.solve(Matrix.identity(3));
            TestMatrix.check(A.times(X), Matrix.identity(3, 3));
            TestMatrix.try_success("CholeskyDecomposition solve()..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("CholeskyDecomposition solve()...", //$NON-NLS-1$
                    "incorrect Choleskydecomposition solve calculation"); //$NON-NLS-1$
        }

        EigenvalueDecomposition Eig = A.eig();
        Matrix D = Eig.getD(), V = Eig.getV();
        try {
            // Check A*V==V*D
            TestMatrix.check(A.times(V), V.times(D));
            TestMatrix.try_success("EigenvalueDecomposition (symmetric)..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("EigenvalueDecomposition (symmetric)...", //$NON-NLS-1$
                    "incorrect symmetric Eigenvalue decomposition calculation"); //$NON-NLS-1$
        }

        try {
            // Eigenvalues
            A = new Matrix(evals);
            Eig = A.eig();
            D = Eig.getD();
            V = Eig.getV();
            TestMatrix.check(A.times(V), V.times(D));
            TestMatrix.try_success("EigenvalueDecomposition (nonsymmetric)..."); //$NON-NLS-1$
        } catch (final RuntimeException e) {
            TestMatrix.try_failure("EigenvalueDecomposition (nonsymmetric)...", //$NON-NLS-1$
                    "incorrect nonsymmetric Eigenvalue decomposition calculation"); //$NON-NLS-1$
        }

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

        System.out.println();
        System.out.println("TestMatrix completed."); //$NON-NLS-1$
        System.out.println("Total errors reported: " + errorCount); //$NON-NLS-1$
        System.out.println("Total warnings reported: " + warningCount); //$NON-NLS-1$
    }
}
