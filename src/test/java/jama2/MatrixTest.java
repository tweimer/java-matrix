package jama2;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;


class MatrixTest {
    final double[][] avals = { { 1D, 4D, 7D, 10D }, { 2D, 5D, 8D, 11D }, { 3D, 6D, 9D, 12D } };
    final double[] columnwise = { 1D, 2D, 3D, 4D, 5D, 6D, 7D, 8D, 9D, 10D, 11D, 12D };
    final double[] rowwise = { 1D, 4D, 7D, 10D, 2D, 5D, 8D, 11D, 3D, 6D, 9D, 12D };

    final double[][] rankdef = avals;


   /* final double[][] pvals = { { 4D, 1D, 1D }, { 1D, 2D, 3D }, { 1D, 3D, 6D } };
    final double[][] ivals = { { 1D, 0D, 0D, 0D }, { 0D, 1D, 0D, 0D }, { 0D, 0D, 1D, 0D } };
    final double[][] sqSolution = { { 13D }, { 15D } };
    final double[][] condmat = { { 1D, 3D }, { 7D, 9D } };*/


    // (raggedr,raggedc) should be out of bounds in ragged array
    //final int raggedr = 0, raggedc = 4;

    // leading dimension of intended test Matrices
    final int validld = 3;

    // leading dimension which is valid, but nonconforming
    final int nonconformld = 4;



    /**
     * Constructors and constructor-like methods:
     *   double[],
     *   int double[][] int, int
     *   int, int, double
     *   int, int, double[][]
     *   constructWithCopy(double[][])
     *   random(int,int)
     *   identity(int)
     **/
    @Test
    @DisplayName("Testing constructors and constructor-like methods...")
    void testConstructor() {

        // should trigger bad shape for construction with val
        final var invalidld = 5;

        assertThrows(IllegalArgumentException.class, () -> new Matrix(columnwise, invalidld),
                "Check that exception is thrown in packed constructor with invalid length");
        

        /** check that exception is thrown in constructWithCopy
            if input array is 'ragged' **/
        final double[][] rvals = {
                { 1D, 4D, 7D },
                { 2D, 5D, 8D, 11D },
                { 3D, 6D, 9D, 12D }
             };
        assertThrows(IllegalArgumentException.class, () -> new Matrix(rvals),
                "check that exception is thrown in default constructor if input array is 'ragged'");
        assertThrows(IllegalArgumentException.class, () -> Matrix.constructWithCopy(rvals),
                "check that exception is thrown in constructWithCopy if input array is 'ragged'");

        final var B = Matrix.constructWithCopy(avals);
        final double tmp = avals[0][0];
        avals[0][0] = 0.0;
        assertEquals(tmp, B.get(0, 0), "copy not effected... data visible outside");
        avals[0][0] = tmp;
        
        final double[][] ivals = {
                { 1D, 0D, 0D, 0D },
                { 0D, 1D, 0D, 0D },
                { 0D, 0D, 1D, 0D }
             };
        final var I = new Matrix(ivals);
        assertEquals(Matrix.identity(3, 4), I, "identity Matrix not successfully created");
    }

    @Test
    @DisplayName("Testing equals")
    void testEquals() {
        var A1 = new Matrix(avals);
        assertEquals(A1, A1, "equals not reflexive");
        assertFalse(A1.equals(null), "equals failed");

        var A2 = new Matrix(avals);
        assertEquals(A1, A2, "equals failed");
        assertEquals(A2, A1, "equals failed");

        A2 = new Matrix(A1);
        assertEquals(A1, A2, "equals failed");
        assertEquals(A2, A1, "equals failed");

        A1 = new Matrix(2, 3);
        A2 = new Matrix(2, 3);
        assertEquals(A1, A2, "equals failed");
        assertEquals(A2, A1, "equals failed");

        A1.set(0, 0, 1);
        assertNotEquals(A1, A2, "equals failed");
        assertNotEquals(A2, A1, "equals failed");

        A1 = new Matrix(2, 3);
        A2 = new Matrix(3, 3);
        assertNotEquals(A1, A2, "equals failed");
    }

    /**
     * Access Methods:
     * getColumnDimension()
     * getRowDimension()
     * getArray()
     * getArrayCopy()
     * getColumnPackedCopy()
     * getRowPackedCopy()
     * get(int,int)
     * getMatrix(int,int,int,int)
     * getMatrix(int,int,int[])
     * getMatrix(int[],int,int)
     * getMatrix(int[],int[])
     * set(int,int,double)
     * setMatrix(int,int,int,int,Matrix)
     * setMatrix(int,int,int[],Matrix)
     * setMatrix(int[],int,int,Matrix)
     * setMatrix(int[],int[],Matrix)
     **/
    @Test
    @DisplayName("Testing access methods")
    void testAccessMethods() {
        final var rows = 3;
        final var cols = 4;
        // Various get methods:
        final var B = new Matrix(avals);
        assertEquals(rows, B.getRowDimension(), "getRowDimension failed ");
        assertEquals(cols, B.getColumnDimension(), "getColumnDimension failed");

        // Check if getArray returns the same array
        var barray = B.getArray();
        assertSame(avals, barray, "getArray: didn't return the same array");

        // Check if getArrayCopy returns a copy of the array
        barray = B.getArrayCopy();
        assertNotSame(avals, barray, "getArrayCopy: data not (deep) copied");
        assertArrayEquals(avals, barray, "getArrayCopy: data not successfully (deep) copied");

        var bpacked = B.getColumnPackedCopy();
        assertArrayEquals(columnwise, bpacked, "getColumnPackedCopy: data not successfully (deep) copied by columns");

        bpacked = B.getRowPackedCopy();
        assertArrayEquals(rowwise, bpacked, "getRowPackedCopy: data not successfully (deep) copied by rows");

        assertThrows(ArrayIndexOutOfBoundsException.class, () -> B.get(B.getRowDimension(), B.getColumnDimension() - 1),
                "get(int,int): OutOfBoundsException expected but not thrown");

        assertThrows(ArrayIndexOutOfBoundsException.class, () -> B.get(B.getRowDimension() - 1, B.getColumnDimension()),
                "get(int,int): OutOfBoundsException expected but not thrown");

        assertEquals(avals[B.getRowDimension() - 1][B.getColumnDimension() - 1],
                B.get(B.getRowDimension() - 1, B.getColumnDimension() - 1),
                "get(int,int): Matrix entry (i,j) not successfully retreived");

        assertThrows(ArrayIndexOutOfBoundsException.class,
                () -> B.set(B.getRowDimension(), B.getColumnDimension() - 1, 0D),
                "set(int,int,double): OutOfBoundsException expected but not thrown (invalid row index)");
        assertThrows(ArrayIndexOutOfBoundsException.class,
                () -> B.set(B.getRowDimension() - 1, B.getColumnDimension(), 0D),
                "set(int,int,double): OutOfBoundsException expected but not thrown (invalid column index)");

        // Various set methods:
        // index ranges for sub Matrix
        final int ib = 1, ie = 2, jb = 1, je = 3;

        final double[][] subavals = { { 5D, 8D, 11D }, { 6D, 9D, 12D } };

        final int[] rowindexset = { 1, 2 }, columnindexset = { 1, 2, 3 }, badrowindexset = { 1, 3 },
                badcolumnindexset = { 1, 2, 4 };

        final var SUB = new Matrix(subavals);

        // Test getMatrix(int,int,int,int)
        assertThrows(ArrayIndexOutOfBoundsException.class, () -> B.getMatrix(ib, ie + B.getRowDimension() + 1, jb, je),
                "getMatrix(int,int,int,int): ArrayIndexOutOfBoundsException expected but not thrown (rows)");
        assertThrows(ArrayIndexOutOfBoundsException.class,
                () -> B.getMatrix(ib, ie, jb, je + B.getColumnDimension() + 1),
                "getMatrix(int,int,int,int): ArrayIndexOutOfBoundsException expected but not thrown (columns)");
        assertEquals(SUB, B.getMatrix(ib, ie, jb, je), "getMatrix: submatrix not successfully retreived");

        // Test getMatrix(int,int,int[])
        assertThrows(ArrayIndexOutOfBoundsException.class, () -> B.getMatrix(ib, ie, badcolumnindexset),
                "getMatrix(int,int,int[]): ArrayIndexOutOfBoundsException expected but not thrown (rows)");
        assertThrows(ArrayIndexOutOfBoundsException.class,
                () -> B.getMatrix(ib, ie + B.getRowDimension() + 1, columnindexset),
                "getMatrix(int,int,int[]): ArrayIndexOutOfBoundsException expected but not thrown (columns)");
        assertEquals(SUB, B.getMatrix(ib, ie, columnindexset),
                "getMatrix(int,int,int[]): submatrix not successfully retreived");

        // Test getMatrix(int[],int,int)
        assertThrows(ArrayIndexOutOfBoundsException.class, () -> B.getMatrix(badrowindexset, jb, je),
                "getMatrix(int[],int,int): ArrayIndexOutOfBoundsException expected but not thrown (rows)");
        assertThrows(ArrayIndexOutOfBoundsException.class,
                () -> B.getMatrix(rowindexset, jb, je + B.getColumnDimension() + 1),
                "getMatrix(int[],int,int): ArrayIndexOutOfBoundsException expected but not thrown (columns)");
        assertEquals(SUB, B.getMatrix(rowindexset, jb, je),
                "getMatrix(int[],int,int): submatrix not successfully retreived");

        // Test getMatrix(int[],int[])
        assertThrows(ArrayIndexOutOfBoundsException.class, () -> B.getMatrix(badrowindexset, columnindexset),
                "getMatrix(int[],int[]): ArrayIndexOutOfBoundsException expected but not thrown (rows)");
        assertThrows(ArrayIndexOutOfBoundsException.class, () -> B.getMatrix(rowindexset, badcolumnindexset),
                "getMatrix(int[],int[]): ArrayIndexOutOfBoundsException expected but not thrown (columns)");
        assertEquals(SUB, B.getMatrix(rowindexset, columnindexset),
                "getMatrix(int[],int[]): submatrix not successfully retreived");

        final Matrix M = new Matrix(2, 3);

        // Test setMatrix(int,int,int,int,Matrix)
        assertThrows(ArrayIndexOutOfBoundsException.class,
                () -> B.setMatrix(ib, ie + B.getRowDimension() + 1, jb, je, M),
                "setMatrix(int,int,int,int,Matrix): ArrayIndexOutOfBoundsException expected but not thrown (rows)");
        assertThrows(ArrayIndexOutOfBoundsException.class,
                () -> B.setMatrix(ib, ie, jb, je + B.getColumnDimension() + 1, M),
                "setMatrix(int,int,int,int,Matrix): ArrayIndexOutOfBoundsException expected but not thrown (columns)");
        B.setMatrix(ib, ie, jb, je, M);
        assertEquals(M, M.minus(B.getMatrix(ib, ie, jb, je)),
                "setMatrix(int,int,int,int,Matrix): submatrix not successfully set");

        // Test setMatrix(int,int,int[],Matrix)
        assertThrows(ArrayIndexOutOfBoundsException.class,
                () -> B.setMatrix(ib, ie + B.getRowDimension() + 1, columnindexset, M),
                "setMatrix(int,int,int[],Matrix): ArrayIndexOutOfBoundsException expected but not thrown (rows)");
        assertThrows(ArrayIndexOutOfBoundsException.class, () -> B.setMatrix(ib, ie, badcolumnindexset, M),
                "setMatrix(int,int,int[],Matrix): ArrayIndexOutOfBoundsException expected but not thrown (columns)");
        B.setMatrix(ib, ie, columnindexset, M);
        assertEquals(M, M.minus(B.getMatrix(ib, ie, columnindexset)),
                "setMatrix(int,int,int[],Matrix): submatrix not successfully set");

        // Test setMatrix(int[],int,int,Matrix)
        assertThrows(ArrayIndexOutOfBoundsException.class,
                () -> B.setMatrix(rowindexset, jb, je + B.getColumnDimension() + 1, M),
                "setMatrix(int[],int,int,Matrix): ArrayIndexOutOfBoundsException expected but not thrown (rows)");
        assertThrows(ArrayIndexOutOfBoundsException.class, () -> B.setMatrix(badrowindexset, jb, je, M),
                "setMatrix(int[],int,int,Matrix): ArrayIndexOutOfBoundsException expected but not thrown (columns)");
        B.setMatrix(rowindexset, jb, je, M);
        assertEquals(M, M.minus(B.getMatrix(rowindexset, jb, je)),
                "setMatrix(int[],int,int,Matrix): submatrix not successfully set");

        // Test setMatrix(int[],int[],Matrix)
        assertThrows(ArrayIndexOutOfBoundsException.class, () -> B.setMatrix(rowindexset, badcolumnindexset, M),
                "setMatrix(int[],int[],Matrix): ArrayIndexOutOfBoundsException expected but not thrown (rows)");
        assertThrows(ArrayIndexOutOfBoundsException.class, () -> B.setMatrix(badrowindexset, columnindexset, M),
                "setMatrix(int[],int[],Matrix): ArrayIndexOutOfBoundsException expected but not thrown (columns)");
        B.setMatrix(rowindexset, columnindexset, M);
        assertEquals(M, M.minus(B.getMatrix(rowindexset, jb, je)),
                "setMatrix(int[],int[],Matrix): submatrix not successfully set");
    }
    
    @Test
    @DisplayName("Testing linear algebra methods...")
    public void testLinearAlgebra() {
        final double columnsummax = 33D, rowsummax = 30D, sumofdiagonals = 15D, sumofsquares = 650D;
        final var A = new Matrix(columnwise, 3);
        
        final double[][] tvals =
              { { 1D, 2D, 3D },
                { 4D, 5D, 6D },
                { 7D, 8D, 9D },
                { 10D, 11D, 12D } };
        final var T = A.transpose();
        assertEquals(new Matrix(tvals), T, "transpose()... transpose unsuccessful");
        
        assertEquals(columnsummax, A.norm1(), 0D, "norm1()... incorrect norm calculation");
        assertEquals(rowsummax, A.normInf(), 0D, "normInf()... incorrect norm calculation");
        assertEquals(Math.sqrt(sumofsquares), A.normF(), 1E-13, "normF()... incorrect norm calculation");
        assertEquals(sumofdiagonals, A.trace(), 0D, "trace()... incorrect norm calculation");
        assertEquals(A.getMatrix(0, A.getRowDimension() - 1, 0, A.getRowDimension() - 1).det(), 0D, 0D,
                "det()...incorrect determinant calculation");
        

        final double[][] square = {
                { 166D, 188D, 210D },
                { 188D, 214D, 240D },
                { 210D, 240D, 270D } };
        
        final var SQ = new Matrix(square);
        assertEquals(SQ, A.times(T), "times(Matrix)... incorrect Matrix-Matrix product calculation");
        
        final double[] columnwise2 = { 2D, 4D, 6D, 8D, 10D, 12D, 14D, 16D, 18D, 20D, 22D, 24D };
        assertEquals(new Matrix(columnwise2, 3), A.times(2), "times(double)... incorrect Matrix-scalar product calculation");
        assertEquals(new Matrix(3, 4), A.times(0), "times(double)... incorrect Matrix-scalar product calculation");
    }
    
    @Test
    @DisplayName("Testing QR decomposition")
    public void testQR() {
        final var A = new Matrix(columnwise, 4);
        final var QR = A.qr();
        final var R = QR.getR();
        final var Q = QR.getQ();
        final var M = Q.times(R);
        assertEquals(M.norm1(), A.norm1(), 1E-10,
                "QRDecomposition... incorrect QR decomposition calculation");
    }
    
    @Test
    @DisplayName("Testing SVD decomposition")
    public void testSVD() {
        final var A = new Matrix(columnwise, 4);
        final var SVD = A.svd();
        final var S = SVD.getS();
        final var V = SVD.getV();
        final var U = SVD.getU();
        
        assertEquals(0, A.minus(U.times(S.times(V.transpose()))).normInf(), 1E-10,
                "SingularValueDecomposition... incorrect singular value decomposition calculation");
    }
    
    @Test
    @DisplayName("Testing Rank")
    public void testRank() {
        final var DEF = new Matrix(rankdef);
        assertEquals(DEF.rank(), Math.min(DEF.getRowDimension(), DEF.getColumnDimension()) - 1, "rank()... incorrect rank calculation");
    }

    @Test
    @DisplayName("Testing condition")
    public void testCond() {
        final double[][] condmat = { { 1D, 3D }, { 7D, 9D } };
        final var B = new Matrix(condmat);
        final var SVD = B.svd();
        final double[] singularvalues = SVD.getSingularValues();
        assertEquals(B.cond(), singularvalues[0] / singularvalues[Math.min(B.getRowDimension(), B.getColumnDimension()) - 1],
                "cond()...incorrect condition number calculation");
    }
    
    @Test
    @DisplayName("Testing LU decomposition")
    public void testLU() {
        var A = new Matrix(columnwise, 4);
        final var n = A.getColumnDimension();
        A = A.getMatrix(0, n - 1, 0, n - 1);
        A.set(0, 0, 0D);

        final var LU = A.lu();
        final var L = LU.getL();
        final var U = LU.getU();
        assertEquals(A.getMatrix(LU.getPivot(), 0, n - 1), L.times(U),
                "LUDecomposition...incorrect LU decomposition calculation");
    }

    @Test
    @DisplayName("Testing inverse")
    public void testInv() {
        var A = new Matrix(columnwise, 4);
        final var n = A.getColumnDimension();
        A = A.getMatrix(0, n - 1, 0, n - 1);
        A.set(0, 0, 0D);
        var X = A.inverse();
        assertEquals(0, Matrix.identity(3).minus(A.times(X)).normInf(), 1E-10, "inverse()... incorrect inverse calculation");
    }

}
