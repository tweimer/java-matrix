/**
 * 
 */
package jama2;

/**
 * This represents an abstract rectangular Matrix which can be used as an anonymous class.
 * 
 * @author Tobias Weimer
 *
 */
public abstract class RectangularMatrix implements SizedMatrix {
    /** Number of rows and columns. */
    protected final int nRows, nColumns;
    
    /**
     * Creates an rectangular Matrix with given number of rows and columns.
     * 
     * @param nRows
     *    Number of rows
     * @param nColumns
     *    Number of columns
     */
    protected RectangularMatrix(final int nRows, final int nColumns) {
        this.nRows = nRows;
        this.nColumns = nColumns;
    }

    @Override
    public final int getRowDimension() {
        return nRows;
    }

    @Override
    public final int getColumnDimension() {
        return nColumns;
    }
    
    @Override
    public final boolean isSquare() {
        return nRows == nColumns;
    }
}
