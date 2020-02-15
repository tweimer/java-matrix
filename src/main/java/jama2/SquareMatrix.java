package jama2;

/**
 * This represents an abstract square Matrix which can be used as an anonymous class.
 * 
 * @author Tobias Weimer
 *
 */
public abstract class SquareMatrix implements SizedMatrix {
    protected final int size;
    
    /**
     * Created a square Matrix.
     * @param size
     *    Size of the square Matrix (number of rows and columns)
     */
    protected SquareMatrix(final int size) {
        this.size = size;
    }

    @Override
    public final int getRowDimension() {
        return size;
    }

    @Override
    public final int getColumnDimension() {
        return size;
    }
    

    /**
     * @return
     *  Always {@code true}
     */
    @Override
    public final boolean isSquare() {
        return true;
    }
}
