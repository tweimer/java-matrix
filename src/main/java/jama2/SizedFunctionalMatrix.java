/**
 * 
 */
package jama2;

/**
 * Adapter to convert a {@code FunctionalMatrix} into a {@code SizedMatrix}
 * with the given size
 * 
 * @author Tobias Weimer
 * @see FunctionalMatrix
 * @see SizedMatrix
 */
public class SizedFunctionalMatrix extends RectangularMatrix {
    private final FunctionalMatrix f;
    
    /**
     * Creates a new SizedFunctionalMatrix
     * 
     * @param nRows
     *            Number of rows.
     * @param nColumns
     *            Number of colums
     * @param f
     *           The FunctionalMatrix
     */
    public SizedFunctionalMatrix(final int nRows, final int nColumns,
            final FunctionalMatrix f) {
        super(nRows, nColumns);
        this.f = f;
    }

    @Override
    public double get(final int r, final int c) {
        return f.get(r, c);
    }
}
