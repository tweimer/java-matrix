package jama2;

/**
 * This represents the unity matrix of size n.
 * 
 * @author Tobias Weimer
 *
 */
public final class UnityMatrix extends SquareMatrix {
    /**
     * Creates a new unity matrix.
     * 
     * @param size
     *   Size of the unity matrix.
     */
    public UnityMatrix(final int size) {
        super(size);
    }

    /**
     * @return
     *  Returns 1 iff {@code i == j}, 0 otherwise (so called "Kronecker delta").
     */
    @Override
    public double get(final int r, final int c) {
        return r == c ? 1.0 : 0.0;
    }
}
