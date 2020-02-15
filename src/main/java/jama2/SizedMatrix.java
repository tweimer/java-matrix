package jama2;

/**
 * This interface represents a sized matrix with row and column dimensions.
 * 
 * @author Tobias Weimer
 * @version 2.0
 *
 */
public interface SizedMatrix extends FunctionalMatrix {
   /**
    * Get row dimension.
    *
    * @return The number of rows of the Matrix.
    */
   int getRowDimension();
    
   /**
   * Get column dimension.
   *
   * @return The number of columns of the Matrix.
   */
   int getColumnDimension();
   
   /**
    * Checks of the Matrix is square
    * @return
    *   {@code true} if it is square, {@code false} otherwise.
    */
   default boolean isSquare() {
       return getRowDimension() == getColumnDimension();
   }

   /**
    * Returns a new Array of the given size with all entries.
    * @return
    *   new array
    */
   default double[][] getArrayCopy() {
       final var nRows = getRowDimension();
       final var nColumns = getColumnDimension();
       final var A = new double[nRows][nColumns];
       for (var r = 0; r < nRows; r++) {
           for (var c = 0; c < nColumns; c++) {
               A[r][c] = get(r, c);
           }
        
       }
       return A;
   }
   
   @Override
   default SizedMatrix transpose() {
      return new SizedMatrix() {
        @Override
        public double get(final int r, final int c) {
            return SizedMatrix.this.get(c, r);
        }
        
        @Override
        public int getRowDimension() {
            return SizedMatrix.this.getColumnDimension();
        }
        
        @Override
        public int getColumnDimension() {
            return SizedMatrix.this.getRowDimension();
        }
      };
   }
   
   @Override
   default SizedMatrix plus(final FunctionalMatrix f) {
       return new SizedMatrix() {
        @Override
        public double get(final int r, final int c) {
            return SizedMatrix.this.get(r, c) + f.get(r, c);
        }
        
        @Override
        public int getRowDimension() {
            return SizedMatrix.this.getRowDimension();
        }
        
        @Override
        public int getColumnDimension() {
            return SizedMatrix.this.getColumnDimension();
        }
    };
   }
}
