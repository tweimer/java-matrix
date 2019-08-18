package jama2;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

class DiagonalMatixTest {

    @Test
    @DisplayName("Test empty diagonal matrix")
    void testConstructors1() {
        final var d = new DiagonalMatrix(5);
        assertEquals(5, d.getSize(), "Size of DiagonalMatrix is wrong");
    }
    
    @Test
    @DisplayName("Test DiagonalMatrix with entries")
    void testConstructors2() {
        final var d = DiagonalMatrix.diag(1d, 2d, 3d);
        assertEquals(3, d.getSize(), "Size of DiagonalMatrix is wrong");
        assertEquals(1, d.get(0), "First diagonal entry wrong");
        assertEquals(1, d.get(0, 0), "First diagonal entry wrong (two arg getter)");
        assertEquals(0, d.get(0, 1), "Entries not on the diagonal shall be 0");
        assertEquals(2, d.get(1), "Second diagonal entry wrong");
        assertEquals(3, d.get(2), "Third diagonal entry wrong");
    }
    
    @Test
    @DisplayName("Test Constructor with lambda")
    void testConstructors3() {
        final var d = new DiagonalMatrix(5, i -> i + 1d);
        assertEquals(5, d.getSize(), "Size of DiagonalMatrix is wrong");
        assertEquals(1, d.get(0, 0), "First diagonal entry wrong (two arg getter)");
        assertEquals(0, d.get(0, 1), "Entries not on the diagonal shall be 0");
        assertEquals(1, d.get(0), "First diagonal entry wrong");
        assertEquals(2, d.get(1), "Second diagonal entry wrong");
        assertEquals(3, d.get(2), "Third diagonal entry wrong");
        assertEquals(4, d.get(3), "Forth diagonal entry wrong");
        assertEquals(5, d.get(4), "Fifth diagonal entry wrong");
        
    }
    
    @Test
    @DisplayName("Test identity")
    void testConstructors4() {
        final var d = DiagonalMatrix.identity(2);
        assertEquals(2, d.getSize(), "Size of DiagonalMatrix is wrong");
        assertEquals(1, d.get(0), "First diagonal entry wrong");
        assertEquals(1, d.get(1), "Second diagonal entry wrong");
        
    }

}
