/// A simple 2D matrix implementation to replace ndarray::Array2<f64>
#[derive(Debug, Clone)]
pub struct ScoringMatrix {
    data: Vec<f64>,
    rows: usize,
    cols: usize,
}

impl ScoringMatrix {
    /// Create a new scoring matrix with the given dimensions
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            data: vec![0.0; rows * cols],
            rows,
            cols,
        }
    }

    /// Create a scoring matrix from existing data
    pub fn from_data(data: Vec<f64>, rows: usize, cols: usize) -> Self {
        assert_eq!(data.len(), rows * cols, "Data length must match rows * cols");
        Self { data, rows, cols }
    }

    /// Get a value at the specified row and column
    #[inline]
    pub fn get(&self, row: usize, col: usize) -> f64 {
        assert!(row < self.rows, "Row index {} out of bounds (max: {})", row, self.rows);
        assert!(col < self.cols, "Column index {} out of bounds (max: {})", col, self.cols);
        self.data[row * self.cols + col]
    }

    /// Set a value at the specified row and column
    #[inline]
    pub fn set(&mut self, row: usize, col: usize, value: f64) {
        assert!(row < self.rows, "Row index {} out of bounds (max: {})", row, self.rows);
        assert!(col < self.cols, "Column index {} out of bounds (max: {})", col, self.cols);
        self.data[row * self.cols + col] = value;
    }

    /// Get the number of columns
    pub fn ncols(&self) -> usize {
        self.cols
    }

    /// Create a matrix filled with zeros
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self::new(rows, cols)
    }


    /// Create a slice view of the matrix (similar to ndarray's slice functionality)
    pub fn slice(&self, row_range: std::ops::Range<usize>, col_range: std::ops::Range<usize>) -> ScoringMatrix {
        assert!(row_range.end <= self.rows, "Row range end {} out of bounds (max: {})", row_range.end, self.rows);
        assert!(col_range.end <= self.cols, "Column range end {} out of bounds (max: {})", col_range.end, self.cols);
        
        let new_rows = row_range.len();
        let new_cols = col_range.len();
        let mut new_data = Vec::with_capacity(new_rows * new_cols);
        
        for row in row_range {
            for col in col_range.clone() {
                new_data.push(self.get(row, col));
            }
        }
        
        ScoringMatrix::from_data(new_data, new_rows, new_cols)
    }
}

/// Implement indexing with [[row, col]] syntax similar to ndarray
impl std::ops::Index<[usize; 2]> for ScoringMatrix {
    type Output = f64;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        let [row, col] = index;
        &self.data[row * self.cols + col]
    }
}

/// Implement mutable indexing with [[row, col]] syntax
impl std::ops::IndexMut<[usize; 2]> for ScoringMatrix {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        let [row, col] = index;
        &mut self.data[row * self.cols + col]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_matrix_creation() {
        let matrix = ScoringMatrix::new(3, 4);
        assert_eq!(matrix.get(0, 0), 0.0);
    }

    #[test]
    fn test_matrix_indexing() {
        let mut matrix = ScoringMatrix::new(2, 2);
        matrix.set(0, 1, 5.0);
        assert_eq!(matrix.get(0, 1), 5.0);
        assert_eq!(matrix[[0, 1]], 5.0);
    }

    #[test]
    fn test_matrix_slice() {
        let mut matrix = ScoringMatrix::new(4, 4);
        matrix.set(1, 1, 10.0);
        matrix.set(2, 2, 20.0);
        
        let slice = matrix.slice(1..3, 1..3);
        assert_eq!(slice.get(0, 0), 10.0);
        assert_eq!(slice.get(1, 1), 20.0);
    }

    #[test]
    fn test_from_data() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let matrix = ScoringMatrix::from_data(data, 2, 2);
        assert_eq!(matrix.get(0, 0), 1.0);
        assert_eq!(matrix.get(0, 1), 2.0);
        assert_eq!(matrix.get(1, 0), 3.0);
        assert_eq!(matrix.get(1, 1), 4.0);
    }
}
