/// A 2D matrix implementation
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

    /// Set a value at the specified row and column
    #[inline]
    pub fn set(&mut self, row: usize, col: usize, value: f64) {
        assert!(
            row < self.rows,
            "Row index {} out of bounds (max: {})",
            row,
            self.rows
        );
        assert!(
            col < self.cols,
            "Column index {} out of bounds (max: {})",
            col,
            self.cols
        );
        self.data[row * self.cols + col] = value;
    }

    /// Create a matrix filled with zeros
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self::new(rows, cols)
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
