use std::cell::RefCell;

// Thread-local buffer pool for alignment matrices to reduce allocations
thread_local! {
    static DYNAMIC_MATRIX_POOL: RefCell<Vec<Vec<Vec<f64>>>> = const { RefCell::new(Vec::new()) };
    static TRACEBACK_MATRIX_POOL: RefCell<Vec<Vec<Vec<u8>>>> = const { RefCell::new(Vec::new()) };
}

/// A reusable matrix buffer that automatically returns to pool when dropped
pub struct MatrixBuffer<T> {
    pub matrix: Vec<Vec<T>>,
    pool: fn(Vec<Vec<T>>),
}

impl<T> std::ops::Deref for MatrixBuffer<T> {
    type Target = Vec<Vec<T>>;

    fn deref(&self) -> &Self::Target {
        &self.matrix
    }
}

impl<T> std::ops::DerefMut for MatrixBuffer<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.matrix
    }
}

impl<T> Drop for MatrixBuffer<T> {
    fn drop(&mut self) {
        // Return matrix to pool for reuse
        let matrix = std::mem::take(&mut self.matrix);
        (self.pool)(matrix);
    }
}

/// Get a reusable dynamic programming matrix (f64)
pub fn get_dynamic_matrix(consensus_len: usize, query_len: usize) -> MatrixBuffer<f64> {
    let matrix = DYNAMIC_MATRIX_POOL.with(|pool| {
        let mut pool = pool.borrow_mut();

        // Try to find a matrix that's large enough
        if let Some(pos) = pool.iter().position(|m| {
            m.len() > consensus_len && m.first().map(|row| row.len() > query_len).unwrap_or(false)
        }) {
            let mut matrix = pool.swap_remove(pos);

            // Resize if needed (should be rare)
            matrix.resize_with(consensus_len + 1, || vec![0.0; query_len + 1]);
            for row in &mut matrix {
                row.resize(query_len + 1, 0.0);
            }

            // Clear the matrix
            for row in &mut matrix {
                for cell in row {
                    *cell = 0.0;
                }
            }

            matrix
        } else {
            // Create new matrix
            vec![vec![0.0; query_len + 1]; consensus_len + 1]
        }
    });

    MatrixBuffer {
        matrix,
        pool: return_dynamic_matrix,
    }
}

/// Get a reusable traceback matrix (u8)
pub fn get_traceback_matrix(consensus_len: usize, query_len: usize) -> MatrixBuffer<u8> {
    let matrix = TRACEBACK_MATRIX_POOL.with(|pool| {
        let mut pool = pool.borrow_mut();

        // Try to find a matrix that's large enough
        if let Some(pos) = pool.iter().position(|m| {
            m.len() > consensus_len && m.first().map(|row| row.len() > query_len).unwrap_or(false)
        }) {
            let mut matrix = pool.swap_remove(pos);

            // Resize if needed (should be rare)
            matrix.resize_with(consensus_len + 1, || vec![0; query_len + 1]);
            for row in &mut matrix {
                row.resize(query_len + 1, 0);
            }

            // Clear the matrix
            for row in &mut matrix {
                for cell in row {
                    *cell = 0;
                }
            }

            matrix
        } else {
            // Create new matrix
            vec![vec![0; query_len + 1]; consensus_len + 1]
        }
    });

    MatrixBuffer {
        matrix,
        pool: return_traceback_matrix,
    }
}

/// Return a dynamic matrix to the pool
fn return_dynamic_matrix(matrix: Vec<Vec<f64>>) {
    DYNAMIC_MATRIX_POOL.with(|pool| {
        let mut pool = pool.borrow_mut();
        // Only keep a reasonable number of matrices in pool to avoid memory bloat
        if pool.len() < 16 {
            pool.push(matrix);
        }
    });
}

/// Return a traceback matrix to the pool
fn return_traceback_matrix(matrix: Vec<Vec<u8>>) {
    TRACEBACK_MATRIX_POOL.with(|pool| {
        let mut pool = pool.borrow_mut();
        // Only keep a reasonable number of matrices in pool to avoid memory bloat
        if pool.len() < 16 {
            pool.push(matrix);
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_matrix_buffer_basic() {
        let buffer = get_dynamic_matrix(10, 20);
        assert_eq!(buffer.len(), 11);
        assert_eq!(buffer[0].len(), 21);
    }

    #[test]
    fn test_matrix_buffer_reuse() {
        // Get a matrix and let it drop
        {
            let _buffer = get_dynamic_matrix(5, 10);
        }

        // Get another matrix - should reuse the previous one
        let buffer = get_dynamic_matrix(3, 8);
        assert_eq!(buffer.len(), 4); // Should be resized appropriately
        assert_eq!(buffer[0].len(), 9);
    }
}
