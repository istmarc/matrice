#ifndef MATRICE_MATRIX_H
#define MATRICE_MATRIX_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

// Data type
typedef enum data_type { knone, kint, kfloat, kdouble } data_type;

// Handling the matrix dat with type, shape and strides
typedef struct matrix_data {
  data_type type;
  uint32_t size;
  uint32_t shape[2];
  uint32_t strides[2];
  void *data;
} matrix_data;

// Create a new matrix_data
matrix_data *matrix_data_make(data_type type, uint32_t shape[2],
                              uint32_t strides[2]);

// Print a matrix
void matrix_data_print(matrix_data *matdata);

// Free a matrix_data
void matrix_data_free(matrix_data *matdata);

// Compute the linear offset
uint32_t compute_offset(uint32_t strides[2], uint32_t row, uint32_t col);

// Get the value at row and colum
void matrix_data_get(matrix_data *matdata, uint32_t row, uint32_t col, void* value);

// Get the value at linear index
void matrix_data_at(matrix_data *matdata, uint32_t index, void* value);

// Set the matrix at row and column to value
void matrix_data_set(matrix_data *matdata, uint32_t row, uint32_t col,
                     void *value);

// Set the matrix at index to value
void matrix_data_set_at(matrix_data *matdata, uint32_t index, void *value);

// Matrix data type
typedef struct matrix {
  bool own_data;
  matrix_data *data;
} matrix;

// Create a new col major matrix from shape
matrix *matrix_make(data_type type, uint32_t shape[2]);

// Free a matrix
void matrix_free(matrix *mat);

// Print a matrix
void matrix_print(matrix *mat);

// Get the rows of the matrix
uint32_t matrix_rows(matrix* mat);

// Get the columns of the matrix
uint32_t matrix_cols(matrix* mat);

// Get the value at row and column
void matrix_get(matrix* mat, uint32_t row, uint32_t col, void* value);

// Get the value at the linear index
void matrix_at(matrix* mat, uint32_t index, void* value);

// Set the matrix at row and column to value
void matrix_set(matrix* mat, uint32_t row, uint32_t col, void* value);

// Set the matrix ar index to value
void matrix_set_at(matrix* mat, uint32_t index, void* value);

// Add two matrices
void matrix_add(matrix* x, matrix* y, matrix* z);

// Substract two matrices
void matrix_sub(matrix* x, matrix* y, matrix* z);

// Elementwise multiply two matrices
void matrix_mul(matrix* x, matrix* y, matrix* z);

// Divide two matrices
void matrix_div(matrix* x, matrix* y, matrix* z);

// Matrix multiply two matrices
void matrix_matmul(matrix* x, matrix* y, matrix* z);

// Transposed a matrix
void matrix_transpose(matrix* x, matrix* y);

// Create a matrix of ones
matrix *matrix_ones(data_type type, uint32_t shape[2]);

// Create a matrix of zeros
matrix *matrix_zeros(data_type type, uint32_t shape[2]);

// Create a matrix filled with a single value
matrix *matrix_fill(data_type, uint32_t shape[2], void *value);

// Create a matrix of a range
matrix *matrix_arange(data_type, uint32_t shape[2], void *value);

bool matrix_equal(matrix* x, matrix* y, double eps);

#endif
