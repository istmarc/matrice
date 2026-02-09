#ifndef MATRICE_MATRIX_H
#define MATRICE_MATRIX_H

#include "data_type.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

// Compute the linear offset
uint32_t compute_offset(const uint32_t strides[2], const uint32_t row, const uint32_t col);

// Matrix data type
typedef struct matrix {
  bool own_data;
  data_type type;
  uint32_t size;
  uint32_t shape[2];
  uint32_t strides[2];
  void *data;
} matrix;

// Create a new col major matrix from shape
matrix *matrix_make(data_type type, uint32_t shape[2]);

// Free a matrix
void matrix_free(matrix *mat);

// Print a matrix
void matrix_print(const matrix *mat);

// Get the rows of the matrix
uint32_t matrix_rows(const matrix* mat);

// Get the columns of the matrix
uint32_t matrix_cols(const matrix* mat);

// Get the value at row and column
void matrix_get(const matrix* mat, uint32_t row, uint32_t col, void* value);

// Get the value at the linear index
void matrix_at(const matrix* mat, uint32_t index, void* value);

// Set the matrix at row and column to value
void matrix_set(matrix* mat, uint32_t row, uint32_t col, void* value);

// Set the matrix at index to value
void matrix_set_at(matrix* mat, uint32_t index, void* value);

// Add two matrices
void matrix_add(const matrix* x, const matrix* y, matrix* z);

// Substract two matrices
void matrix_sub(const matrix* x, const matrix* y, matrix* z);

// Elementwise multiply two matrices
void matrix_mul(const matrix* x, const matrix* y, matrix* z);

// Divide two matrices
void matrix_div(const matrix* x, const matrix* y, matrix* z);

// Matrix multiply two matrices
void matrix_matmul(const matrix* x, const matrix* y, matrix* z);

// Transposed a matrix
void matrix_transpose(const matrix* x, matrix* y);

// Create a matrix of ones
matrix *matrix_ones(data_type type, uint32_t shape[2]);

// Create a matrix of zeros
matrix *matrix_zeros(data_type type, uint32_t shape[2]);

// Create a matrix filled with a single value
matrix *matrix_fill(data_type, uint32_t shape[2], void *value);

// Create a matrix of a range
matrix *matrix_arange(data_type, uint32_t shape[2], void *value);

// Test if the values of two matrices are close
bool matrix_are_close(const matrix* x, const matrix* y, double eps);

// Test if two matrices are equal
bool matrix_equal(const matrix* x, const matrix* y);

#endif
