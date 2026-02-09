#ifndef MATRICE_VECTOR_H
#define MATRICE_VECTOR_H

#include "data_type.h"

#ifdef __cplusplus
   #include <cstdint>
   #include <cstdlib>
#else
   #include <stdbool.h>
   #include <stdint.h>
   #include <stdlib.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Vector data type
typedef struct vector {
   bool own_data;
   data_type type;
   uint32_t size;
   uint32_t capacity;
   void *data;
} vector;

// Create a new col major matrix from shape
vector *vector_make(data_type type, uint32_t size);

// Free a vector
void vector_free(vector *vec);

// Print a vector
void vector_print(const vector *vec);

// Get the size of a vector
uint32_t vector_size(const vector* vec);

// Get the value at index
void vector_at(const vector* vec, uint32_t index, void* value);

// Set the vector at index to value
void vector_set(vector* vec, uint32_t index, void* value);

// Add two vectors
void vector_add(const vector* x, const vector* y, vector* z);

// Substract two vectors
void vector_sub(const vector* x, const vector* y, vector* z);

// Elementwise multiply two vectors
void vector_mul(const vector* x, const vector* y, vector* z);

// Divide two vectors
void vector_div(const vector* x, const vector* y, vector* z);

// Create a vector of ones
vector *vector_ones(data_type type, uint32_t size);

// Create a vector of zeros
vector *vector_zeros(data_type type, uint32_t size);

// Create a vector filled with a single value
vector *vector_fill(data_type, uint32_t size, void *value);

// Create a vector of a range
vector *vector_arange(data_type, uint32_t size, void *value);

// Test if two vectors are equal
bool vector_are_close(const vector* x, const vector* y, double eps);

// Test if two vectors are equal
bool vector_equal(const vector* x, const vector* y);

#ifdef __cplusplus
}
#endif

#endif
