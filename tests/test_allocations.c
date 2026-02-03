#include "../inc/matrix.h"

void allocations() {
  printf("Allocations\n");
  uint32_t shape[2] = {2, 3};
  uint32_t strides[2] = {6, 1};
  // Test allocations of matrix data
  matrix_data *data = matrix_data_make(kfloat, shape, strides);
  matrix_data_print(data);
  matrix_data_free(data);
  // Test allocations of matrix
  matrix *mat = matrix_make(kint, shape);
  matrix_print(mat);
  matrix_free(mat);
}

void creation() {
  printf("Creattion\n");
  uint32_t shape[2] = {2, 3};

  matrix *mat_int = matrix_ones(kint, shape);
  matrix_print(mat_int);
  matrix_free(mat_int);

  matrix *mat_float = matrix_ones(kfloat, shape);
  matrix_print(mat_float);
  matrix_free(mat_float);
}

int main() {
  allocations();
  creation();
}
