#include "../inc/matrix.h"

#include <stdio.h>

void allocations() {
	uint32_t shape[2] = {2, 3};
	uint32_t strides[2] = {6, 1};
	// Test allocations of uninitialized matrix data
	printf("Allocations of matrix_data\n");
	matrix_data *data = matrix_data_make(kfloat, shape, strides);
	matrix_data_print(data);
	matrix_data_free(data);
	// Test allocations of uninitialized matrix
	printf("Allocations of matrix\n");
	matrix *mat = matrix_make(kfloat, shape);
	matrix_print(mat);
	matrix_free(mat);
}

void creation_ones() {
	uint32_t shape[2] = {2, 3};
	printf("Creattion matrix_ones(kint, shape)\n");
	matrix *mat_int = matrix_ones(kint, shape);
	matrix_print(mat_int);
	matrix_free(mat_int);
	printf("Creattion matrix_ones(kfloat, shape)\n");
	matrix *mat_float = matrix_ones(kfloat, shape);
	matrix_print(mat_float);
	matrix_free(mat_float);
	printf("Creattion matrix_ones(kdouble, shape)\n");
	matrix *mat_double = matrix_ones(kdouble, shape);
	matrix_print(mat_double);
	matrix_free(mat_double);
}

void creation_zeros() {
	uint32_t shape[2] = {2, 3};
	printf("Creattion matrix_zeros(kint, shape)\n");
	matrix *mat_int = matrix_zeros(kint, shape);
	matrix_print(mat_int);
	matrix_free(mat_int);
	printf("Creattion matrix_zeros(kfloat, shape)\n");
	matrix *mat_float = matrix_zeros(kfloat, shape);
	matrix_print(mat_float);
	matrix_free(mat_float);
	printf("Creattion matrix_zeros(kdouble, shape)\n");
	matrix *mat_double = matrix_zeros(kdouble, shape);
	matrix_print(mat_double);
	matrix_free(mat_double);
}

void creation_fill() {
	uint32_t shape[2] = {2, 3};
	printf("Creattion matrix_fill(kint, shape, value)\n");
	int value_int = 10;
	matrix *mat_int = matrix_fill(kint, shape, &value_int);
	matrix_print(mat_int);
	matrix_free(mat_int);
	printf("Creattion matrix_fill(kfloat, shape, value)\n");
	float value_float = 10.f;
	matrix *mat_float = matrix_fill(kfloat, shape, &value_float);
	matrix_print(mat_float);
	matrix_free(mat_float);
	printf("Creattion matrix_fill(kdouble, shape, value)\n");
	double value_double = 10.;
	matrix *mat_double = matrix_fill(kdouble, shape, &value_double);
	matrix_print(mat_double);
	matrix_free(mat_double);
}

void creation_arange() {
	uint32_t shape[2] = {2, 3};
	printf("Creattion matrix_arange(kint, shape, value)\n");
	int value_int = 1;
	matrix *mat_int = matrix_arange(kint, shape, &value_int);
	matrix_print(mat_int);
	matrix_free(mat_int);
	printf("Creattion matrix_arange(kfloat, shape, value)\n");
	float value_float = 1.f;
	matrix *mat_float = matrix_arange(kfloat, shape, &value_float);
	matrix_print(mat_float);
	matrix_free(mat_float);
	printf("Creattion matrix_arange(kdouble, shape, value)\n");
	double value_double = 1.;
	matrix *mat_double = matrix_arange(kdouble, shape, &value_double);
	matrix_print(mat_double);
	matrix_free(mat_double);
}

int main() {
	allocations();
	creation_ones();
	creation_zeros();
	creation_fill();
	creation_arange();
}
