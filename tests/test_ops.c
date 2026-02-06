#include "../inc/matrix.h"

#include <stdio.h>

void add() {
	uint32_t shape[2] = {3, 3};
	float value = 1.f;
	matrix* x = matrix_arange(kfloat, shape, &value);
	matrix* y = matrix_arange(kfloat, shape, &value);
	matrix* z = matrix_make(kfloat, shape);
	matrix_add(x, y, z);
	matrix_print(z);
	matrix_free(x);
	matrix_free(y);
	matrix_free(z);
}

void sub() {
	uint32_t shape[2] = {3, 3};
	float value = 1.f;
	matrix* x = matrix_arange(kfloat, shape, &value);
	value = 2.f;
	matrix* y = matrix_arange(kfloat, shape, &value);
	matrix* z = matrix_make(kfloat, shape);
	matrix_sub(x, y, z);
	matrix_print(z);
	matrix_free(x);
	matrix_free(y);
	matrix_free(z);
}

void mul() {
	uint32_t shape[2] = {3, 3};
	float value = 1.f;
	matrix* x = matrix_arange(kfloat, shape, &value);
	matrix* y = matrix_arange(kfloat, shape, &value);
	matrix* z = matrix_make(kfloat, shape);
	matrix_mul(x, y, z);
	matrix_print(z);
	matrix_free(x);
	matrix_free(y);
	matrix_free(z);
}

void division() {
	uint32_t shape[2] = {3, 3};
	float value = 1.f;
	matrix* x = matrix_arange(kfloat, shape, &value);
	matrix* y = matrix_arange(kfloat, shape, &value);
	matrix* z = matrix_make(kfloat, shape);
	matrix_div(x, y, z);
	matrix_print(z);
	matrix_free(x);
	matrix_free(y);
	matrix_free(z);
}

int main() {
	add();
	sub();
	mul();
	division();
}
