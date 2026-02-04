#include "../inc/matrix.h"

#include <stdio.h>

void set_get_value() {
	uint32_t shape[2] = {3, 3};
	matrix *mat = matrix_make(kfloat, shape);
	float x = 1.f;
	for (uint32_t j = 0; j < matrix_cols(mat); j++) {
		for (uint32_t i = 0; i < matrix_rows(mat); i++) {
			matrix_set(mat, i, j, &x);
			x += 1.f;
		}
	}
	matrix_print(mat);
	for (uint32_t j = 0; j < matrix_cols(mat); j++) {
		for (uint32_t i = 0; i < matrix_rows(mat); i++) {
			float value;
			matrix_get(mat, i, j, &value);
			printf("m[%i, %i] = %f\n", i, j, value);
		}
	}
	matrix_free(mat);
}

int main() {
	set_get_value();
}
