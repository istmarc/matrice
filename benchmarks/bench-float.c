#include "../inc/matrix.h"

#include <cblas.h>

#include <stdio.h>
#include <time.h>

int main() {
	uint32_t sizes[10] = {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
	clock_t begin, end;
	printf("size;matrice;blas\n");
	for (uint32_t i = 0; i < 10; i++) {
		uint32_t size = sizes[i];
		uint32_t shape[2] = {size, size};
		matrix* x = matrix_make(kfloat, shape);
		matrix* y = matrix_make(kfloat, shape);
		float start = 1.f;
		for (uint32_t idx = 0; idx < size*size; idx++) {
			matrix_set_at(x, idx, &start);
			matrix_set_at(y, idx, &start);
			start += 1.f;
		}
		matrix* z = matrix_zeros(kfloat, shape);
		matrix* zblas = matrix_make(kfloat, shape);
		begin = clock();
		matrix_matmul(x, y, z);
		end = clock();
		double tmatrice = (double)(end - begin) / CLOCKS_PER_SEC;
		int32_t m = (int32_t) size;
		int32_t n = (int32_t) size;
		int32_t k = (int32_t) size;
		float alpha = 1.0f;
		float beta = 0.0f;
		const float* casted_x = (const float*) x->data->data;
		const float* casted_y = (const float*) y->data->data;
		float* casted_zblas = (float*) zblas->data->data;
		int32_t ld = (int32_t)size;
		begin = clock();
		cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha,
			casted_x, ld, casted_y, ld, beta, casted_zblas, ld);
		end = clock();
		double tblas = (double)(end - begin) / CLOCKS_PER_SEC;
		printf("%i;%f;%f\n", size, tmatrice, tblas);
		matrix_free(x);
		matrix_free(y);
		matrix_free(z);
		matrix_free(zblas);
	}
	return 0;
}
