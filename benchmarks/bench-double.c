#include "../inc/matrix.h"

#include <cblas.h>

#include <stdio.h>
#include <time.h>

int main() {
	uint32_t sizes[10] = {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
	clock_t begin, end;
	uint32_t K = 10;
	printf("size;matrice;blas\n");
	for (uint32_t i = 0; i < 10; i++) {
		uint32_t size = sizes[i];
		uint32_t shape[2] = {size, size};
		matrix* x = matrix_make(kdouble, shape);
		matrix* y = matrix_make(kdouble, shape);
		double start = 1.;
		for (uint32_t idx = 0; idx < size*size; idx++) {
			matrix_set_at(x, idx, &start);
			matrix_set_at(y, idx, &start);
			start += 1.;
		}
		matrix* z = matrix_zeros(kdouble, shape);
		matrix* zblas = matrix_zeros(kdouble, shape);
		begin = clock();
		matrix_matmul(x, y, z);
		end = clock();
		double tmatrice = (double)(end - begin) / CLOCKS_PER_SEC;
		int32_t m = (int32_t) size;
		int32_t n = (int32_t) size;
		int32_t k = (int32_t) size;
		double alpha = 1.;
		double beta = 0.;
		const double* casted_x = (const double*) x->data->data;
		const double* casted_y = (const double*) y->data->data;
		double* casted_zblas = (double*) zblas->data->data;
		int32_t ld = (int32_t)size;
		begin = clock();
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha,
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
