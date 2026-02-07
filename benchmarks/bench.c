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
		float start = 1.f;
		matrix* x = matrix_arange(kfloat, shape, &start);
		matrix* y = matrix_arange(kfloat, shape, &start);
		matrix* z = matrix_make(kfloat, shape);
		matrix* zblas = matrix_make(kfloat, shape);
		begin = clock();
		matrix_matmul(x, y, z);
		end = clock();
		double tmatrice = (double)(end - begin) / CLOCKS_PER_SEC;
		uint32_t m = size;
		uint32_t n = size;
		uint32_t k = size;
		float alpha = 1.0f;
		float beta = 0.0f;
		float* casted_x = (float*) x->data->data;
		float* casted_y = (float*) y->data->data;
		float* casted_zblas = (float*) zblas->data->data;
		begin = clock();
		cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha,
			casted_x, size, casted_y, size, beta, casted_zblas, size);
		end = clock();
		double tblas = (double)(end - begin) / CLOCKS_PER_SEC;
		printf("%i;%f;%f\n", size, tmatrice, tblas);
		bool are_equal = matrix_equal(z, zblas, 1e-3);
		if (!are_equal) {
			fprintf(stderr, "Different output matrices.\n");
		}
		matrix_free(x);
		matrix_free(y);
		matrix_free(z);
		matrix_free(zblas);
		if (!are_equal) {
			exit(EXIT_FAILURE);
		}
	}
	return 0;
}
