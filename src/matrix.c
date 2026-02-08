#include "../inc/matrix.h"

#include <stdio.h>
#include <math.h>

matrix_data *matrix_data_make(data_type type, uint32_t shape[2],
                              uint32_t strides[2]) {
	// Compute the size
	uint32_t size = shape[0] * shape[1];
	if (size == 0) {
		fprintf(stderr, "Error size must be positive.\n");
		exit(EXIT_FAILURE);
	}
	// Allocate the data
	void *data = NULL;
	if (type == kint) {
		data = (int *)malloc(size * sizeof(int));
	} else if (type == kfloat) {
		data = (float *)malloc(size * sizeof(float));
	} else if (type == kdouble) {
		data = (double *)malloc(size * sizeof(double));
	} else {
		fprintf(stderr, "Error unknown data type.\n");
		exit(EXIT_FAILURE);
	}

	if (!data) {
		fprintf(stderr, "Error allocating data.\n");
		exit(EXIT_FAILURE);
	}

	matrix_data *mat = malloc(sizeof(matrix_data));
	if (!mat) {
		fprintf(stderr, "Error allocating matrix_data.\n");
		exit(EXIT_FAILURE);
	}

	mat->data = data;
	mat->size = size;
	mat->shape[0] = shape[0];
	mat->shape[1] = shape[1];
	mat->strides[0] = strides[0];
	mat->strides[1] = strides[1];
	mat->type = type;

	return mat;
}

uint32_t compute_offset(uint32_t strides[2], uint32_t row, uint32_t col) {
	return row * strides[0] + col * strides[1];
}

void matrix_data_print(matrix_data *matdata) {
	uint32_t rows = matdata->shape[0];
	uint32_t cols = matdata->shape[1];
	printf("matrix[");
	if (matdata->type == kint) {
		printf("int");
	} else if (matdata->type == kfloat) {
		printf("float");
	} else if (matdata->type == kdouble) {
		printf("double");
	}
	printf("]\n");
	for (uint32_t i = 0; i < rows; i++) {
		for (uint32_t j = 0; j < cols; j++) {
			uint32_t offset = compute_offset(matdata->strides, i, j);
			if (matdata->type == kint) {
				int *casted_data = (int *)matdata->data;
				printf("%i", casted_data[offset]);
			} else if (matdata->type == kfloat) {
				float *casted_data = (float *)matdata->data;
				printf("%f", casted_data[offset]);
			} else if (matdata->type == kdouble) {
				double *casted_data = (double *)matdata->data;
				printf("%f", casted_data[offset]);
			} else {
				fprintf(stderr, "Error unknown data type.\n");
				exit(EXIT_FAILURE);
			}
			if (j + 1 < cols) {
				printf(" ");
			}
		}
		printf("\n");
	}
}

void matrix_data_free(matrix_data *matdata) {
	if (matdata) {
		free(matdata->data);
		free(matdata);
	}
}

void matrix_data_get(matrix_data *matdata, uint32_t row, uint32_t col, void* value) {
	if (!matdata) {
		fprintf(stderr, "Invalid pointer to matrix_data.\n");
		exit(EXIT_FAILURE);
	}
	uint32_t offset = compute_offset(matdata->strides, row, col);
	if (offset >= matdata->size) {
		fprintf(stderr, "Index out of range.\n");
		exit(EXIT_FAILURE);
	}
	if (matdata->type == kint) {
		int* data = (int *)matdata->data;
		int* value_int = (int*) value;
		*value_int = data[offset];
	} else if (matdata->type == kfloat) {
		float *data = (float *)matdata->data;
		float* value_float = (float*) value;
		*value_float = data[offset];
	} else if (matdata->type == kdouble) {
		double* data = (double *)matdata->data;
		double* value_double = (double*) value;
		*value_double = data[offset];
	} else {
		fprintf(stderr, "Unknown data type.\n");
		exit(EXIT_FAILURE);
	}
}

void matrix_data_at(matrix_data *matdata, uint32_t index, void* value) {
	if (!matdata) {
		fprintf(stderr, "Invalid pointer to matrix_data.\n");
		exit(EXIT_FAILURE);
	}
	if (index >= matdata->size) {
		fprintf(stderr, "Index out of range.\n");
		exit(EXIT_FAILURE);
	}
	if (matdata->type == kint) {
		int* data = (int *)matdata->data;
		int* value_int = (int*) value;
		*value_int = data[index];
	} else if (matdata->type == kfloat) {
		float *data = (float *)matdata->data;
		float* value_float = (float*) value;
		*value_float = data[index];
	} else if (matdata->type == kdouble) {
		double* data = (double *)matdata->data;
		double* value_double = (double*) value;
		*value_double = data[index];
	} else {
		fprintf(stderr, "Unknown data type.\n");
		exit(EXIT_FAILURE);
	}
}

void matrix_data_set(matrix_data *matdata, uint32_t row, uint32_t col,
                     void *value) {
	if (!matdata) {
		fprintf(stderr, "Invalid pointer to matrix_data.\n");
		exit(EXIT_FAILURE);
	}
	uint32_t offset = compute_offset(matdata->strides, row, col);
	if (offset >= matdata->size) {
		fprintf(stderr, "Index [row, col] out of range.\n");
	}
	if (matdata->type == kint) {
		int *ptr = (int *)matdata->data;
		ptr[offset] = *(int *)value;
	} else if (matdata->type == kfloat) {
		float *ptr = (float *)matdata->data;
		ptr[offset] = *(float *)value;
	} else if (matdata->type == kdouble) {
		double *ptr = (double *)matdata->data;
		ptr[offset] = *(double *)value;
	}
}

void matrix_data_set_at(matrix_data *matdata, uint32_t index, void *value) {
	if (!matdata) {
		fprintf(stderr, "Invalid pointer to matrix_data.\n");
		exit(EXIT_FAILURE);
	}
	if (index >= matdata->size) {
		fprintf(stderr, "Index [row, col] out of range.\n");
	}
	if (matdata->type == kint) {
		int *ptr = (int *)matdata->data;
		ptr[index] = *(int *)value;
	} else if (matdata->type == kfloat) {
		float *ptr = (float *)matdata->data;
		ptr[index] = *(float *)value;
	} else if (matdata->type == kdouble) {
		double *ptr = (double *)matdata->data;
		ptr[index] = *(double *)value;
	}
}

matrix *matrix_make(data_type type, uint32_t shape[2]) {
	matrix *mat = malloc(sizeof(matrix));
	if (!mat) {
		fprintf(stderr, "Error allocating matrix.\n");
		exit(EXIT_FAILURE);
	}
	mat->own_data = true;
	uint32_t strides[2] = {1, shape[0]};
	mat->data = matrix_data_make(type, shape, strides);
	return mat;
}

void matrix_free(matrix *mat) {
	if (mat) {
		matrix_data_free(mat->data);
		free(mat);
	}
}

void matrix_print(const matrix *mat) {
	if (!mat) {
		fprintf(stderr, "Error cannot print a NULLL matrix.\n");
		exit(EXIT_FAILURE);
	}
	matrix_data_print(mat->data);
}

uint32_t matrix_rows(const matrix* mat) {
	if (!mat) {
		fprintf(stderr, "Error cannot get rows of a NULLL matrix.\n");
		exit(EXIT_FAILURE);
	}
	return mat->data->shape[0];
}

uint32_t matrix_cols(const matrix* mat) {
	if (!mat) {
		fprintf(stderr, "Error cannot get columns of a NULLL matrix.\n");
		exit(EXIT_FAILURE);
	}
	return mat->data->shape[1];
}

void matrix_get(const matrix* mat, uint32_t row, uint32_t col, void* value) {
	if (!mat) {
		fprintf(stderr, "Error cannot get value from a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	matrix_data_get(mat->data, row, col, value);
}

void matrix_at(const matrix* mat, uint32_t index, void* value) {
	if (!mat) {
		fprintf(stderr, "Error cannot get value from a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	matrix_data_at(mat->data, index, value);
}

void matrix_set(matrix* mat, uint32_t row, uint32_t col, void* value) {
	if (!mat) {
		fprintf(stderr, "Error cannot get value from a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	matrix_data_set(mat->data, row, col, value);
}

void matrix_set_at(matrix* mat, uint32_t index, void* value) {
	if (!mat) {
		fprintf(stderr, "Error cannot get value from a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	matrix_data_set_at(mat->data, index, value);
}

void matrix_add(const matrix* x, const matrix* y, matrix* z) {
	if (!x) {
		fprintf(stderr, "Error matrix_add x is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!y) {
		fprintf(stderr, "Error matrix_add y is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!z) {
		fprintf(stderr, "Error matrix_add z is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}

	if (x->data->type != y->data->type) {
		fprintf(stderr, "Error matrix_add different types.");
		exit(EXIT_FAILURE);
	}
	if (x->data->type != z->data->type) {
		fprintf(stderr, "Error matrix_add different output type.");
		exit(EXIT_FAILURE);
	}

	if ((x->data->shape[0] != y->data->shape[0]) || (x->data->shape[1] != y->data->shape[1])) {
		fprintf(stderr, "Error matrix_add different sizes.");
		exit(EXIT_FAILURE);
	}
	if ((x->data->shape[0] != z->data->shape[0]) || (x->data->shape[1] != z->data->shape[1])) {
		fprintf(stderr, "Error matrix_add different output size.");
		exit(EXIT_FAILURE);
	}

	uint32_t size = x->data->size;
	data_type type = x->data->type;
	if (type == kint) {
		int* casted_x = (int*) x->data->data;
		int* casted_y = (int*) y->data->data;
		int* casted_z = (int*) z->data->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] + casted_y[i];
		}
	} else if (type == kfloat) {
		float* casted_x = (float*) x->data->data;
		float* casted_y = (float*) y->data->data;
		float* casted_z = (float*) z->data->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] + casted_y[i];
		}
	} else if (type == kdouble) {
		double* casted_x = (double*) x->data->data;
		double* casted_y = (double*) y->data->data;
		double* casted_z = (double*) z->data->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] + casted_y[i];
		}
	}
}

void matrix_sub(const matrix* x, const matrix* y, matrix* z) {
	if (!x) {
		fprintf(stderr, "Error matrix_sub x is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!y) {
		fprintf(stderr, "Error matrix_sub y is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!z) {
		fprintf(stderr, "Error matrix_sub z is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}

	if (x->data->type != y->data->type) {
		fprintf(stderr, "Error matrix_sub different types.");
		exit(EXIT_FAILURE);
	}
	if (x->data->type != z->data->type) {
		fprintf(stderr, "Error matrix_sub different output type.");
		exit(EXIT_FAILURE);
	}

	if ((x->data->shape[0] != y->data->shape[0]) || (x->data->shape[1] != y->data->shape[1])) {
		fprintf(stderr, "Error matrix_sub different sizes.");
		exit(EXIT_FAILURE);
	}
	if ((x->data->shape[0] != z->data->shape[0]) || (x->data->shape[1] != z->data->shape[1])) {
		fprintf(stderr, "Error matrix_sub different output size.");
		exit(EXIT_FAILURE);
	}

	uint32_t size = x->data->size;
	data_type type = x->data->type;
	if (type == kint) {
		int* casted_x = (int*) x->data->data;
		int* casted_y = (int*) y->data->data;
		int* casted_z = (int*) z->data->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] - casted_y[i];
		}
	} else if (type == kfloat) {
		float* casted_x = (float*) x->data->data;
		float* casted_y = (float*) y->data->data;
		float* casted_z = (float*) z->data->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] - casted_y[i];
		}
	} else if (type == kdouble) {
		double* casted_x = (double*) x->data->data;
		double* casted_y = (double*) y->data->data;
		double* casted_z = (double*) z->data->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] - casted_y[i];
		}
	}
}

void matrix_mul(const matrix* x, const matrix* y, matrix* z) {
	if (!x) {
		fprintf(stderr, "Error matrix_mul x is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!y) {
		fprintf(stderr, "Error matrix_mul y is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!z) {
		fprintf(stderr, "Error matrix_mul z is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}

	if (x->data->type != y->data->type) {
		fprintf(stderr, "Error matrix_mul different types.");
		exit(EXIT_FAILURE);
	}
	if (x->data->type != z->data->type) {
		fprintf(stderr, "Error matrix_mul different output type.");
		exit(EXIT_FAILURE);
	}

	if ((x->data->shape[0] != y->data->shape[0]) || (x->data->shape[1] != y->data->shape[1])) {
		fprintf(stderr, "Error matrix_mul different sizes.");
		exit(EXIT_FAILURE);
	}
	if ((x->data->shape[0] != z->data->shape[0]) || (x->data->shape[1] != z->data->shape[1])) {
		fprintf(stderr, "Error matrix_mul different output size.");
		exit(EXIT_FAILURE);
	}

	uint32_t size = x->data->size;
	data_type type = x->data->type;
	if (type == kint) {
		int* casted_x = (int*) x->data->data;
		int* casted_y = (int*) y->data->data;
		int* casted_z = (int*) z->data->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] * casted_y[i];
		}
	} else if (type == kfloat) {
		float* casted_x = (float*) x->data->data;
		float* casted_y = (float*) y->data->data;
		float* casted_z = (float*) z->data->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] * casted_y[i];
		}
	} else if (type == kdouble) {
		double* casted_x = (double*) x->data->data;
		double* casted_y = (double*) y->data->data;
		double* casted_z = (double*) z->data->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] * casted_y[i];
		}
	}
}

void matrix_div(const matrix* x, const matrix* y, matrix* z) {
	if (!x) {
		fprintf(stderr, "Error matrix_div x is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!y) {
		fprintf(stderr, "Error matrix_div y is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!z) {
		fprintf(stderr, "Error matrix_div z is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}

	if (x->data->type != y->data->type) {
		fprintf(stderr, "Error matrix_div different types.");
		exit(EXIT_FAILURE);
	}
	if (x->data->type != z->data->type) {
		fprintf(stderr, "Error matrix_div different output type.");
		exit(EXIT_FAILURE);
	}

	if ((x->data->shape[0] != y->data->shape[0]) || (x->data->shape[1] != y->data->shape[1])) {
		fprintf(stderr, "Error matrix_div different sizes.");
		exit(EXIT_FAILURE);
	}
	if ((x->data->shape[0] != z->data->shape[0]) || (x->data->shape[1] != z->data->shape[1])) {
		fprintf(stderr, "Error matrix_div different output size.");
		exit(EXIT_FAILURE);
	}

	uint32_t size = x->data->size;
	data_type type = x->data->type;
	if (type == kint) {
		int* casted_x = (int*) x->data->data;
		int* casted_y = (int*) y->data->data;
		int* casted_z = (int*) z->data->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] / casted_y[i];
		}
	} else if (type == kfloat) {
		float* casted_x = (float*) x->data->data;
		float* casted_y = (float*) y->data->data;
		float* casted_z = (float*) z->data->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] / casted_y[i];
		}
	} else if (type == kdouble) {
		double* casted_x = (double*) x->data->data;
		double* casted_y = (double*) y->data->data;
		double* casted_z = (double*) z->data->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] / casted_y[i];
		}
	}
}

void matrix_matmul(const matrix* x, const matrix* y, matrix* z) {
	if (!x) {
		fprintf(stderr, "Error matrix multiplication x is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!y) {
		fprintf(stderr, "Error matrix multiplication y is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!z) {
		fprintf(stderr, "Error matrix multiplication z is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}

	if (x->data->type != y->data->type) {
		fprintf(stderr, "Error matrix multiplication different types.");
		exit(EXIT_FAILURE);
	}
	if (x->data->type != z->data->type) {
		fprintf(stderr, "Error matrix multiplication different output type.");
		exit(EXIT_FAILURE);
	}

	if (x->data->shape[1] != y->data->shape[0]) {
		fprintf(stderr, "Error matrix multiplication incompatible input sizes.");
		exit(EXIT_FAILURE);
	}
	if ((x->data->shape[0] != z->data->shape[0]) || (y->data->shape[1] != z->data->shape[1])) {
		fprintf(stderr, "Error matrix multiplication different output size.");
		exit(EXIT_FAILURE);
	}

	// [m, p] x [p, n] -> [m, n]
	uint32_t m = x->data->shape[0];
	uint32_t p = x->data->shape[1];
	uint32_t n = y->data->shape[1];
	data_type type = x->data->type;
	uint32_t y_strides_rows = y->data->strides[0];
	uint32_t y_strides_cols = y->data->strides[1];
	uint32_t z_strides_rows = z->data->strides[0];
	uint32_t z_strides_cols = z->data->strides[1];
	if (type == kint) {
		// Transpose x
		uint32_t tshape[2] = {p, m};
		matrix* tx = matrix_make(type, tshape);
		matrix_transpose(x, tx);
		// Compute x*y in z
		int* casted_tx = (int*) tx->data->data;
		int* casted_y = (int*) y->data->data;
		int* casted_z = (int*) z->data->data;
		uint32_t tx_strides_rows = tx->data->strides[0];
		uint32_t tx_strides_cols = tx->data->strides[1];
		// Initialize z to zero
		#pragma clang loop vectorize(enable)
		for (uint32_t idx = 0; idx < m*n; idx++) {
				casted_z[idx] = 0.f;
		}
		for (uint32_t i = 0; i < m; i++) {
			for (uint32_t j = 0; j < n; j++) {
				#pragma clang loop vectorize(enable)
				for (uint32_t k = 0; k < p; k++) {
					casted_z[i*z_strides_rows + j*z_strides_cols]
						+= casted_tx[k*tx_strides_rows+ i*tx_strides_cols] * casted_y[k*y_strides_rows + j*y_strides_cols];
				}
			}
		}
		matrix_free(tx);
	} else if (type == kfloat) {
		// Transpose x
		uint32_t tshape[2] = {p, m};
		matrix* tx = matrix_make(type, tshape);
		matrix_transpose(x, tx);
		// Compute x*y in z
		float* casted_tx = (float*) tx->data->data;
		float* casted_y = (float*) y->data->data;
		float* casted_z = (float*) z->data->data;
		uint32_t tx_strides_rows = tx->data->strides[0];
		uint32_t tx_strides_cols = tx->data->strides[1];
		// Initialize z to zero
		#pragma clang loop vectorize(enable)
		for (uint32_t idx = 0; idx < m*n; idx++) {
				casted_z[idx] = 0.f;
		}
		for (uint32_t i = 0; i < m; i++) {
			for (uint32_t j = 0; j < n; j++) {
				#pragma clang loop vectorize(enable)
				for (uint32_t k = 0; k < p; k++) {
					casted_z[i*z_strides_rows + j*z_strides_cols] 
						+= casted_tx[k*tx_strides_rows + i*tx_strides_cols] * casted_y[k*y_strides_rows + j*y_strides_cols];
				}
			}
		}
		matrix_free(tx);
	} else if (type == kdouble) {
		// Transpose x
		uint32_t tshape[2] = {p, m};
		matrix* tx = matrix_make(type, tshape);
		matrix_transpose(x, tx);
		// Compute x*y in z
		double* casted_tx = (double*) tx->data->data;
		double* casted_y = (double*) y->data->data;
		double* casted_z = (double*) z->data->data;
		uint32_t tx_strides_rows = tx->data->strides[0];
		uint32_t tx_strides_cols = tx->data->strides[1];
		// Initialize z to zero
		#pragma clang loop vectorize(enable)
		for (uint32_t idx = 0; idx < m*n; idx++) {
				casted_z[idx] = 0.;
		}
		for (uint32_t i = 0; i < m; i++) {
			for (uint32_t j = 0; j < n; j++) {
				#pragma clang loop vectorize(enable)
				for (uint32_t k = 0; k < p; k++) {
					casted_z[i*z_strides_rows + j * z_strides_cols]
						+= casted_tx[k*tx_strides_rows+ i*tx_strides_cols] * casted_y[k*y_strides_rows + j*y_strides_cols];
				}
			}
		}
		matrix_free(tx);
	} else {
		fprintf(stderr, "Matrix multiplication unknown data type.\n");
		exit(EXIT_FAILURE);
	}
}

void matrix_transpose(const matrix* x, matrix* y) {
	if (!x) {
		fprintf(stderr, "Matrix transpose, x is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!y) {
		fprintf(stderr, "Matrix transpose, y is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (x->data->type != y->data->type) {
		fprintf(stderr, "Matrix transpose, different data types.\n");
		exit(EXIT_FAILURE);
	}
	if ((x->data->shape[0] != y->data->shape[1]) || (x->data->shape[1] != y->data->shape[0])) {
		fprintf(stderr, "Matrix transpose incompatible shapes.\n");
		exit(EXIT_FAILURE);
	}
	data_type type = x->data->type;
	uint32_t m = x->data->shape[0];
	uint32_t n = x->data->shape[1];
	if (type == kint) {
		for (uint32_t i = 0; i < m; i++) {
			for (uint32_t j = 0; j < n; j++) {
				int value;
				matrix_get(x, i, j, &value);
				matrix_set(y, j, i, &value);
			}
		}
	} else if (type == kfloat) {
		for (uint32_t i = 0; i < m; i++) {
			for (uint32_t j = 0; j < n; j++) {
				float value;
				matrix_get(x, i, j, &value);
				matrix_set(y, j, i, &value);
			}
		}
	} else if (type == kdouble) {
		for (uint32_t i = 0; i < m; i++) {
			for (uint32_t j = 0; j < n; j++) {
				double value;
				matrix_get(x, i, j, &value);
				matrix_set(y, j, i, &value);
			}
		}
	}
}

matrix *matrix_ones(data_type type, uint32_t shape[2]) {
	matrix *mat = matrix_make(type, shape);
	if (type == kint) {
		int value = 1;
		for (uint32_t i = 0; i < mat->data->size; i++) {
			matrix_data_set_at(mat->data, i, &value);
		}
	} else if (type == kfloat) {
		float value = 1.;
		for (uint32_t i = 0; i < mat->data->size; i++) {
			matrix_data_set_at(mat->data, i, &value);
		}
	} else if (type == kdouble) {
		double value = 1.;
		for (uint32_t i = 0; i < mat->data->size; i++) {
			matrix_data_set_at(mat->data, i, &value);
		}
	}
	return mat;
}

matrix *matrix_zeros(data_type type, uint32_t shape[2]) {
	matrix *mat = matrix_make(type, shape);
	if (type == kint) {
		int value = 0;
		for (uint32_t i = 0; i < mat->data->size; i++) {
			matrix_data_set_at(mat->data, i, &value);
		}
	} else if (type == kfloat) {
		float value = 0.f;
		for (uint32_t i = 0; i < mat->data->size; i++) {
			matrix_data_set_at(mat->data, i, &value);
		}
	} else if (type == kdouble) {
		double value = 0.;
		for (uint32_t i = 0; i < mat->data->size; i++) {
			matrix_data_set_at(mat->data, i, &value);
		}
	}
	return mat;
}

matrix *matrix_fill(data_type type, uint32_t shape[2], void* value) {
	matrix *mat = matrix_make(type, shape);
	if (type == kint) {
		for (uint32_t i = 0; i < mat->data->size; i++) {
			matrix_data_set_at(mat->data, i, value);
		}
	} else if (type == kfloat) {
		for (uint32_t i = 0; i < mat->data->size; i++) {
			matrix_data_set_at(mat->data, i, value);
		}
	} else if (type == kdouble) {
		for (uint32_t i = 0; i < mat->data->size; i++) {
			matrix_data_set_at(mat->data, i, value);
		}
	}
	return mat;
}

matrix *matrix_arange(data_type type, uint32_t shape[2], void* value) {
	matrix *mat = matrix_make(type, shape);
	if (type == kint) {
		matrix_data_set_at(mat->data, 0, value);
		for (uint32_t i = 1; i < mat->data->size; i++) {
			int next_value = ((int*)mat->data->data)[i-1] + 1;
			matrix_data_set_at(mat->data, i, &next_value);
		}
	} else if (type == kfloat) {
		matrix_data_set_at(mat->data, 0, value);
		for (uint32_t i = 1; i < mat->data->size; i++) {
			float next_value = ((float*)mat->data->data)[i-1] + 1.f;
			matrix_data_set_at(mat->data, i, &next_value);
		}
	} else if (type == kdouble) {
		matrix_data_set_at(mat->data, 0, value);
		for (uint32_t i = 1; i < mat->data->size; i++) {
			double next_value = ((double*)mat->data->data)[i-1] + 1.;
			matrix_data_set_at(mat->data, i, &next_value);
		}
	}
	return mat;
}

bool matrix_are_close(const matrix* x, const matrix *y, double eps) {
	if (!x) {
		fprintf(stderr, "Error x is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!y) {
		fprintf(stderr, "Error y is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (x->data->type != y->data->type) {
		fprintf(stderr, "Error different types.");
		exit(EXIT_FAILURE);
	}
	if ((x->data->shape[0] != y->data->shape[0]) || (x->data->shape[1] != y->data->shape[1])) {
		fprintf(stderr, "Eror different shapes.\n");
		return false;
	}
	uint32_t size = x->data->size;
	data_type type = x->data->type;
	bool ret = true;
	if (type == kint) {
		int* casted_x = (int*) x->data->data;
		int* casted_y = (int*) y->data->data;
		for (uint32_t i = 0; i < size; i++) {
			if (casted_x[i] != casted_y[i]) {
				fprintf(stderr, "Error different values at index %i, x = %i and y = %i\n", i, casted_x[i], casted_y[i]);
				ret = false;
				break;
			}
		}
	} else if (type == kfloat) {
		float* casted_x = (float*) x->data->data;
		float* casted_y = (float*) y->data->data;
		for (uint32_t i = 0; i < size; i++) {
			if (fabs(casted_x[i] - casted_y[i]) >= eps) {
				fprintf(stderr, "Error different values at index %i, x = %f and y = %f\n", i, casted_x[i], casted_y[i]);
				ret = false;
				break;
			}
		}
	} else if (type == kdouble) {
		double* casted_x = (double*) x->data->data;
		double* casted_y = (double*) y->data->data;
		for (uint32_t i = 0; i < size; i++) {
			if (fabs(casted_x[i] - casted_y[i]) >= eps) {
				fprintf(stderr, "Error different values at index %i, x = %f and y = %f\n", i, casted_x[i], casted_y[i]);
				ret = false;
				break;
			}
		}
	}
	return ret;
}

bool matrix_equal(const matrix* x, const matrix *y) {
	if (!x) {
		fprintf(stderr, "Error x is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!y) {
		fprintf(stderr, "Error y is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (x->data->type != y->data->type) {
		fprintf(stderr, "Error different types.");
		exit(EXIT_FAILURE);
	}
	if ((x->data->shape[0] != y->data->shape[0]) || (x->data->shape[1] != y->data->shape[1])) {
		fprintf(stderr, "Eror different shapes.\n");
		return false;
	}
	return x->data == y->data;
}

