#include "../inc/matrix.h"

#include <stdio.h>
#include <math.h>

uint32_t compute_offset(const uint32_t strides[2], const uint32_t row, const uint32_t col) {
	return row * strides[0] + col * strides[1];
}

bool is_valid_matrix_pointer(const matrix* mat) {
	if (!mat) {
		return false;
	} else {
		return mat->data;
	}
}

matrix *matrix_make(data_type type, uint32_t shape[2]) {
	// Compute the size
	if (shape[0] == 0 || shape[1] == 0) {
		fprintf(stderr, "Error shape must be > 0.\n");
		exit(EXIT_FAILURE);
	}
	uint32_t size = shape[0] * shape[1];

	// Allocate the matrix
	matrix *mat = malloc(sizeof(matrix));
	if (!mat) {
		fprintf(stderr, "Error allocating matrix.\n");
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
		fprintf(stderr, "Error allocating matrix data.\n");
		exit(EXIT_FAILURE);
	}


	mat->own_data = true;
	mat->type = type;
	mat->size = size;
	mat->shape[0] = shape[0];
	mat->shape[1] = shape[1];
	mat->strides[0] = 1;
	mat->strides[1] = shape[0];
	mat->data = data;

	return mat;
}

void matrix_free(matrix *mat) {
	if (mat) {
		if (mat->data)
			free(mat->data);
		free(mat);
	}
}

void matrix_print(const matrix *mat) {
	if (!is_valid_matrix_pointer(mat)) {
		fprintf(stderr, "Error invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->type == knone) {
		fprintf(stderr, "Error unknown data type.\n");
		exit(EXIT_FAILURE);
	}
	uint32_t rows = mat->shape[0];
	uint32_t cols = mat->shape[1];
	printf("matrix of ");
	if (mat->type == kint) {
		printf("int");
	} else if (mat->type == kfloat) {
		printf("float");
	} else if (mat->type == kdouble) {
		printf("double");
	}
	printf("\n");
	for (uint32_t i = 0; i < rows; i++) {
		for (uint32_t j = 0; j < cols; j++) {
			uint32_t offset = compute_offset(mat->strides, i, j);
			if (mat->type == kint) {
				int *casted_data = (int *)mat->data;
				printf("%i", casted_data[offset]);
			} else if (mat->type == kfloat) {
				float *casted_data = (float *)mat->data;
				printf("%f", casted_data[offset]);
			} else if (mat->type == kdouble) {
				double *casted_data = (double *)mat->data;
				printf("%f", casted_data[offset]);
			}
			if (j + 1 < cols) {
				printf(" ");
			}
		}
		printf("\n");
	}
}

uint32_t matrix_rows(const matrix* mat) {
	if (!is_valid_matrix_pointer(mat)) {
		fprintf(stderr, "Error cannot get rows of an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	return mat->shape[0];
}

uint32_t matrix_cols(const matrix* mat) {
	if (!is_valid_matrix_pointer(mat)) {
		fprintf(stderr, "Error cannot get columns of an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	return mat->shape[1];
}

void matrix_get(const matrix* mat, uint32_t row, uint32_t col, void* value) {
	if (!is_valid_matrix_pointer(mat)) {
		fprintf(stderr, "Error cannot get value from an ivalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	uint32_t offset = compute_offset(mat->strides, row, col);
	if (offset >= mat->size) {
		fprintf(stderr, "Index out of range.\n");
		exit(EXIT_FAILURE);
	}
	if (mat->type == kint) {
		int* data = (int *)mat->data;
		int* value_int = (int*) value;
		*value_int = data[offset];
	} else if (mat->type == kfloat) {
		float *data = (float *)mat->data;
		float* value_float = (float*) value;
		*value_float = data[offset];
	} else if (mat->type == kdouble) {
		double* data = (double *)mat->data;
		double* value_double = (double*) value;
		*value_double = data[offset];
	} else {
		fprintf(stderr, "Unknown data type.\n");
		exit(EXIT_FAILURE);
	}
}

void matrix_at(const matrix* mat, uint32_t index, void* value) {
	if (!is_valid_matrix_pointer(mat)) {
		fprintf(stderr, "Error cannot get value from an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (index >= mat->size) {
		fprintf(stderr, "Index out %i of range [0, %i].\n", index, mat->size);
		exit(EXIT_FAILURE);
	}
	if (mat->type == kint) {
		int* data = (int *)mat->data;
		int* value_int = (int*) value;
		*value_int = data[index];
	} else if (mat->type == kfloat) {
		float *data = (float *)mat->data;
		float* value_float = (float*) value;
		*value_float = data[index];
	} else if (mat->type == kdouble) {
		double* data = (double *)mat->data;
		double* value_double = (double*) value;
		*value_double = data[index];
	} else {
		fprintf(stderr, "Unknown data type.\n");
		exit(EXIT_FAILURE);
	}
}

void matrix_set(matrix* mat, uint32_t row, uint32_t col, void* value) {
	if (!is_valid_matrix_pointer(mat)) {
		fprintf(stderr, "Error cannot get value from an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	uint32_t offset = compute_offset(mat->strides, row, col);
	if (offset >= mat->size) {
		fprintf(stderr, "Index [%i, %i] out of range.\n", row, col);
	}
	if (mat->type == kint) {
		int *ptr = (int *)mat->data;
		ptr[offset] = *(int *)value;
	} else if (mat->type == kfloat) {
		float *ptr = (float *)mat->data;
		ptr[offset] = *(float *)value;
	} else if (mat->type == kdouble) {
		double *ptr = (double *)mat->data;
		ptr[offset] = *(double *)value;
	}
}

void matrix_set_at(matrix* mat, uint32_t index, void* value) {
	if (!is_valid_matrix_pointer(mat)) {
		fprintf(stderr, "Error invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (index >= mat->size) {
		fprintf(stderr, "Index %idx out of range 0-%i.\n", index, mat->size);
	}
	if (mat->type == kint) {
		int *ptr = (int *)mat->data;
		ptr[index] = *(int *)value;
	} else if (mat->type == kfloat) {
		float *ptr = (float *)mat->data;
		ptr[index] = *(float *)value;
	} else if (mat->type == kdouble) {
		double *ptr = (double *)mat->data;
		ptr[index] = *(double *)value;
	} else {
		fprintf(stderr, "Error unknown data type.\n");
		exit(EXIT_FAILURE);
	}
}

void matrix_add(const matrix* x, const matrix* y, matrix* z) {
	if (!is_valid_matrix_pointer(x)) {
		fprintf(stderr, "Error x is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_matrix_pointer(y)) {
		fprintf(stderr, "Error y is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_matrix_pointer(z)) {
		fprintf(stderr, "Error z is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}

	if (x->type != y->type) {
		fprintf(stderr, "Error matrix_add different types.");
		exit(EXIT_FAILURE);
	}
	if (x->type != z->type) {
		fprintf(stderr, "Error matrix_add different output type.");
		exit(EXIT_FAILURE);
	}

	if ((x->shape[0] != y->shape[0]) || (x->shape[1] != y->shape[1])) {
		fprintf(stderr, "Error different input shapes.");
		exit(EXIT_FAILURE);
	}
	if ((x->shape[0] != z->shape[0]) || (x->shape[1] != z->shape[1])) {
		fprintf(stderr, "Error different output shape.");
		exit(EXIT_FAILURE);
	}

	uint32_t size = x->size;
	data_type type = x->type;
	if (type == kint) {
		int* casted_x = (int*) x->data;
		int* casted_y = (int*) y->data;
		int* casted_z = (int*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] + casted_y[i];
		}
	} else if (type == kfloat) {
		float* casted_x = (float*) x->data;
		float* casted_y = (float*) y->data;
		float* casted_z = (float*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] + casted_y[i];
		}
	} else if (type == kdouble) {
		double* casted_x = (double*) x->data;
		double* casted_y = (double*) y->data;
		double* casted_z = (double*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] + casted_y[i];
		}
	}
}

void matrix_sub(const matrix* x, const matrix* y, matrix* z) {
	if (!is_valid_matrix_pointer(x)) {
		fprintf(stderr, "Error x is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_matrix_pointer(y)) {
		fprintf(stderr, "Error y is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_matrix_pointer(z)) {
		fprintf(stderr, "Error z is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}

	if (x->type != y->type) {
		fprintf(stderr, "Error matrix_add different types.");
		exit(EXIT_FAILURE);
	}
	if (x->type != z->type) {
		fprintf(stderr, "Error matrix_add different output type.");
		exit(EXIT_FAILURE);
	}

	if ((x->shape[0] != y->shape[0]) || (x->shape[1] != y->shape[1])) {
		fprintf(stderr, "Error different input shapes.");
		exit(EXIT_FAILURE);
	}
	if ((x->shape[0] != z->shape[0]) || (x->shape[1] != z->shape[1])) {
		fprintf(stderr, "Error different output shape.");
		exit(EXIT_FAILURE);
	}

	uint32_t size = x->size;
	data_type type = x->type;
	if (type == kint) {
		int* casted_x = (int*) x->data;
		int* casted_y = (int*) y->data;
		int* casted_z = (int*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] - casted_y[i];
		}
	} else if (type == kfloat) {
		float* casted_x = (float*) x->data;
		float* casted_y = (float*) y->data;
		float* casted_z = (float*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] - casted_y[i];
		}
	} else if (type == kdouble) {
		double* casted_x = (double*) x->data;
		double* casted_y = (double*) y->data;
		double* casted_z = (double*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] - casted_y[i];
		}
	}
}

void matrix_mul(const matrix* x, const matrix* y, matrix* z) {
	if (!is_valid_matrix_pointer(x)) {
		fprintf(stderr, "Error x is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_matrix_pointer(y)) {
		fprintf(stderr, "Error y is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_matrix_pointer(z)) {
		fprintf(stderr, "Error z is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}

	if (x->type != y->type) {
		fprintf(stderr, "Error matrix_add different types.");
		exit(EXIT_FAILURE);
	}
	if (x->type != z->type) {
		fprintf(stderr, "Error matrix_add different output type.");
		exit(EXIT_FAILURE);
	}

	if ((x->shape[0] != y->shape[0]) || (x->shape[1] != y->shape[1])) {
		fprintf(stderr, "Error different input shapes.");
		exit(EXIT_FAILURE);
	}
	if ((x->shape[0] != z->shape[0]) || (x->shape[1] != z->shape[1])) {
		fprintf(stderr, "Error different output shape.");
		exit(EXIT_FAILURE);
	}

	uint32_t size = x->size;
	data_type type = x->type;
	if (type == kint) {
		int* casted_x = (int*) x->data;
		int* casted_y = (int*) y->data;
		int* casted_z = (int*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] * casted_y[i];
		}
	} else if (type == kfloat) {
		float* casted_x = (float*) x->data;
		float* casted_y = (float*) y->data;
		float* casted_z = (float*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] * casted_y[i];
		}
	} else if (type == kdouble) {
		double* casted_x = (double*) x->data;
		double* casted_y = (double*) y->data;
		double* casted_z = (double*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] * casted_y[i];
		}
	}
}

void matrix_div(const matrix* x, const matrix* y, matrix* z) {
	if (!is_valid_matrix_pointer(x)) {
		fprintf(stderr, "Error x is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_matrix_pointer(y)) {
		fprintf(stderr, "Error y is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_matrix_pointer(z)) {
		fprintf(stderr, "Error z is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}

	if (x->type != y->type) {
		fprintf(stderr, "Error matrix_add different types.");
		exit(EXIT_FAILURE);
	}
	if (x->type != z->type) {
		fprintf(stderr, "Error matrix_add different output type.");
		exit(EXIT_FAILURE);
	}

	if ((x->shape[0] != y->shape[0]) || (x->shape[1] != y->shape[1])) {
		fprintf(stderr, "Error different input shapes.");
		exit(EXIT_FAILURE);
	}
	if ((x->shape[0] != z->shape[0]) || (x->shape[1] != z->shape[1])) {
		fprintf(stderr, "Error different output shape.");
		exit(EXIT_FAILURE);
	}

	uint32_t size = x->size;
	data_type type = x->type;
	if (type == kint) {
		int* casted_x = (int*) x->data;
		int* casted_y = (int*) y->data;
		int* casted_z = (int*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] / casted_y[i];
		}
	} else if (type == kfloat) {
		float* casted_x = (float*) x->data;
		float* casted_y = (float*) y->data;
		float* casted_z = (float*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] / casted_y[i];
		}
	} else if (type == kdouble) {
		double* casted_x = (double*) x->data;
		double* casted_y = (double*) y->data;
		double* casted_z = (double*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] / casted_y[i];
		}
	}
}

void matrix_matmul(const matrix* x, const matrix* y, matrix* z) {
	if (!is_valid_matrix_pointer(x)) {
		fprintf(stderr, "Error x in an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_matrix_pointer(y)) {
		fprintf(stderr, "Error y in an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_matrix_pointer(z)) {
		fprintf(stderr, "Error z in an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}

	if (x->type != y->type) {
		fprintf(stderr, "Error different input types.");
		exit(EXIT_FAILURE);
	}
	if (x->type != z->type) {
		fprintf(stderr, "Error matrix different output type.");
		exit(EXIT_FAILURE);
	}

	if (x->shape[1] != y->shape[0]) {
		fprintf(stderr, "Error incompatible input sizes.");
		exit(EXIT_FAILURE);
	}
	if ((x->shape[0] != z->shape[0]) || (y->shape[1] != z->shape[1])) {
		fprintf(stderr, "Error incompatible output size.");
		exit(EXIT_FAILURE);
	}

	// [m, p] x [p, n] -> [m, n]
	uint32_t m = x->shape[0];
	uint32_t p = x->shape[1];
	uint32_t n = y->shape[1];
	data_type type = x->type;
	uint32_t y_strides_rows = y->strides[0];
	uint32_t y_strides_cols = y->strides[1];
	uint32_t z_strides_rows = z->strides[0];
	uint32_t z_strides_cols = z->strides[1];
	if (type == kint) {
		// Transpose x
		uint32_t tshape[2] = {p, m};
		matrix* tx = matrix_make(type, tshape);
		matrix_transpose(x, tx);
		// Compute x*y in z
		int* casted_tx = (int*) tx->data;
		int* casted_y = (int*) y->data;
		int* casted_z = (int*) z->data;
		uint32_t tx_strides_rows = tx->strides[0];
		uint32_t tx_strides_cols = tx->strides[1];
		// Initialize z to zero
		#pragma clang loop vectorize(enable)
		for (uint32_t idx = 0; idx < m*n; idx++) {
				casted_z[idx] = 0;
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
		float* casted_tx = (float*) tx->data;
		float* casted_y = (float*) y->data;
		float* casted_z = (float*) z->data;
		uint32_t tx_strides_rows = tx->strides[0];
		uint32_t tx_strides_cols = tx->strides[1];
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
		const double* casted_tx = (const double*) tx->data;
		const double* casted_y = (const double*) y->data;
		double* casted_z = (double*) z->data;
		uint32_t tx_strides_rows = tx->strides[0];
		uint32_t tx_strides_cols = tx->strides[1];
		#pragma clang loop vectorize(enable)
		for (uint32_t idx = 0; idx < m*n; idx++) {
				casted_z[idx] = 0.;
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
	} else {
		fprintf(stderr, "Matrix multiplication unknown data type.\n");
		exit(EXIT_FAILURE);
	}
}

void matrix_transpose(const matrix* x, matrix* y) {
	if (!is_valid_matrix_pointer(x)) {
		fprintf(stderr, "Matrix transpose x is an invalid matrix pinter.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_matrix_pointer(y)) {
		fprintf(stderr, "Matrix transpose y is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (x->type != y->type) {
		fprintf(stderr, "Matrix transpose, different data types.\n");
		exit(EXIT_FAILURE);
	}
	if ((x->shape[0] != y->shape[1]) || (x->shape[1] != y->shape[0])) {
		fprintf(stderr, "Matrix transpose incompatible shapes.\n");
		exit(EXIT_FAILURE);
	}
	data_type type = x->type;
	uint32_t m = x->shape[0];
	uint32_t n = x->shape[1];
	uint32_t x_strides_rows = x->strides[0];
	uint32_t x_strides_cols = x->strides[1];
	uint32_t y_strides_rows = y->strides[0];
	uint32_t y_strides_cols = y->strides[1];
	if (type == kint) {
		int* casted_x = (int*)x->data;
		int* casted_y = (int*)y->data;
		for (uint32_t i = 0; i < m; i++) {
			for (uint32_t j = 0; j < n; j++) {
				casted_y[j*y_strides_rows + i*y_strides_cols] = casted_x[i*x_strides_rows + j*x_strides_cols];
			}
		}
	} else if (type == kfloat) {
		float* casted_x = (float*)x->data;
		float* casted_y = (float*)y->data;
		for (uint32_t i = 0; i < m; i++) {
			for (uint32_t j = 0; j < n; j++) {
				casted_y[j*y_strides_rows + i*y_strides_cols] = casted_x[i*x_strides_rows + j*x_strides_cols];
			}
		}
	} else if (type == kdouble) {
		double* casted_x = (double*)x->data;
		double* casted_y = (double*)y->data;
		for (uint32_t i = 0; i < m; i++) {
			for (uint32_t j = 0; j < n; j++) {
				casted_y[j*y_strides_rows + i*y_strides_cols] = casted_x[i*x_strides_rows + j*x_strides_cols];
			}
		}
	}
}

matrix *matrix_ones(data_type type, uint32_t shape[2]) {
	matrix *mat = matrix_make(type, shape);
	if (type == kint) {
		int value = 1;
		for (uint32_t i = 0; i < mat->size; i++) {
			matrix_set_at(mat, i, &value);
		}
	} else if (type == kfloat) {
		float value = 1.f;
		for (uint32_t i = 0; i < mat->size; i++) {
			matrix_set_at(mat, i, &value);
		}
	} else if (type == kdouble) {
		double value = 1.;
		for (uint32_t i = 0; i < mat->size; i++) {
			matrix_set_at(mat, i, &value);
		}
	}
	return mat;
}

matrix *matrix_zeros(data_type type, uint32_t shape[2]) {
	matrix *mat = matrix_make(type, shape);
	if (type == kint) {
		int value = 0;
		for (uint32_t i = 0; i < mat->size; i++) {
			matrix_set_at(mat, i, &value);
		}
	} else if (type == kfloat) {
		float value = 0.f;
		for (uint32_t i = 0; i < mat->size; i++) {
			matrix_set_at(mat, i, &value);
		}
	} else if (type == kdouble) {
		double value = 0.;
		for (uint32_t i = 0; i < mat->size; i++) {
			matrix_set_at(mat, i, &value);
		}
	}
	return mat;
}

matrix *matrix_fill(data_type type, uint32_t shape[2], void* value) {
	matrix *mat = matrix_make(type, shape);
	if (type == kint) {
		for (uint32_t i = 0; i < mat->size; i++) {
			matrix_set_at(mat, i, value);
		}
	} else if (type == kfloat) {
		for (uint32_t i = 0; i < mat->size; i++) {
			matrix_set_at(mat, i, value);
		}
	} else if (type == kdouble) {
		for (uint32_t i = 0; i < mat->size; i++) {
			matrix_set_at(mat, i, value);
		}
	}
	return mat;
}

matrix *matrix_arange(data_type type, uint32_t shape[2], void* value) {
	matrix *mat = matrix_make(type, shape);
	if (type == kint) {
		matrix_set_at(mat, 0, value);
		for (uint32_t i = 1; i < mat->size; i++) {
			int next_value = ((int*)mat->data)[i-1] + 1;
			matrix_set_at(mat, i, &next_value);
		}
	} else if (type == kfloat) {
		matrix_set_at(mat, 0, value);
		for (uint32_t i = 1; i < mat->size; i++) {
			float next_value = ((float*)mat->data)[i-1] + 1.f;
			matrix_set_at(mat, i, &next_value);
		}
	} else if (type == kdouble) {
		matrix_set_at(mat, 0, value);
		for (uint32_t i = 1; i < mat->size; i++) {
			double next_value = ((double*)mat->data)[i-1] + 1.;
			matrix_set_at(mat, i, &next_value);
		}
	}
	return mat;
}

bool matrix_are_close(const matrix* x, const matrix *y, double eps) {
	if (!is_valid_matrix_pointer(x)) {
		fprintf(stderr, "Error x is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_matrix_pointer(y)) {
		fprintf(stderr, "Error y is an invalid matrix pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (x->type != y->type) {
		fprintf(stderr, "Error different types.");
		exit(EXIT_FAILURE);
	}
	if ((x->shape[0] != y->shape[0]) || (x->shape[1] != y->shape[1])) {
		fprintf(stderr, "Eror different shapes.\n");
		return false;
	}
	uint32_t size = x->size;
	data_type type = x->type;
	bool ret = true;
	if (type == kint) {
		int* casted_x = (int*) x->data;
		int* casted_y = (int*) y->data;
		for (uint32_t i = 0; i < size; i++) {
			if (casted_x[i] != casted_y[i]) {
				fprintf(stderr, "Error different values at index %i, x = %i and y = %i\n", i, casted_x[i], casted_y[i]);
				ret = false;
				break;
			}
		}
	} else if (type == kfloat) {
		float* casted_x = (float*) x->data;
		float* casted_y = (float*) y->data;
		for (uint32_t i = 0; i < size; i++) {
			if (fabs(casted_x[i] - casted_y[i]) >= eps) {
				fprintf(stderr, "Error different values at index %i, x = %f and y = %f\n", i, casted_x[i], casted_y[i]);
				ret = false;
				break;
			}
		}
	} else if (type == kdouble) {
		double* casted_x = (double*) x->data;
		double* casted_y = (double*) y->data;
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
	if (!is_valid_matrix_pointer(x)) {
		fprintf(stderr, "Error x is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_matrix_pointer(y)) {
		fprintf(stderr, "Error y is a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	if (x->type != y->type) {
		fprintf(stderr, "Error different types.");
		exit(EXIT_FAILURE);
	}
	if ((x->shape[0] != y->shape[0]) || (x->shape[1] != y->shape[1])) {
		fprintf(stderr, "Eror different shapes.\n");
		return false;
	}
	return x->data == y->data;
}

