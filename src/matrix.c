#include "../inc/matrix.h"

#include <stdio.h>

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

void matrix_print(matrix *mat) {
	if (!mat) {
		fprintf(stderr, "Error cannot print a NULLL matrix.\n");
		exit(EXIT_FAILURE);
	}
	matrix_data_print(mat->data);
}

uint32_t matrix_rows(matrix* mat) {
	if (!mat) {
		fprintf(stderr, "Error cannot get rows of a NULLL matrix.\n");
		exit(EXIT_FAILURE);
	}
	return mat->data->shape[0];
}

uint32_t matrix_cols(matrix* mat) {
	if (!mat) {
		fprintf(stderr, "Error cannot get columns of a NULLL matrix.\n");
		exit(EXIT_FAILURE);
	}
	return mat->data->shape[1];
}

void matrix_get(matrix* mat, uint32_t row, uint32_t col, void* value) {
	if (!mat) {
		fprintf(stderr, "Error cannot get value from a NULL matrix.\n");
		exit(EXIT_FAILURE);
	}
	matrix_data_get(mat->data, row, col, value);
}

void matrix_at(matrix* mat, uint32_t index, void* value) {
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


