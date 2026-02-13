#include "../inc/vector.h"

#include <stdio.h>
#include <math.h>

bool is_valid_vector_pointer(const vector* vec) {
	if (!vec) {
		return false;
	} else {
		return vec->data;
	}
}

vector *vector_make(data_type type, uint32_t size) {
	// Compute the size
	if (size == 0) {
		fprintf(stderr, "Error size must be > 0.\n");
		exit(EXIT_FAILURE);
	}
	// Allocate the vector
	vector *vec = malloc(sizeof(vector));
	if (!vec) {
		fprintf(stderr, "Error allocating vector.\n");
		exit(EXIT_FAILURE);
	}

	// Allocate the data
	void *data = NULL;
	if (type == kint || type == kint32) {
		data = (int *)malloc(size * sizeof(int));
	} else if (type == kfloat || type == kfloat32) {
		data = (float *)malloc(size * sizeof(float));
	} else if (type == kdouble || type == kfloat64) {
		data = (double *)malloc(size * sizeof(double));
	} else if (type == kint64) {
		data = (int64_t *)malloc(size * sizeof(int64_t));
	} else {
		fprintf(stderr, "Error unknown data type.\n");
		exit(EXIT_FAILURE);
	}

	if (!data) {
		fprintf(stderr, "Error allocating vector data.\n");
		exit(EXIT_FAILURE);
	}

	vec->own_data = true;
	vec->type = type;
	vec->size = size;
	vec->capacity = size;
	vec->data = data;

	return vec;
}

void vector_free(vector *vec) {
	if (vec) {
		if (vec->data)
			free(vec->data);
		free(vec);
	}
}

void vector_print(const vector *vec) {
	if (!is_valid_vector_pointer(vec)) {
		fprintf(stderr, "Error print invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if(vec->type == knone) {
		fprintf(stderr, "Error unknown data type.\n");
		exit(EXIT_FAILURE);
	}
	printf("vector of ");
	if (vec->type == kint || vec->type == kint32) {
		printf("int32");
	} else if (vec->type == kfloat || vec->type == kfloat32) {
		printf("float32");
	} else if (vec->type == kdouble || vec->type == kfloat64) {
		printf("float64");
	} else if (vec->type == kint64) {
		printf("int64");
	}
	printf("\n[");
	for (uint32_t i = 0; i < vec->size; i++) {
		if (vec->type == kint || vec->type == kint32) {
			int *casted_data = (int *)vec->data;
			printf("%i", casted_data[i]);
		} else if (vec->type == kfloat || vec->type == kfloat32) {
			float *casted_data = (float *)vec->data;
			printf("%f", casted_data[i]);
		} else if (vec->type == kdouble || vec->type == kfloat64) {
			double *casted_data = (double *)vec->data;
			printf("%f", casted_data[i]);
		} else if (vec->type == kint64) {
			int64_t *casted_data = (int64_t *)vec->data;
			printf("%li", casted_data[i]);
		}

		if (i + 1 < vec->size) {
			printf(", ");
		}
	}
	printf("]\n");
}

data_type vector_type(const vector* vec) {
	if (!is_valid_vector_pointer(vec)) {
		fprintf(stderr, "Error cannot get the type of an invalid vector pointer");
		exit(EXIT_FAILURE);
	}
	return vec->type;
}

uint32_t vector_size(const vector* vec) {
	if (!is_valid_vector_pointer(vec)) {
		fprintf(stderr, "Error cannot get the size of an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	return vec->size;
}

void vector_at(const vector* vec, uint32_t index, void* value) {
	if (!is_valid_vector_pointer(vec)) {
		fprintf(stderr, "Error cannot get value from an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (index >= vec->size) {
		fprintf(stderr, "Index out %i of range 0-%i.\n", index, vec->size);
		exit(EXIT_FAILURE);
	}
	if (vec->type == kint || vec->type == kint32) {
		int* data = (int *)vec->data;
		int* value_int = (int*) value;
		*value_int = data[index];
	} else if (vec->type == kfloat || vec->type == kfloat32) {
		float *data = (float *)vec->data;
		float* value_float = (float*) value;
		*value_float = data[index];
	} else if (vec->type == kdouble || vec->type == kfloat64) {
		double* data = (double *)vec->data;
		double* value_double = (double*) value;
		*value_double = data[index];
	} else if (vec->type == kint64) {
		int64_t* data = (int64_t *)vec->data;
		int64_t* value_int = (int64_t*) value;
		*value_int = data[index];
	} else {
		fprintf(stderr, "Unknown data type.\n");
		exit(EXIT_FAILURE);
	}
}

void vector_set(vector* vec, uint32_t index, void* value) {
	if (!is_valid_vector_pointer(vec)) {
		fprintf(stderr, "Error invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (index >= vec->size) {
		fprintf(stderr, "Index %idx out of range 0-%i.\n", index, vec->size);
	}
	if (vec->type == kint || vec->type == kint32) {
		int *ptr = (int *)vec->data;
		ptr[index] = *(int *)value;
	} else if (vec->type == kfloat || vec->type == kfloat32) {
		float *ptr = (float *)vec->data;
		ptr[index] = *(float *)value;
	} else if (vec->type == kdouble || vec->type == kfloat64) {
		double *ptr = (double *)vec->data;
		ptr[index] = *(double *)value;
	} else if (vec->type == kint64) {
		int64_t *ptr = (int64_t*)vec->data;
		ptr[index] = *(int64_t *)value;
	} else {
		fprintf(stderr, "Error unknown data type.\n");
		exit(EXIT_FAILURE);
	}
}

void vector_add(const vector* x, const vector* y, vector* z) {
	if (!is_valid_vector_pointer(x)) {
		fprintf(stderr, "Error x is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_vector_pointer(y)) {
		fprintf(stderr, "Error y is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_vector_pointer(z)) {
		fprintf(stderr, "Error z is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}

	if (x->type != y->type) {
		fprintf(stderr, "Error vector_add different types.");
		exit(EXIT_FAILURE);
	}
	if (x->type != z->type) {
		fprintf(stderr, "Error vector_add different output type.");
		exit(EXIT_FAILURE);
	}

	if (x->size != y->size) {
		fprintf(stderr, "Error different input shapes.");
		exit(EXIT_FAILURE);
	}
	if (x->size != z->size) {
		fprintf(stderr, "Error different output shape.");
		exit(EXIT_FAILURE);
	}

	uint32_t size = x->size;
	data_type type = x->type;
	if (type == kint || type == kint32) {
		int* casted_x = (int*) x->data;
		int* casted_y = (int*) y->data;
		int* casted_z = (int*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] + casted_y[i];
		}
	} else if (type == kfloat || type == kfloat32) {
		float* casted_x = (float*) x->data;
		float* casted_y = (float*) y->data;
		float* casted_z = (float*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] + casted_y[i];
		}
	} else if (type == kdouble || type == kfloat64) {
		double* casted_x = (double*) x->data;
		double* casted_y = (double*) y->data;
		double* casted_z = (double*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] + casted_y[i];
		}
	} else if (type == kint64) {
		int64_t* casted_x = (int64_t*) x->data;
		int64_t* casted_y = (int64_t*) y->data;
		int64_t* casted_z = (int64_t*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] + casted_y[i];
		}
	}
}

void vector_sub(const vector* x, const vector* y, vector* z) {
	if (!is_valid_vector_pointer(x)) {
		fprintf(stderr, "Error x is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_vector_pointer(y)) {
		fprintf(stderr, "Error y is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_vector_pointer(z)) {
		fprintf(stderr, "Error z is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}

	if (x->type != y->type) {
		fprintf(stderr, "Error vector_add different types.");
		exit(EXIT_FAILURE);
	}
	if (x->type != z->type) {
		fprintf(stderr, "Error vector_add different output type.");
		exit(EXIT_FAILURE);
	}

	if (x->size != y->size) {
		fprintf(stderr, "Error different input shapes.");
		exit(EXIT_FAILURE);
	}
	if (x->size != z->size) {
		fprintf(stderr, "Error different output shape.");
		exit(EXIT_FAILURE);
	}

	uint32_t size = x->size;
	data_type type = x->type;
	if (type == kint || type == kint32) {
		int* casted_x = (int*) x->data;
		int* casted_y = (int*) y->data;
		int* casted_z = (int*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] - casted_y[i];
		}
	} else if (type == kfloat || type == kfloat32) {
		float* casted_x = (float*) x->data;
		float* casted_y = (float*) y->data;
		float* casted_z = (float*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] - casted_y[i];
		}
	} else if (type == kdouble || type == kfloat64) {
		double* casted_x = (double*) x->data;
		double* casted_y = (double*) y->data;
		double* casted_z = (double*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] - casted_y[i];
		}
	} else if (type == kint64) {
		int64_t* casted_x = (int64_t*) x->data;
		int64_t* casted_y = (int64_t*) y->data;
		int64_t* casted_z = (int64_t*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] - casted_y[i];
		}
	}
}

void vector_mul(const vector* x, const vector* y, vector* z) {
	if (!is_valid_vector_pointer(x)) {
		fprintf(stderr, "Error x is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_vector_pointer(y)) {
		fprintf(stderr, "Error y is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_vector_pointer(z)) {
		fprintf(stderr, "Error z is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}

	if (x->type != y->type) {
		fprintf(stderr, "Error vector_add different types.");
		exit(EXIT_FAILURE);
	}
	if (x->type != z->type) {
		fprintf(stderr, "Error vector_add different output type.");
		exit(EXIT_FAILURE);
	}

	if (x->size != y->size) {
		fprintf(stderr, "Error different input shapes.");
		exit(EXIT_FAILURE);
	}
	if (x->size != z->size) {
		fprintf(stderr, "Error different output shape.");
		exit(EXIT_FAILURE);
	}

	uint32_t size = x->size;
	data_type type = x->type;
	if (type == kint || type == kint32) {
		int* casted_x = (int*) x->data;
		int* casted_y = (int*) y->data;
		int* casted_z = (int*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] * casted_y[i];
		}
	} else if (type == kfloat || type == kfloat32) {
		float* casted_x = (float*) x->data;
		float* casted_y = (float*) y->data;
		float* casted_z = (float*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] * casted_y[i];
		}
	} else if (type == kdouble || type == kfloat64) {
		double* casted_x = (double*) x->data;
		double* casted_y = (double*) y->data;
		double* casted_z = (double*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] * casted_y[i];
		}
	} else if (type == kint64) {
		int64_t* casted_x = (int64_t*) x->data;
		int64_t* casted_y = (int64_t*) y->data;
		int64_t* casted_z = (int64_t*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] * casted_y[i];
		}
	}

}

void vector_div(const vector* x, const vector* y, vector* z) {
	if (!is_valid_vector_pointer(x)) {
		fprintf(stderr, "Error x is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_vector_pointer(y)) {
		fprintf(stderr, "Error y is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_vector_pointer(z)) {
		fprintf(stderr, "Error z is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}

	if (x->type != y->type) {
		fprintf(stderr, "Error vector_add different types.");
		exit(EXIT_FAILURE);
	}
	if (x->type != z->type) {
		fprintf(stderr, "Error vector_add different output type.");
		exit(EXIT_FAILURE);
	}

	if (x->size != y->size) {
		fprintf(stderr, "Error different input shapes.");
		exit(EXIT_FAILURE);
	}
	if (x->size != z->size) {
		fprintf(stderr, "Error different output shape.");
		exit(EXIT_FAILURE);
	}

	uint32_t size = x->size;
	data_type type = x->type;
	if (type == kint || type == kint32) {
		int* casted_x = (int*) x->data;
		int* casted_y = (int*) y->data;
		int* casted_z = (int*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] / casted_y[i];
		}
	} else if (type == kfloat || type == kfloat32) {
		float* casted_x = (float*) x->data;
		float* casted_y = (float*) y->data;
		float* casted_z = (float*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] / casted_y[i];
		}
	} else if (type == kdouble || type == kfloat64) {
		double* casted_x = (double*) x->data;
		double* casted_y = (double*) y->data;
		double* casted_z = (double*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] / casted_y[i];
		}
	} else if (type == kint64) {
		int64_t* casted_x = (int64_t*) x->data;
		int64_t* casted_y = (int64_t*) y->data;
		int64_t* casted_z = (int64_t*) z->data;
		#pragma clang loop vectorize(enable)
		for (uint32_t i = 0; i < size; i++) {
			casted_z[i] = casted_x[i] / casted_y[i];
		}
	}

}

vector *vector_ones(data_type type, uint32_t size) {
	vector *vec = vector_make(type, size);
	if (type == kint || type == kint32) {
		int value = 1;
		for (uint32_t i = 0; i < vec->size; i++) {
			vector_set(vec, i, &value);
		}
	} else if (type == kfloat || type == kfloat32) {
		float value = 1.f;
		for (uint32_t i = 0; i < vec->size; i++) {
			vector_set(vec, i, &value);
		}
	} else if (type == kdouble || type == kfloat64) {
		double value = 1.;
		for (uint32_t i = 0; i < vec->size; i++) {
			vector_set(vec, i, &value);
		}
	} else if (type == kint64) {
		int64_t value = 1.;
		for (uint32_t i = 0; i < vec->size; i++) {
			vector_set(vec, i, &value);
		}
	}
	return vec;
}

vector *vector_zeros(data_type type, uint32_t size) {
	vector *vec = vector_make(type, size);
	if (type == kint || type == kint32) {
		int value = 0;
		for (uint32_t i = 0; i < vec->size; i++) {
			vector_set(vec, i, &value);
		}
	} else if (type == kfloat || type == kfloat32) {
		float value = 0.f;
		for (uint32_t i = 0; i < vec->size; i++) {
			vector_set(vec, i, &value);
		}
	} else if (type == kdouble || type == kfloat64) {
		double value = 0.;
		for (uint32_t i = 0; i < vec->size; i++) {
			vector_set(vec, i, &value);
		}
	} else if (type == kint64) {
		int64_t value = 0.;
		for (uint32_t i = 0; i < vec->size; i++) {
			vector_set(vec, i, &value);
		}
	}

	return vec;
}

vector *vector_fill(data_type type, uint32_t size, void* value) {
	vector *vec = vector_make(type, size);
	if (type == kint || type == kint32) {
		for (uint32_t i = 0; i < vec->size; i++) {
			vector_set(vec, i, value);
		}
	} else if (type == kfloat || type == kfloat32) {
		for (uint32_t i = 0; i < vec->size; i++) {
			vector_set(vec, i, value);
		}
	} else if (type == kdouble || type == kfloat64) {
		for (uint32_t i = 0; i < vec->size; i++) {
			vector_set(vec, i, value);
		}
	} else if (type == kint64) {
		for (uint32_t i = 0; i < vec->size; i++) {
			vector_set(vec, i, value);
		}
	}

	return vec;
}

vector *vector_arange(data_type type, uint32_t size, void* value) {
	vector *vec = vector_make(type, size);
	if (type == kint || type == kint32) {
		vector_set(vec, 0, value);
		for (uint32_t i = 1; i < vec->size; i++) {
			int next_value = ((int*)vec->data)[i-1] + 1;
			vector_set(vec, i, &next_value);
		}
	} else if (type == kfloat || type == kfloat32) {
		vector_set(vec, 0, value);
		for (uint32_t i = 1; i < vec->size; i++) {
			float next_value = ((float*)vec->data)[i-1] + 1.f;
			vector_set(vec, i, &next_value);
		}
	} else if (type == kdouble || type == kfloat64) {
		vector_set(vec, 0, value);
		for (uint32_t i = 1; i < vec->size; i++) {
			double next_value = ((double*)vec->data)[i-1] + 1.;
			vector_set(vec, i, &next_value);
		}
	} else if (type == kint64) {
		vector_set(vec, 0, value);
		for (uint32_t i = 1; i < vec->size; i++) {
			int64_t next_value = ((int64_t*)vec->data)[i-1] + 1.;
			vector_set(vec, i, &next_value);
		}
	}
	return vec;
}

bool vector_are_close(const vector* x, const vector *y, double eps) {
	if (!is_valid_vector_pointer(x)) {
		fprintf(stderr, "Error x is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_vector_pointer(y)) {
		fprintf(stderr, "Error y is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (x->type != y->type) {
		fprintf(stderr, "Error different vector types.\n");
		exit(EXIT_FAILURE);
	}
	if (x->size != y->size) {
		fprintf(stderr, "Eror different vector sizes.\n");
		return false;
	}
	uint32_t size = x->size;
	data_type type = x->type;
	bool ret = true;
	if (type == kint || type == kint32) {
		int* casted_x = (int*) x->data;
		int* casted_y = (int*) y->data;
		for (uint32_t i = 0; i < size; i++) {
			if (casted_x[i] != casted_y[i]) {
				ret = false;
				break;
			}
		}
	} else if (type == kfloat || type == kfloat32) {
		float* casted_x = (float*) x->data;
		float* casted_y = (float*) y->data;
		for (uint32_t i = 0; i < size; i++) {
			if (fabs(casted_x[i] - casted_y[i]) >= eps) {
				ret = false;
				break;
			}
		}
	} else if (type == kdouble || type == kfloat64) {
		double* casted_x = (double*) x->data;
		double* casted_y = (double*) y->data;
		for (uint32_t i = 0; i < size; i++) {
			if (fabs(casted_x[i] - casted_y[i]) >= eps) {
				ret = false;
				break;
			}
		}
	} else if (type == kint64) {
		int64_t* casted_x = (int64_t*) x->data;
		int64_t* casted_y = (int64_t*) y->data;
		for (uint32_t i = 0; i < size; i++) {
			if (casted_x[i] != casted_y[i]) {
				ret = false;
				break;
			}
		}
	}
	return ret;
}

bool vector_equal(const vector* x, const vector *y) {
	if (!is_valid_vector_pointer(x)) {
		fprintf(stderr, "Error x is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (!is_valid_vector_pointer(y)) {
		fprintf(stderr, "Error y is an invalid vector pointer.\n");
		exit(EXIT_FAILURE);
	}
	if (x->type != y->type) {
		fprintf(stderr, "Error different types.");
		exit(EXIT_FAILURE);
	}
	if (x->size != y->size) {
		fprintf(stderr, "Eror different sizes.\n");
		return false;
	}
	return x->data == y->data;
}

