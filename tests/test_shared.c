#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>

#include "../inc/vector.h"

int main() {
	uint32_t size = 10;
	float value = 1.f;
	vector* x = vector_arange(kfloat, size, &value);
	vector* y = vector_arange(kfloat, size, &value);
	vector* z = vector_make(kfloat, size);
	vector_print(z);
	vector_free(x);
	vector_free(y);
	vector_free(z);
}
