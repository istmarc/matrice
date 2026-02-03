# matrice

Column major matrix data type in pure C with bindings

## Features
- [x] Basic matrix type for `int`, `float` and `double`
- [x] Factory functions: `matrix_ones`, `matrix_zeros`, `matrix_fill`  and `matrix_arange`

## Roadmap
- [] Matrix operations (`matrix_add`, `matrix_sub`, `matrix_mul`, `matrix_div`)
- [] Matrix view
- [] Matrix transpose
- [] Vector type
- [] Bindings

## Examples

Right now only a basic example like the following is working. For more examples see the `tests/` directory.

- Uninitialized matrix

```c
#include "matrix.h"

int main() {
    uint32_t shape[2] = {2, 3};
    matrix* mat = matrix_make(kfloat, shape);
    matrix_print(mat);
    matrix_free(mat);
}
```

- Matrix initialized with range

```c
#include "matrix.h"

int main() {
    uint32_t shape[2] = {2, 3};
    float value = 1.f;
    matrix* mat = matrix_arange(kfloat, shape, &value);
    matrix_print(mat);
    matrix_free(mat);
}
```

