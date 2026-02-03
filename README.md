# matrice

Matrix data type in pure C with bindings

## Roadmap

- [x] Basic matrix type for `int`, `float` and `double`
- [] Matrix operations (`matrix_add`, `matrix_sub`, `matrix_mul`, `matrix_div`)
- [] Matrix view
- [] Matrix transpose
- [] Vector type
- [] Bindings

## Examples

Right now only a basic example like the following is working. For more examples see the `tests/` directory.

```c
#include "matrix.h"

int main() {
    uint32_t shape[2] = {2, 3};
    matrix* mat = matrix_make(kfloat, shape);
    matrix_print(mat);
    matrix_free(mat);
}
```

