# matrice

Column major matrix data type in pure C with bindings

## Features
- [x] Basic matrix type for `int`, `float` and `double`
- [x] Factory functions: `matrix_ones`, `matrix_zeros`, `matrix_fill`  and `matrix_arange`
- [x] Matrix operations (`matrix_add`, `matrix_sub`, `matrix_mul`, `matrix_div`)
- [x] Matrix multiplication using auto vectorization

## Roadmap

- [] Matrix view
- [] Matrix transpose
- [] Vector type
- [] Bindings

## Data types and operations

- Data types

| matrice | numpy | pytorch | Eigen |
|---------|-------|--------|-----------|
| int     | np.int | torch.int32; torch.int | int |
| float   | np.float32 | torch.float32; torch.float | float |
| double | np.float64 | torch.float64 | double

- Empty, ones and zeros

| matrice | numpy | pytorch | Eigen |
|---------|-------|--------|-----------|
| matrix_make     | np.empty | torch.empty |  |
| matrix_zeros   | np.zeros | torch.zeros | |
| matrix_ones | np.ones | torch.ones | 

## Examples

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

- Examples of using `matrix_set` and `matrix_get`

```c
#include "matrix.h"

int main() {
   uint32_t shape[2] = {3, 3};
   matrix *mat = matrix_make(kfloat, shape);
   float x = 1.f;
   for (uint32_t j = 0; j < matrix_cols(mat); j++) {
      for (uint32_t i = 0; i < matrix_rows(mat); i++) {
         matrix_set(mat, i, j, &x);
         x += 1.f;
      }
   }
   matrix_print(mat);
   for (uint32_t j = 0; j < matrix_cols(mat); j++) {
      for (uint32_t i = 0; i < matrix_rows(mat); i++) {
         float value;
         matrix_get(mat, i, j, &value);
         printf("m[%i, %i] = %f\n", i, j, value);
      }
   }
   matrix_free(mat);
}
```

This should print:

```shell
matrix[float]
1.000000 4.000000 7.000000
2.000000 5.000000 8.000000
3.000000 6.000000 9.000000
m[0, 0] = 1.000000
m[1, 0] = 2.000000
m[2, 0] = 3.000000
m[0, 1] = 4.000000
m[1, 1] = 5.000000
m[2, 1] = 6.000000
m[0, 2] = 7.000000
m[1, 2] = 8.000000
m[2, 2] = 9.000000
```


