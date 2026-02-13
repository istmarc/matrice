# matrice

Column major matrix data type in pure C with bindings

## Features
- [x] Basic vector and matrix type for `int` (`int32`), `float` (`float32`) and `double` (`float64`).
- [x] Factory functions: `vector_ones`,`vector_zeros`, `vector_fill`, `vector_arange`, `matrix_ones`, `matrix_zeros`, `matrix_fill` and `matrix_arange`.
- [x] Elementwise vector and matrix operations `vector_add`, `vector_sub`, `vector_mul`, `vector_div`, `matrix_add`, `matrix_sub`, `matrix_mul`, and `matrix_div`.
- [x] Matrix multiplication using auto vectorization
- [x] Bindings for C++ with classes.
- [x] Bindings for OCaml

## Roadmap
- [] Support for `int64` data type
- [] Matrix view
- [] Matrix transpose
- [] Bindings for Python.
- [] Matrix decomposition methods : QR, LU, Cholesky
- [] Determinant
- [] Inverse of a matrix

## Data types and operations

- Data types

| matrice | numpy | pytorch | Eigen |
|---------|-------|--------|-----------|
| `int`; `int32` | np.int | torch.int32; torch.int | int |
| `float`; `float32` | np.float32 | torch.float32; torch.float | float |
| `double`; `float64` | np.float64 | torch.float64 | double

- Empty, ones and zeros

| matrice | numpy | pytorch | Eigen |
|---------|-------|--------|-----------|
| `matrix_make(type, shape)`;`vector_make(type, size)` | np.empty | torch.empty |  |
| `matrix_zeros(type, shape)`;`vector_zeros(type, size)` | np.zeros | torch.zeros | |
| `matrix_ones(type, shape)`;`vector_ones(type, size)`| np.ones | torch.ones | 

## Build and install

The following command will build and install the C shared library and header files.

```shell
make main
make install
```

The C++ bindings can be installed by running

```shell
make install-cpp
```

The OCaml bindings can be installed by typing

```shell
cd bindings/ocaml/matrice
dune build
dune install
```

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

- Example of using `matrix_set` and `matrix_get`

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

## C++ Bindings

In C++, the classes `matrice::vector` and `matrice::matrix` are defined.

A basic example is the following:

```cpp
#include <matrice/matrice.hxx>

int main() {
    using namespace matrice;
    matrix x = arange(data_type::kfloat, 4, 4, 1.0f);
    matrix y = arange(data_type::kfloat, 4, 4, 1.0f);
    auto z = matmul(x, y);
    std::cout << "shape = " << z.rows() << "x" << z.cols() << std::endl;
    std::cout << "size = " << z.size() << std::endl;
    std::cout << z;
    // Access indices with at
    std::cout << z.at<float>(0) << std::endl;
    // Modify a value
    z.mutable_data()[0] = 1.0f;
    // Access indices with the pointer to the data
    std::cout << z.data()[0] << std::endl;
}
```

## OCaml Bindins

In OCaml, two records `Vector.vector` and `Matrix.matrix` are defined.

A working example is the following:

```ocaml
open Matrice

(*Create an uninitialized matrix of float*)
let rows = 3;;
let cols = 4;;
let a = Matrix.make Float32 rows cols;;
Matrix.print_endline a;;

(*Create an uninitialized vector of int*)
let b = Vector.make Int32 10;;
(*Set the values*)
for i = 0 to 9 do
    Vector.set_int32 b (Int32.of_int i)
done
Vector.print_endline b;;

(*Matrix multiplication*)
let x = Matrix.arange_float 4 4 1.0;;
let y = Matrix.arange_float 4 4 1.0;;
Matrix.print_endline matmul x y;;
```

