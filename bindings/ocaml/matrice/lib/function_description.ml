open Ctypes

module Types = Types_generated

module Functions (F: Ctypes.FOREIGN) = struct
open F

   (*Vector bindings*)
   let vector_ptr_make = foreign "vector_make" (uint32_t @-> uint32_t @-> returning Types.vector_ptr)

   let vector_ptr_free = foreign "vector_free" (Types.vector_ptr @-> returning void)

   let vector_ptr_print = foreign "vector_print" (Types.vector_ptr @-> returning void)

   let vector_ptr_size = foreign "vector_size" (Types.vector_ptr @-> returning uint32_t)

   let vector_ptr_at = foreign "vector_at" (Types.vector_ptr @-> uint32_t @-> Types.void_ptr @-> returning void)

   let vector_ptr_set = foreign "vector_set" (Types.vector_ptr @-> uint32_t @-> Types.void_ptr @-> returning void)

   let vector_ptr_add = foreign "vector_add" (Types.vector_ptr @-> Types.vector_ptr @-> Types.vector_ptr @-> returning void)

   let vector_ptr_sub = foreign "vector_sub" (Types.vector_ptr @-> Types.vector_ptr @-> Types.vector_ptr @-> returning void)

   let vector_ptr_mul = foreign "vector_mul" (Types.vector_ptr @-> Types.vector_ptr @-> Types.vector_ptr @-> returning void)

   let vector_ptr_div = foreign "vector_div" (Types.vector_ptr @-> Types.vector_ptr @-> Types.vector_ptr @-> returning void)

   let vector_ptr_ones = foreign "vector_ones" (uint32_t @-> uint32_t @-> returning Types.vector_ptr)

   let vector_ptr_zeros = foreign "vector_zeros" (uint32_t @-> uint32_t @-> returning Types.vector_ptr)

   let vector_ptr_fill = foreign "vector_fill" (uint32_t @-> uint32_t @-> Types.void_ptr @-> returning Types.vector_ptr)

   let vector_ptr_arange = foreign "vector_arange" (uint32_t @-> uint32_t @-> Types.void_ptr @-> returning Types.vector_ptr)

   let vector_ptr_are_close = foreign "vector_are_close" (Types.vector_ptr @-> Types.vector_ptr @-> double @-> returning bool)

   let vector_ptr_equal = foreign "vector_equal" (Types.vector_ptr @-> Types.vector_ptr @-> returning bool)

   (*Matrix bindings*)
   let matrix_ptr_make = foreign "matrix_make" (uint32_t @-> Types.uint32_ptr @-> returning Types.matrix_ptr)

   let matrix_ptr_free = foreign "matrix_free" (Types.matrix_ptr @-> returning void)

   let matrix_ptr_print = foreign "matrix_print" (Types.matrix_ptr @-> returning void)

   let matrix_ptr_rows = foreign "matrix_rows" (Types.matrix_ptr @-> returning uint32_t)

   let matrix_ptr_cols = foreign "matrix_cols" (Types.matrix_ptr @-> returning uint32_t)

   let matrix_ptr_get = foreign "matrix_get" (Types.matrix_ptr @-> uint32_t @-> uint32_t @-> Types.void_ptr @-> returning void)

   let matrix_ptr_at = foreign "matrix_at" (Types.matrix_ptr @-> uint32_t @-> Types.void_ptr @-> returning void)

   let matrix_ptr_set = foreign "matrix_set" (Types.matrix_ptr @-> uint32_t @-> uint32_t @-> Types.void_ptr @-> returning void)

   let matrix_ptr_set = foreign "matrix_set_at" (Types.matrix_ptr @-> uint32_t @-> Types.void_ptr @-> returning void)

   let matrix_ptr_add = foreign "matrix_add" (Types.matrix_ptr @-> Types.matrix_ptr @-> Types.matrix_ptr @-> returning void)

   let matrix_ptr_sub = foreign "matrix_sub" (Types.matrix_ptr @-> Types.matrix_ptr @-> Types.matrix_ptr @-> returning void)

   let matrix_ptr_mul = foreign "matrix_mul" (Types.matrix_ptr @-> Types.matrix_ptr @-> Types.matrix_ptr @-> returning void)

   let matrix_ptr_div = foreign "matrix_div" (Types.matrix_ptr @-> Types.matrix_ptr @-> Types.matrix_ptr @-> returning void)

   let matrix_ptr_matmul = foreign "matrix_matmul" (Types.matrix_ptr @-> Types.matrix_ptr @-> Types.matrix_ptr @-> returning void)

   let matrix_ptr_ones = foreign "matrix_ones" (uint32_t @-> Types.uint32_ptr @-> returning Types.matrix_ptr)

   let matrix_ptr_zeros = foreign "matrix_zeros" (uint32_t @-> Types.uint32_ptr @-> returning Types.matrix_ptr)

   let matrix_ptr_fill = foreign "matrix_fill" (uint32_t @-> Types.uint32_ptr @-> Types.void_ptr @-> returning Types.matrix_ptr)

   let matrix_ptr_arange = foreign "matrix_arange" (uint32_t @-> Types.uint32_ptr @-> Types.void_ptr @-> returning Types.matrix_ptr)

   let matrix_ptr_are_close = foreign "matrix_are_close" (Types.matrix_ptr @-> Types.matrix_ptr @-> double @-> returning bool)

   let matrix_ptr_equal = foreign "matrix_equal" (Types.matrix_ptr @-> Types.matrix_ptr @-> returning bool)

end
