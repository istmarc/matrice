open Ctypes

module Types = Types_generated

module Functions (F: Ctypes.FOREIGN) = struct
open F

   (*Vector bindings*)
   let vector_ptr_make = foreign "vector_make" (uint32_t @-> uint32_t @-> returning Types.vector_ptr)

   let vector_ptr_free = foreign "vector_free" (Types.vector_ptr @-> returning void)

   let vector_ptr_print = foreign "vector_print" (Types.vector_ptr @-> returning void)

   let vector_ptr_size = foreign "vector_size" (Types.vector_ptr @-> returning uint32_t)


end
