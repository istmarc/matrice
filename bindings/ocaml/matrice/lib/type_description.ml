open Ctypes

module Types (F: Ctypes.TYPE) = struct
   open F

   (*Enum data_type*)
   let kint = constant "kint" uint32_t
   let kfloat = constant "kfloat" uint32_t
   let kdouble = constant "kdouble" uint32_t
   let knone = constant "knone" uint32_t

   (*struct vector*)
   type vector_ptr = unit ptr
   let vector_ptr : vector_ptr typ = ptr void

   (*struct matrix*)
   type matrix_ptr = unit ptr
   let matrix_ptr : matrix_ptr typ = ptr void

end

