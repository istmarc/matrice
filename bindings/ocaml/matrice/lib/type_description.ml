open Ctypes

module Types (F: Ctypes.TYPE) = struct
   open F

   (*enum data_type*)
   let knone = constant "knone" uint32_t
   let kint = constant "kint" uint32_t
   let kfloat = constant "kfloat" uint32_t
   let kdouble = constant "kdouble" uint32_t

   (*float32 type*)
   let float32 : float typ = float

   (*struct vector*)
   type vector_ptr = unit ptr
   let vector_ptr : vector_ptr typ = ptr void

   (*struct matrix*)
   type matrix_ptr = unit ptr
   let matrix_ptr : matrix_ptr typ = ptr void

   (*void pointer*)
   type void_ptr = unit ptr
   let void_ptr : void_ptr typ = ptr void

   (*Array of uint32_t*)
   let uint32_ptr = ptr uint32_t

end

