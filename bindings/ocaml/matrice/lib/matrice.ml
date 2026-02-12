module Types = Types_generated

module Functions = Function_description.Functions (Matrice__c_generated_functions__Function_description__Functions)

open Ctypes

(*data_type*)
type dtype = Knone | Kint | Kfloat | Kdouble

let to_int = function
   | Knone -> Types.knone
   | Kint -> Types.kint
   | Kfloat -> Types.kfloat
   | Kdouble -> Types.kdouble

(*Vector data type*)
type vector = {
   mutable own_data : bool;
   mutable ptr : Types.vector_ptr;
}

(*Create an initialized vector*)
let make_vector (t : dtype) (size : int) : vector= 
   {own_data = true; ptr = Functions.vector_ptr_make (to_int t) (Unsigned.UInt32.of_int size)}


(*Get the length/size of a vector*)
let length(v : vector) : int = 
  Unsigned.UInt32.to_int (Functions.vector_ptr_size v.ptr)

(*Print a vector*)
let print_vector(v : vector) = 
   Functions.vector_ptr_print(v.ptr)

(*Matrix data type*)
type matrix = {
   mutable own_data : bool;
   mutable ptr : Types.matrix_ptr;
}

