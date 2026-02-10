module Types = Types_generated

module Functions = Function_description.Functions (Matrice__c_generated_functions__Function_description__Functions)

open Ctypes

(*Vector data type*)
type vector = {
   mutable own_data : bool;
   mutable ptr : Types.vector_ptr;
}

(*Get the length/size of a vector*)
let rec length(v : vector) : int = 
  Unsigned.UInt32.to_int (Functions.vector_ptr_size v.ptr)

(*Matrix data type*)
type matrix = {
   mutable own_data : bool;
   mutable ptr : Types.matrix_ptr;
}

