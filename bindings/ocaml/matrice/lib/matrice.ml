module Types = Types_generated

module Functions = Function_description.Functions (Matrice__c_generated_functions__Function_description__Functions)

open Ctypes

(*data_type*)
type dtype = None | Int32 | Float32 | Float64

type value = ValueInt32 of int32 | ValueFloat32 of float | ValueFloat64 of float

let dtype_to_string(x : dtype) =
   match x with
      | None    -> "None"
      | Int32   -> "Int32"
      | Float32 -> "Float32"
      | Float64 -> "Float64"

let to_string(x : value) =
   match x with
   | ValueInt32 i -> Base.Int32.to_string i
   | ValueFloat32 f -> Base.Float.to_string f
   | ValueFloat64 d -> Base.Float.to_string d

let to_int = function
   | None -> Types.knone
   | Int32 -> Types.kint
   | Float32 -> Types.kfloat
   | Float64 -> Types.kdouble

let to_dtype = function
   | 0 -> None
   | 1 -> Int32
   | 2 -> Float32
   | 3 -> Float64
   | _ -> failwith "Unsupported type.\n"

module Vector = struct
   (*Vector data type*)
   type vector = {
      mutable ptr : Types.vector_ptr;
   }

   (*Create an initialized vector*)
   let make (t : dtype) (size : int) : vector= 
      {ptr = Functions.vector_ptr_make (to_int t) (Unsigned.UInt32.of_int size)}

   (*Get the length/size of a vector*)
   let length(v : vector) : int = 
     Unsigned.UInt32.to_int (Functions.vector_ptr_size v.ptr)

   (*Print a vector*)
   (*TODO Define a function to_string and rework print_endline*)
   let print_endline(v : vector) = 
      Functions.vector_ptr_print(v.ptr)

   (*Get the element at idx*)
   let get_int32 (v:vector) (idx) : int32 = 
      let xptr = Ctypes.allocate Ctypes.int32_t 0l in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.vector_ptr_at v.ptr idx voidxptr;
      Ctypes.(!@) xptr

   let get_float32 (v:vector) (idx) : float = 
      let xptr = Ctypes.allocate Ctypes.float 0. in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.vector_ptr_at v.ptr idx voidxptr;
      let ret = Ctypes.(!@) xptr in
      ret

   let get_float64 (v:vector) (idx) : float = 
      let xptr = Ctypes.allocate Ctypes.double 0. in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.vector_ptr_at v.ptr idx voidxptr;
      Ctypes.(!@) xptr

   (*Get the element at idx typed type*)
   let get (vec : vector) (t : dtype) (idx: int) : value  =
      let idx_uint32 = Unsigned.UInt32.of_int idx in
      match t with
         | None -> failwith "Unknown data type.\n"
         | Int32 -> ValueInt32(get_int32 vec idx_uint32)
         | Float32 -> ValueFloat32(get_float32 vec idx_uint32)
         | Float64 -> ValueFloat64(get_float64 vec idx_uint32)

   (*Set the element at idx to value*)
   let set_int32 (vec : vector) (idx : int) (v : int32) =
      let idx_uint32 = Unsigned.UInt32.of_int idx in
      let xptr = Ctypes.allocate Ctypes.int32_t v in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.vector_ptr_set vec.ptr idx_uint32 voidxptr

   let set_float32 (vec : vector) (idx : int) (v : float) =
      let idx_uint32 = Unsigned.UInt32.of_int idx in
      let xptr = Ctypes.allocate Ctypes.float v in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.vector_ptr_set vec.ptr idx_uint32 voidxptr

   let set_float64 (vec : vector) (idx : int) (v : float) =
      let idx_uint32 = Unsigned.UInt32.of_int idx in
      let xptr = Ctypes.allocate Ctypes.double v in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.vector_ptr_set vec.ptr idx_uint32 voidxptr

   let set_float (vec : vector) (idx : int) (v : float) =
      set_float64 vec idx v

   (*Create a vector of ones*)
   let ones (t : dtype) (size : int) : vector = 
      let idx_uint32 = Unsigned.UInt32.of_int size in
      let t_uint32 = to_int t in
      {ptr=Functions.vector_ptr_ones t_uint32 idx_uint32}

   (*Create a vector of zeros*)
   let zeros (t : dtype) (size : int) : vector = 
      let idx_uint32 = Unsigned.UInt32.of_int size in
      let t_uint32 = to_int t in
      {ptr=Functions.vector_ptr_zeros t_uint32 idx_uint32}

   (*Create a vector of arange*)
   let arange_int32 (size : int) (start : int) : vector =
      let n = Unsigned.UInt32.of_int size in
      let startint32 = Signed.Int32.of_int start in
      let xptr = Ctypes.allocate Ctypes.int32_t startint32 in
      let voidxptr = Ctypes.to_voidp xptr in
      let tuint32 = to_int Int32 in
      {ptr = Functions.vector_ptr_arange tuint32 n voidxptr}

   let arange_float32 (size : int) (start:float) : vector =
      let n = Unsigned.UInt32.of_int size in
      let xptr = Ctypes.allocate Ctypes.float start in
      let voidxptr = Ctypes.to_voidp xptr in
      let tuint32 = to_int Float32 in
      {ptr = Functions.vector_ptr_arange tuint32 n voidxptr}

   let arange_float64 (size : int) (start:float) : vector =
      let n = Unsigned.UInt32.of_int size in
      let xptr = Ctypes.allocate Ctypes.double start in
      let voidxptr = Ctypes.to_voidp xptr in
      let tuint32 = to_int Float64 in
      {ptr = Functions.vector_ptr_arange tuint32 n voidxptr}

   let arange_float (size : int) (start : float) : vector =
      arange_float64 size start

   (*Get the type*)
   let dtype(vec : vector) : dtype =
      let ty = Functions.vector_ptr_type vec.ptr in
      to_dtype (Unsigned.UInt32.to_int ty)

   (*Add two vectors*)
   let add (x : vector) (y : vector) : vector =
      let ty = dtype x in
      let size = length x in
      let z = make ty size in
      Functions.vector_ptr_add x.ptr y.ptr z.ptr;
      z

   (*Substract two vectors*)
   let sub (x : vector) (y : vector) : vector =
      let ty = dtype x in
      let size = length x in
      let z = make ty size in
      Functions.vector_ptr_sub x.ptr y.ptr z.ptr;
      z

   (*ELementwise multiply two vectors*)
   let mul (x : vector) (y : vector) : vector =
      let ty = dtype x in
      let size = length x in
      let z = make ty size in
      Functions.vector_ptr_mul x.ptr y.ptr z.ptr;
      z

   (*Elementwise divide two vectors*)
   let div (x : vector) (y : vector) : vector =
      let ty = dtype x in
      let size = length x in
      let z = make ty size in
      Functions.vector_ptr_div x.ptr y.ptr z.ptr;
      z

   (*Test if the values of two vectors are close to a precision eps*)
   let are_close (x:vector) (y:vector) (eps : float) : bool =
      Functions.vector_ptr_are_close x.ptr y.ptr eps

   (*Test if two vectors are equal, that is their pointer to the data are equal in memory*)
   let equal (x : vector) (y : vector) : bool =
      Functions.vector_ptr_equal x.ptr y.ptr
end

module Matrix = struct 

   (*Matrix data type*)
   type matrix = {
      mutable ptr : Types.matrix_ptr;
   }

   (*Create an initialized matrix*)
   let make (t : dtype) (rows : int) (cols : int) : matrix= 
      let tint = to_int t in
      let rowsint = Unsigned.UInt32.of_int rows in
      let colsint = Unsigned.UInt32.of_int cols in
      let shape = Ctypes.allocate_n Ctypes.uint32_t ~count:2 in
      shape <-@ rowsint;
      (shape +@ 1) <-@ colsint;
      {ptr = Functions.matrix_ptr_make tint shape}

   (*Get the number of rows of a matrix*)
   let rows(m : matrix) : int = 
     Unsigned.UInt32.to_int (Functions.matrix_ptr_rows m.ptr)

   (*Get the number of cols of a matrix*)
   let cols(m : matrix) : int = 
     Unsigned.UInt32.to_int (Functions.matrix_ptr_cols m.ptr)

   (*Print a matrix*)
   (*TODO Create a function to_string to convert a matrix to string and then print_endline*)
   let print_endline (m : matrix) = 
      Functions.matrix_ptr_print(m.ptr)

   (*Get the element at idx*)
   let get_int32 (v:matrix) (rows) (cols) : int32 = 
      let xptr = Ctypes.allocate Ctypes.int32_t 0l in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.matrix_ptr_get v.ptr rows cols voidxptr;
      Ctypes.(!@) xptr

   let get_float32 (v:matrix) (rows) (cols) : float = 
      let xptr = Ctypes.allocate Ctypes.float 0. in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.matrix_ptr_get v.ptr rows cols voidxptr;
      let ret = Ctypes.(!@) xptr in
      ret

   let get_float64 (v:matrix) (rows) (cols) : float = 
      let xptr = Ctypes.allocate Ctypes.double 0. in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.matrix_ptr_get v.ptr rows cols voidxptr;
      Ctypes.(!@) xptr

   (*Get the element at idx typed type*)
   let get (m : matrix) (t : dtype) (rows : int) (cols : int) : value  =
      let rows_uint32 = Unsigned.UInt32.of_int rows in
      let cols_uint32 = Unsigned.UInt32.of_int cols in
      match t with
         | None -> failwith "Unknown data type.\n"
         | Int32 -> ValueInt32(get_int32 m rows_uint32 cols_uint32)
         | Float32 -> ValueFloat32(get_float32 m rows_uint32 cols_uint32)
         | Float64 -> ValueFloat64(get_float64 m rows_uint32 cols_uint32)

   (*Get the element at the linear index*)
   let at_int32 (v:matrix) (index) : int32 = 
      let xptr = Ctypes.allocate Ctypes.int32_t 0l in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.matrix_ptr_at v.ptr index voidxptr;
      Ctypes.(!@) xptr

   let at_float32 (v:matrix) (index) : float = 
      let xptr = Ctypes.allocate Ctypes.float 0. in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.matrix_ptr_at v.ptr index voidxptr;
      let ret = Ctypes.(!@) xptr in
      ret

   let at_float64 (v:matrix) (index) : float = 
      let xptr = Ctypes.allocate Ctypes.double 0. in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.matrix_ptr_at v.ptr index voidxptr;
      Ctypes.(!@) xptr

   (*Get the element at idx typed type*)
   let at (m : matrix) (t : dtype) (index : int) : value  =
      let index_uint32 = Unsigned.UInt32.of_int index in
      match t with
         | None -> failwith "Unknown data type.\n"
         | Int32 -> ValueInt32(at_int32 m index_uint32)
         | Float32 -> ValueFloat32(at_float32 m index_uint32)
         | Float64 -> ValueFloat64(at_float64 m index_uint32)

   (*Set the element at row and col to value*)
   let set_int32 (m : matrix) (row : int) (col: int) (v : int32) =
      let row_uint32 = Unsigned.UInt32.of_int row in
      let col_uint32 = Unsigned.UInt32.of_int col in
      let xptr = Ctypes.allocate Ctypes.int32_t v in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.matrix_ptr_set m.ptr row_uint32 col_uint32 voidxptr

   let set_float32 (m : matrix) (row : int) (col: int) (v : float) =
      let row_uint32 = Unsigned.UInt32.of_int row in
      let col_uint32 = Unsigned.UInt32.of_int col in
      let xptr = Ctypes.allocate Ctypes.float v in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.matrix_ptr_set m.ptr row_uint32 col_uint32 voidxptr

   let set_float64 (m : matrix) (row : int) (col: int) (v : float) =
      let row_uint32 = Unsigned.UInt32.of_int row in
      let col_uint32 = Unsigned.UInt32.of_int col in
      let xptr = Ctypes.allocate Ctypes.double v in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.matrix_ptr_set m.ptr row_uint32 col_uint32 voidxptr

   let set_float (m : matrix) (row : int) (col : int) (v : float) =
      set_float64 m row col v

   (*Set the element at row and col to value*)
   let set_at_int32 (m : matrix) (idx : int) (v : int32) =
      let idx_uint32 = Unsigned.UInt32.of_int idx in
      let xptr = Ctypes.allocate Ctypes.int32_t v in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.matrix_ptr_set_at m.ptr idx_uint32 voidxptr

   let set_at_float32 (m : matrix) (idx : int) (v : float) =
      let idx_uint32 = Unsigned.UInt32.of_int idx in
      let xptr = Ctypes.allocate Ctypes.float v in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.matrix_ptr_set_at m.ptr idx_uint32 voidxptr

   let set_at_float64 (m : matrix) (idx : int) (v : float) =
      let idx_uint32 = Unsigned.UInt32.of_int idx in
      let xptr = Ctypes.allocate Ctypes.double v in
      let voidxptr = Ctypes.to_voidp xptr in
      Functions.matrix_ptr_set_at m.ptr idx_uint32 voidxptr

   let set_at_float (m : matrix) (idx : int) (v : float) =
      set_at_float64 m idx v

   (*Create a matrix of ones*)
   let ones (t : dtype) (rows : int) (cols : int): matrix = 
      let rows_uint32 = Unsigned.UInt32.of_int rows in
      let cols_uint32 = Unsigned.UInt32.of_int cols in
      let t_uint32 = to_int t in
      let shape = Ctypes.allocate_n Ctypes.uint32_t ~count:2 in
      shape <-@ rows_uint32;
      (shape +@ 1) <-@ cols_uint32;
      {ptr=Functions.matrix_ptr_ones t_uint32 shape}

   (*Create a matrix of zeros*)
   let zeros (t : dtype) (rows : int) (cols : int): matrix = 
      let rows_uint32 = Unsigned.UInt32.of_int rows in
      let cols_uint32 = Unsigned.UInt32.of_int cols in
      let t_uint32 = to_int t in
      let shape = Ctypes.allocate_n Ctypes.uint32_t ~count:2 in
      shape <-@ rows_uint32;
      (shape +@ 1) <-@ cols_uint32;
      {ptr=Functions.matrix_ptr_zeros t_uint32 shape}

   (*Create a matrix of arange*)
   let arange_int32 (rows : int) (cols : int) (start : int) : matrix =
      let rows_uint32 = Unsigned.UInt32.of_int rows in
      let cols_uint32 = Unsigned.UInt32.of_int cols in
      let startint32 = Signed.Int32.of_int start in
      let xptr = Ctypes.allocate Ctypes.int32_t startint32 in
      let voidxptr = Ctypes.to_voidp xptr in
      let tuint32 = to_int Int32 in
      let shape = Ctypes.allocate_n Ctypes.uint32_t ~count:2 in
      shape <-@ rows_uint32;
      (shape +@ 1) <-@ cols_uint32;
      {ptr = Functions.matrix_ptr_arange tuint32 shape voidxptr}

   (*Create a matrix of arange*)
   let arange_float32 (rows : int) (cols : int) (start : float) : matrix =
      let rows_uint32 = Unsigned.UInt32.of_int rows in
      let cols_uint32 = Unsigned.UInt32.of_int cols in
      let xptr = Ctypes.allocate Ctypes.float start in
      let voidxptr = Ctypes.to_voidp xptr in
      let tuint32 = to_int Float32 in
      let shape = Ctypes.allocate_n Ctypes.uint32_t ~count:2 in
      shape <-@ rows_uint32;
      (shape +@ 1) <-@ cols_uint32;
      {ptr = Functions.matrix_ptr_arange tuint32 shape voidxptr}

   (*Create a matrix of arange*)
   let arange_float64 (rows : int) (cols : int) (start : float) : matrix =
      let rows_uint32 = Unsigned.UInt32.of_int rows in
      let cols_uint32 = Unsigned.UInt32.of_int cols in
      let xptr = Ctypes.allocate Ctypes.double start in
      let voidxptr = Ctypes.to_voidp xptr in
      let tuint32 = to_int Float64 in
      let shape = Ctypes.allocate_n Ctypes.uint32_t ~count:2 in
      shape <-@ rows_uint32;
      (shape +@ 1) <-@ cols_uint32;
      {ptr = Functions.matrix_ptr_arange tuint32 shape voidxptr}

   let arange_float (rows : int) (cols : int) (start : float) : matrix =
      arange_float64 rows cols start

   (*Get the type*)
   let dtype(m : matrix) : dtype =
      let ty = Functions.matrix_ptr_type m.ptr in
      to_dtype (Unsigned.UInt32.to_int ty)

   (*Add two matrices*)
   let add (x : matrix) (y : matrix) : matrix =
      let ty = dtype x in
      let nrows = rows x in
      let ncols = cols x in
      let z = make ty nrows ncols in
      Functions.matrix_ptr_add x.ptr y.ptr z.ptr;
      z

   (*Substract two matrices*)
   let sub (x : matrix) (y : matrix) : matrix =
      let ty = dtype x in
      let nrows = rows x in
      let ncols = cols x in
      let z = make ty nrows ncols in
      Functions.matrix_ptr_sub x.ptr y.ptr z.ptr;
      z

   (*Elementwise multiply two matrices*)
   let mul (x : matrix) (y : matrix) : matrix =
      let ty = dtype x in
      let nrows = rows x in
      let ncols = cols x in
      let z = make ty nrows ncols in
      Functions.matrix_ptr_mul x.ptr y.ptr z.ptr;
      z

   (*Divide two matrices*)
   let div (x : matrix) (y : matrix) : matrix =
      let ty = dtype x in
      let nrows = rows x in
      let ncols = cols x in
      let z = make ty nrows ncols in
      Functions.matrix_ptr_div x.ptr y.ptr z.ptr;
      z

   (*Matrix multiplication*)
   let matmul (x : matrix) (y : matrix) : matrix = 
      let ty = dtype x in
      let nrows = rows x in
      let ncols = cols y in
      let z = make ty nrows ncols in
      Functions.matrix_ptr_matmul x.ptr y.ptr z.ptr;
      z

   (*Test if the values of two matrixs are close to a precision eps*)
   let are_close (x:matrix) (y:matrix) (eps : float) : bool =
      Functions.matrix_ptr_are_close x.ptr y.ptr eps

   (*Test if two matrixs are equal, that is their pointer to the data are equal in memory*)
   let equal (x : matrix) (y : matrix) : bool =
      Functions.matrix_ptr_equal x.ptr y.ptr
end

