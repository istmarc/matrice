open Matrice

(*Test for vector of int*)
let iv = Vector.make Int32 10;;
let niv = Vector.length iv;;

print_endline "Values of the vector=";;

for i = 0 to niv-1 do
   Vector.set_int32 iv i (Int32.of_int i)
done;;

for i = 0 to niv-1 do
   let ivi = Vector.get iv Int32 i in
   print_endline (to_string ivi)
done;;

Printf.printf "Length = %s\n" (Int.to_string niv);;
print_endline "Print the vector";;
Vector.print_endline iv;;

(*Test for vector of float*)
let fv = Vector.make Float32 10;;
let nfv = Vector.length fv;;

print_endline "Values of the vector=";;

for i = 0 to nfv-1 do
   Vector.set_float32 fv i (Float.of_int i)
done;;

for i = 0 to nfv-1 do
   let fvi = Vector.get fv Float32 i in
   print_endline (to_string fvi)
done;;

Printf.printf "Length = %s\n" (Int.to_string nfv);;
print_endline "Print the vector";;
Vector.print_endline fv;;

(*Test for vector of double*)
let dv = Vector.make Float64 10;;
let ndv = Vector.length dv;;

print_endline "Values of the vector=";;

for i = 0 to ndv-1 do
   Vector.set_float64 dv i (Float.of_int i)
done;;

for i = 0 to ndv-1 do
   let dvi = Vector.get dv Float64 i in
   print_endline (to_string dvi)
done;;

Printf.printf "Length = %s\n" (Int.to_string ndv);;
print_endline "Print the vector";;
Vector.print_endline dv;;

(*Test vector of ones*)
let fones = Vector.ones Float32 10;;
Vector.print_endline fones;;

(*Test vector of range*)
let iarange = Vector.arange_int32 10 1;;
Vector.print_endline iarange;;

let farange = Vector.arange_float32 10 1.;;
Vector.print_endline farange;;

let darange = Vector.arange_float 10 1.;;
Vector.print_endline darange;;

(*Test data type*)
print_endline (dtype_to_string (Vector.dtype iarange));
print_endline (dtype_to_string (Vector.dtype farange));
print_endline (dtype_to_string (Vector.dtype darange));

(*Test add two vectors*)
Vector.print_endline (Vector.add farange farange);;
Vector.print_endline (Vector.sub farange farange);;
Vector.print_endline (Vector.mul farange farange);;
Vector.print_endline (Vector.div farange farange);;

(*Test if they are close*)
print_endline (Bool.to_string (Vector.are_close farange farange 1e-3));;
print_endline (Bool.to_string (Vector.equal farange farange));;


(*Test a matrix*)
let im = Matrix.make Int32 2 3;;
Matrix.print_endline im;;

let fm = Matrix.make Float32 2 3;;
Matrix.print_endline fm;;

let dm = Matrix.make Float64 2 3;;
Matrix.print_endline dm;;

(*Set and get values*)
let rows = 2;;
let cols = 3;;

for i = 0 to rows-1 do
   for j = 0 to cols-1 do
      let v = Float.of_int (i+j) in
      Matrix.set_float32 fm i j v
   done
done;;

print_endline "Value after set";;
for i = 0 to rows-1 do
   for j = 0 to cols-1 do
      let v = Matrix.get fm Float32 i j in
      Printf.printf "%s\n" (to_string v)
   done
done;;

Matrix.print_endline fm;;


(*Test matrix of ones*)
let imones = Matrix.ones Int32 3 4;;
Matrix.print_endline imones;;

let fmones = Matrix.ones Float32 3 4;;
Matrix.print_endline fmones;;

let dmones = Matrix.ones Float64 3 4;;
Matrix.print_endline dmones;;

(*Matrix of arange*)
let imarange = Matrix.arange_int32 2 3 1;;
Matrix.print_endline imarange;;

let fmarange = Matrix.arange_float32 2 3 1.0;;
Matrix.print_endline fmarange;;

let dmarange = Matrix.arange_float 2 3 1.0;;
Matrix.print_endline dmarange;;

(*Matrix operations*)
Matrix.print_endline (Matrix.add fmarange fmarange);;
Matrix.print_endline (Matrix.sub fmarange fmarange);;
Matrix.print_endline (Matrix.mul fmarange fmarange);;
Matrix.print_endline (Matrix.div fmarange fmarange);;

(*Matrix multiplication*)
let a = Matrix.arange_float 4 4 1.0;;
let b = Matrix.arange_float 4 4 1.0;;
Matrix.print_endline (Matrix.matmul a b);;

