#ifndef MATRICE_HXX
#define MATRICE_HXX

#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <string>

namespace matrice_c {
#include <matrice/data_type.h>
#include <matrice/matrix.h>
#include <matrice/vector.h>
} // namespace matrice_c

namespace matrice {

using data_type = matrice_c::data_type;

template<class T>
auto to_data_type() -> decltype(auto) {
   /*
   static_assert(std::is_same_v<T, int> || std::is_same_v<T, int64_t>,
      || std::is_same_v<T, float> || std::is_same_v<T, double>, "Unsupported type.\n");
   */

   if (std::is_same_v<T, int>)
      return data_type::kint32;
   else if (std::is_same_v<T, float>)
      return data_type::kfloat32;
   else if (std::is_same_v<T, double>)
      return data_type::kfloat64;
   else if (std::is_same_v<T, int64_t>)
      return data_type::kint64;

}

// Bindings for vector data type
template<class T>
class vector {
public:
   /*static_assert(std::is_same_v<T, int> || std::is_same_v<T, int64_t>,
      || std::is_same_v<T, float> || std::is_same_v<T, double>, "Unsupported data type");
   */

 private:
   bool own_vec = true;
   matrice_c::vector *vec = nullptr;

 public:
   // Construct a vector from type and size
   vector(uint32_t size) : own_vec(true) {
      data_type type = to_data_type<T>();
      vec = matrice_c::vector_make(type, size);
   }

   vector(matrice_c::vector *ptr) : own_vec(true), vec(ptr) {}

   // Constructor from an initializer_list of T
   vector(std::initializer_list<T> values) : own_vec(true) {
      uint32_t size = values.size();
      data_type type = to_data_type<T>();
      vec = matrice_c::vector_make(type, size);
      T *data = (T *)vec->data;
      for (uint32_t i = 0; i < size; i++) {
         data[i] = *(values.begin() + i);
      }
   }

   ~vector() {
      if (own_vec && vec)
         matrice_c::vector_free(vec);
   }

   // Copy constructor
   vector(const vector &v) {
      if (own_vec && vec) {
         matrice_c::vector_free(vec);
      }
      own_vec = false;
      vec = v.vec;
   }

   // Move constructor
   vector(vector &&v) {
      if (own_vec && vec) {
         matrice_c::vector_free(vec);
      }
      own_vec = std::move(v.own_vec);
      vec = std::move(v.vec);
   }

   // = oprator
   vector &operator=(const vector &v) {
      if (own_vec && vec) {
         matrice_c::vector_free(vec);
      }
      own_vec = false;
      vec = v.vec;
      return *this;
   }

   // = operator
   vector &operator=(vector &&v) {
      if (own_vec && vec) {
         matrice_c::vector_free(vec);
      }
      own_vec = std::move(v.own_vec);
      vec = std::move(v.vec);
      return *this;
   }

   // Get the type
   data_type type() const { return vec->type; }

   // Get the size
   uint32_t size() const { return vec->size; }

   // Get the ptr to vector
   matrice_c::vector *ptr() const { return vec; }

   // Get the pointer to the data
   const T *data() const {
      return (const T *)vec->data;
   }

   // Get the mutable pointer to the data
   T *mutable_data() {
      return (T *)vec->data;
   }

   // Get the value at index
   const T &at(uint32_t index) const {
      T *data = (T *)vec->data;
      return data[index];
   }

   // Get the value at index
   T &at(uint32_t index) {
      T *data = (T *)vec->data;
      return data[index];
   }

   // Overload the [] operator
   const T &operator[](uint32_t index) const {
      return at(index);
   }
   T &operator[](uint32_t index) {
      return at(index);
   }

   // Overload the () operator
   const T &operator()(uint32_t index) const {
      return at(index);
   }
   T &operator()(uint32_t index) {
      return at(index);
   }

   vector operator+(const vector &v) {
      vector out(vec->size);
      matrice_c::vector_add(vec, v.ptr(), out.ptr());
      return out;
   }

   vector operator-(const vector &v) {
      vector out(vec->size);
      matrice_c::vector_sub(vec, v.ptr(), out.ptr());
      return out;
   }

   vector operator*(const vector &v) {
      vector out(vec->size);
      matrice_c::vector_mul(vec, v.ptr(), out.ptr());
      return out;
   }

   vector operator/(const vector &v) {
      vector out(vec->size);
      matrice_c::vector_div(vec, v.ptr(), out.ptr());
      return out;
   }

   template<class Ty>
   friend std::ostream &operator<<(std::ostream &os, const vector<Ty> &);
};

template<class T>
std::ostream &operator<<(std::ostream &os, const vector<T> &v) {
   matrice_c::vector_print(v.vec);
   return os;
}

// Instantiation
template class vector<int>;
template class vector<float>;
template class vector<double>;
template class vector<int64_t>;

// Create a vector of ones
template<class T>
vector<T> ones(uint32_t size) {
   data_type type = to_data_type<T>();
   matrice_c::vector *ptr = matrice_c::vector_ones(type, size);
   return vector<T>(ptr);
}

// Create a vector of zeros
template<class T>
vector<T> zeros(uint32_t size) {
   data_type type = to_data_type<T>();
   matrice_c::vector *ptr = matrice_c::vector_zeros(type, size);
   return vector<T>(ptr);
}

// Create a vector of arange
template<class T>
vector<T> arange(uint32_t size, T value = T(0)) {
   data_type type = to_data_type<T>();
   matrice_c::vector *ptr =
       matrice_c::vector_arange(type, size, &value);
   return vector<T>(ptr);
}

// Binding for matrix type
template<class T>
class matrix {
 private:
   bool own_mat = true;
   matrice_c::matrix *mat = nullptr;

 public:
   // Create a matrix from type, rows and columns
   matrix(uint32_t rows, uint32_t cols) : own_mat(true) {
      data_type type = to_data_type<T>();
      uint32_t shape[2] = {rows, cols};
      mat = matrice_c::matrix_make(type, shape);
   }

   // Create a matrix from type and shape
   matrix(uint32_t shape[2]) : own_mat(true) {
      data_type type = to_data_type<T>();
      mat = matrice_c::matrix_make(type, shape);
   }

   matrix(matrice_c::matrix *ptr) : own_mat(true) {
      if (own_mat && mat)
         matrice_c::matrix_free(mat);
      mat = ptr;
   }

   // Constructor from an initializer_list<initializer_list> of T
   // The list are interpreted as the columns of the matrix
   matrix(std::initializer_list<std::initializer_list<T>> values)
       : own_mat(true) {
      uint32_t cols = values.size();
      uint32_t rows = values.begin()->size();
      for (uint32_t i = 1; i < cols; i++) {
         auto iter = values.begin() + i;
         if (iter->size() != rows) {
            throw std::runtime_error("Invalid matrix row size.\n");
         }
      }
      uint32_t shape[2] = {rows, cols};
      data_type type = to_data_type<T>();
      mat = matrice_c::matrix_make(type, shape);
      int *data = (int *)mat->data;
      for (uint32_t j = 0; j < cols; j++) {
         for (uint32_t i = 0; i < rows; i++) {
            auto iter_col = values.begin() + j;
            auto iter_row = (*iter_col).begin() + i;
            data[i * mat->strides[0] + j * mat->strides[1]] = *iter_row;
         }
      }
   }

   // Get the pointer to the matrix
   matrice_c::matrix *ptr() const { return mat; }

   // Get the data type
   data_type type() const { return mat->type; }

   // Get the size
   uint32_t size() const { return mat->size; }

   // Get the shape
   uint32_t *shape() const { return mat->shape; }

   // Get the number of rows
   uint32_t rows() const { return mat->shape[0]; }

   // Get the number of columns
   uint32_t cols() const { return mat->shape[1]; }

   ~matrix() {
      if (own_mat && mat) {
         matrice_c::matrix_free(mat);
      }
   }

   // Copy constructor
   matrix(const matrix &m) {
      if (own_mat && mat) {
         matrice_c::matrix_free(mat);
      }
      own_mat = false;
      mat = m.mat;
   }

   // Move constructor
   matrix(matrix &&m) {
      if (own_mat && mat) {
         matrice_c::matrix_free(mat);
      }
      own_mat = std::move(m.own_mat);
      mat = std::move(m.mat);
   }

   // = operator
   matrix &operator=(const matrix &m) {
      if (own_mat && mat) {
         matrice_c::matrix_free(mat);
      }
      own_mat = false;
      mat = m.mat;
      return *this;
   }

   // = operator
   matrix &operator=(matrix &&m) {
      if (own_mat && mat) {
         matrice_c::matrix_free(mat);
      }
      own_mat = std::move(m.own_mat);
      mat = std::move(m.mat);
      return *this;
   }

   // Get the pointer to the data casted to a type T
   const T *data() const {
      return (const T *)mat->data;
   }

   // Get the pointer to the mutable data casted to a type T
   T *mutable_data() {
      return (T *)mat->data;
   }

   // Get the value at index casted to a type T
   const T &at(uint32_t index) const {
      T *data = (T *)mat->data;
      return data[index];
   }

   // Get the value at index casted to a type T
   T &at(uint32_t index) {
      T *data = (T *)mat->data;
      return data[index];
   }

public:
   const T&operator()(uint32_t row, uint32_t col) const {
      uint32_t offset = matrice_c::matrix_compute_offset(mat->strides, row, col);
      return at(offset);
   }

   T&operator()(uint32_t row, uint32_t col) {
      uint32_t offset = matrice_c::matrix_compute_offset(mat->strides, row, col);
      return at(offset);
   }

   const T& operator[](uint32_t index) const { return at(index);}
   T& operator[](uint32_t index) {return at(index);}

   // Elementwise matrix addition
   matrix operator+(const matrix &m) {
      matrix out(mat->shape);
      matrice_c::matrix_add(mat, m.ptr(), out.ptr());
      return out;
   }

   // Elementwise matrix substraction
   matrix operator-(const matrix &m) {
      matrix out(mat->shape);
      matrice_c::matrix_sub(mat, m.ptr(), out.ptr());
      return out;
   }

   // Elementwise matrix multiplication
   matrix operator*(const matrix &m) {
      matrix out(mat->shape);
      matrice_c::matrix_mul(mat, m.ptr(), out.ptr());
      return out;
   }

   // Elementwise matrix division
   matrix operator/(const matrix &m) {
      matrix out(mat->shape);
      matrice_c::matrix_div(mat, m.ptr(), out.ptr());
      return out;
   }

   template<class Ty>
   friend std::ostream &operator<<(std::ostream &os, const matrix<Ty> &);
};

template<class T>
std::ostream &operator<<(std::ostream &os, const matrix<T> &mat) {
   matrice_c::matrix_print(mat.ptr());
   return os;
}

// Instantiation
template class matrix<float>;
template class matrix<int>;
template class matrix<double>;
template class matrix<int64_t>;

// Matrix multiplication
template<class T>
matrix<T> matmul(const matrix<T> &x, const matrix<T> &y) {
   if (x.cols() != y.rows()) {
      throw std::runtime_error(
          "Matrix multiplication error, incompatible matrice shapes.\n");
   }
   uint32_t rows = x.rows();
   uint32_t cols = y.cols();
   matrix<T> z(rows, cols);
   matrice_c::matrix_matmul(x.ptr(), y.ptr(), z.ptr());
   return z;
}

// Create a matrix of ones
template<class T>
matrix<T> ones(uint32_t rows, uint32_t cols) {
   data_type type = to_data_type<T>();
   uint32_t shape[2] = {rows, cols};
   matrice_c::matrix *ptr = matrice_c::matrix_ones(type, shape);
   return matrix<T>(ptr);
}

// Create a matrix of zeros
template<class T>
matrix<T> zeros(uint32_t rows, uint32_t cols) {
   data_type type = to_data_type<T>();
   uint32_t shape[2] = {rows, cols};
   matrice_c::matrix *ptr = matrice_c::matrix_zeros(type, shape);
   return matrix<T>(ptr);
}

// Create a matrix of arange
template<class T>
matrix<T> arange(uint32_t rows, uint32_t cols, T value = T(0)) {
   uint32_t shape[2] = {rows, cols};
   data_type type = to_data_type<T>();
   matrice_c::matrix *ptr =
       matrice_c::matrix_arange(type, shape, &value);
   return matrix<T>(ptr);
}

} // namespace matrice

#endif
