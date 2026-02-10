#ifndef MATRICE_HXX
#define MATRICE_HXX

#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <stdexcept>

namespace matrice_c {
#include <matrice/data_type.h>
#include <matrice/matrix.h>
#include <matrice/vector.h>
} // namespace matrice_c

namespace matrice {

using data_type = matrice_c::data_type;

// Bindings for vector data type
class vector {
 private:
   bool own_vec = true;
   matrice_c::vector *vec = nullptr;

 public:
   // Construct a vector from type and size
   vector(data_type type, uint32_t size) : own_vec(true) {
      vec = matrice_c::vector_make(type, size);
   }

   vector(matrice_c::vector *ptr) : own_vec(true), vec(ptr) {}

   // Constructor from an initializer_list of int
   vector(std::initializer_list<int> values) : own_vec(true) {
      uint32_t size = values.size();
      vec = matrice_c::vector_make(data_type::kint, size);
      int *data = (int *)vec->data;
      for (uint32_t i = 0; i < size; i++) {
         data[i] = *(values.begin() + i);
      }
   }

   // Constructor from an initializer_list of float
   vector(std::initializer_list<float> values) : own_vec(true) {
      uint32_t size = values.size();
      vec = matrice_c::vector_make(data_type::kfloat, size);
      float *data = (float *)vec->data;
      for (uint32_t i = 0; i < size; i++) {
         data[i] = *(values.begin() + i);
      }
   }

   // Constructor from an initializer_list of double
   vector(std::initializer_list<double> values) : own_vec(true) {
      uint32_t size = values.size();
      vec = matrice_c::vector_make(data_type::kdouble, size);
      double *data = (double *)vec->data;
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
   template <class T> const T *data() const {
      static_assert(std::is_same_v<T, int> || std::is_same_v<T, float> ||
                        std::is_same_v<T, double>,
                    "Incompatible output type.");
      return (const T *)vec->data;
   }

   // Get the mutable pointer to the data
   template <class T> T *mutable_data() {
      static_assert(std::is_same_v<T, int> || std::is_same_v<T, float> ||
                        std::is_same_v<T, double>,
                    "Unsupported type.\n");
      return (T *)vec->data;
   }

   // Get the value at index
   template <class T> const T &at(uint32_t index) const {
      static_assert(std::is_same_v<T, int> || std::is_same_v<T, float> ||
                        std::is_same_v<T, double>,
                    "Unsupported type.\n");
      T *data = (T *)vec->data;
      return data[index];
   }

   // Get the value at index
   template <class T> T &at(uint32_t index) {
      static_assert(std::is_same_v<T, int> || std::is_same_v<T, float> ||
                        std::is_same_v<T, double>,
                    "Unsupported type.\n");
      T *data = (T *)vec->data;
      return data[index];
   }

   vector operator+(const vector &v) {
      vector out(vec->type, vec->size);
      matrice_c::vector_add(vec, v.ptr(), out.ptr());
      return out;
   }

   vector operator-(const vector &v) {
      vector out(vec->type, vec->size);
      matrice_c::vector_sub(vec, v.ptr(), out.ptr());
      return out;
   }

   vector operator*(const vector &v) {
      vector out(vec->type, vec->size);
      matrice_c::vector_mul(vec, v.ptr(), out.ptr());
      return out;
   }

   vector operator/(const vector &v) {
      vector out(vec->type, vec->size);
      matrice_c::vector_div(vec, v.ptr(), out.ptr());
      return out;
   }

   friend std::ostream &operator<<(std::ostream &os, const vector &);
};

std::ostream &operator<<(std::ostream &os, const vector &v) {
   matrice_c::vector_print(v.vec);
   return os;
}

// Create a vector of ones
vector ones(data_type type, uint32_t size) {
   matrice_c::vector *ptr = matrice_c::vector_ones(type, size);
   return vector(ptr);
}

// Create a vector of zeros
vector zeros(data_type type, uint32_t size) {
   matrice_c::vector *ptr = matrice_c::vector_zeros(type, size);
   return vector(ptr);
}

// Create a vector of arange
vector arange(uint32_t size, int value = 0) {
   matrice_c::vector *ptr =
       matrice_c::vector_arange(data_type::kint, size, &value);
   return vector(ptr);
}
vector arange(uint32_t size, float value = 0.f) {
   matrice_c::vector *ptr =
       matrice_c::vector_arange(data_type::kfloat, size, &value);
   return vector(ptr);
}
vector arange(uint32_t size, double value = 0.) {
   matrice_c::vector *ptr =
       matrice_c::vector_arange(data_type::kdouble, size, &value);
   return vector(ptr);
}

// Binding for matrix type
class matrix {
 private:
   bool own_mat = true;
   matrice_c::matrix *mat = nullptr;

 public:
   // Create a matrix from type, rows and columns
   matrix(data_type type, uint32_t rows, uint32_t cols) : own_mat(true) {
      uint32_t shape[2] = {rows, cols};
      mat = matrice_c::matrix_make(type, shape);
   }

   // Create a matrix from type and shape
   matrix(data_type type, uint32_t shape[2]) : own_mat(true) {
      mat = matrice_c::matrix_make(type, shape);
   }

   matrix(matrice_c::matrix *ptr) : own_mat(true) {
      if (own_mat && mat)
         matrice_c::matrix_free(mat);
      mat = ptr;
   }

   // Constructor from an initializer_list<initializer_list> of int
   // The list are interpreted as the columns of the matrix
   matrix(std::initializer_list<std::initializer_list<int>> values)
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
      mat = matrice_c::matrix_make(data_type::kint, shape);
      int *data = (int *)mat->data;
      for (uint32_t j = 0; j < cols; j++) {
         for (uint32_t i = 0; i < rows; i++) {
            auto iter_col = values.begin() + j;
            auto iter_row = (*iter_col).begin() + i;
            data[i * mat->strides[0] + j * mat->strides[1]] = *iter_row;
         }
      }
   }

   // Constructor from an initializer_list<initializer_list> of float
   // The list are interpreted as the columns of the matrix
   matrix(std::initializer_list<std::initializer_list<float>> values)
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
      mat = matrice_c::matrix_make(data_type::kfloat, shape);
      float *data = (float *)mat->data;
      for (uint32_t j = 0; j < cols; j++) {
         for (uint32_t i = 0; i < rows; i++) {
            auto iter_col = values.begin() + j;
            auto iter_row = (*iter_col).begin() + i;
            data[i * mat->strides[0] + j * mat->strides[1]] = *iter_row;
         }
      }
   }

   // Constructor from an initializer_list<initializer_list> of double
   // The list are interpreted as the columns of the matrix
   matrix(std::initializer_list<std::initializer_list<double>> values)
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
      mat = matrice_c::matrix_make(data_type::kdouble, shape);
      double *data = (double *)mat->data;
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
   template <class T> const T *data() const {
      static_assert(std::is_same_v<T, int> || std::is_same_v<T, float> ||
                        std::is_same_v<T, double>,
                    "Incompatible output type.");
      return (const T *)mat->data;
   }

   // Get the pointer to the mutable data casted to a type T
   template <class T> T *mutable_data() {
      static_assert(std::is_same_v<T, int> || std::is_same_v<T, float> ||
                        std::is_same_v<T, double>,
                    "Incompatible output type.");
      return (T *)mat->data;
   }

   // Get the value at index casted to a type T
   template <class T> const T &at(uint32_t index) const {
      static_assert(std::is_same_v<T, int> || std::is_same_v<T, float> ||
                        std::is_same_v<T, double>,
                    "Incompatible output type.");
      T *data = (T *)mat->data;
      return mat[index];
   }

   // Get the value at index casted to a type T
   template <class T> T &at(uint32_t index) {
      static_assert(std::is_same_v<T, int> || std::is_same_v<T, float> ||
                        std::is_same_v<T, double>,
                    "Incompatible output type.");
      T *data = (T *)mat->data;
      return data[index];
   }

 private:
   uint32_t compute_offset(uint32_t row, uint32_t col) const {
      return row * mat->strides[0] + col * mat->strides[1];
   }

 public:
   // Get the value at row and col casted to a type T
   template <class T> const T &at(uint32_t row, uint32_t col) const {
      static_assert(std::is_same_v<T, int> || std::is_same_v<T, float> ||
                        std::is_same_v<T, double>,
                    "Incompatible output type.");
      T *data = (T *)mat->data;
      return mat[compute_offset(row, col)];
   }

   // Get the value at row and col casted to a type T
   template <class T> T &at(uint32_t row, uint32_t col) {
      static_assert(std::is_same_v<T, int> || std::is_same_v<T, float> ||
                        std::is_same_v<T, double>,
                    "Incompatible output type.");
      T *data = (T *)mat->data;
      return data[compute_offset(row, col)];
   }

   // Elementwise matrix addition
   matrix operator+(const matrix &m) {
      matrix out(mat->type, mat->shape);
      matrice_c::matrix_add(mat, m.ptr(), out.ptr());
      return out;
   }

   // Elementwise matrix substraction
   matrix operator-(const matrix &m) {
      matrix out(mat->type, mat->shape);
      matrice_c::matrix_sub(mat, m.ptr(), out.ptr());
      return out;
   }

   // Elementwise matrix multiplication
   matrix operator*(const matrix &m) {
      matrix out(mat->type, mat->shape);
      matrice_c::matrix_mul(mat, m.ptr(), out.ptr());
      return out;
   }

   // Elementwise matrix division
   matrix operator/(const matrix &m) {
      matrix out(mat->type, mat->shape);
      matrice_c::matrix_div(mat, m.ptr(), out.ptr());
      return out;
   }

   friend std::ostream &operator<<(std::ostream &os, const matrix &);
};

std::ostream &operator<<(std::ostream &os, const matrix &mat) {
   matrice_c::matrix_print(mat.ptr());
   return os;
}

// Matrix multiplication
matrix matmul(const matrix &x, const matrix &y) {
   if (x.type() != y.type()) {
      throw std::runtime_error(
          "Matrix multiplication error, different types.\n");
   }
   if (x.cols() != y.rows()) {
      throw std::runtime_error(
          "Matrix multiplication error, incompatible matrice shapes.\n");
   }
   uint32_t rows = x.rows();
   uint32_t cols = y.cols();
   matrix z(x.type(), rows, cols);
   matrice_c::matrix_matmul(x.ptr(), y.ptr(), z.ptr());
   return z;
}

// Create a matrix of ones
matrix ones(data_type type, uint32_t rows, uint32_t cols) {
   uint32_t shape[2] = {rows, cols};
   matrice_c::matrix *ptr = matrice_c::matrix_ones(type, shape);
   return matrix(ptr);
}

// Create a matrix of zeros
matrix zeros(data_type type, uint32_t rows, uint32_t cols) {
   uint32_t shape[2] = {rows, cols};
   matrice_c::matrix *ptr = matrice_c::matrix_zeros(type, shape);
   return matrix(ptr);
}

// Create a matrix of arange
matrix arange(uint32_t rows, uint32_t cols, int value = 0) {
   uint32_t shape[2] = {rows, cols};
   matrice_c::matrix *ptr =
       matrice_c::matrix_arange(data_type::kint, shape, &value);
   return matrix(ptr);
}
matrix arange(uint32_t rows, uint32_t cols, float value = 0.f) {
   uint32_t shape[2] = {rows, cols};
   matrice_c::matrix *ptr =
       matrice_c::matrix_arange(data_type::kfloat, shape, &value);
   return matrix(ptr);
}
matrix arange(uint32_t rows, uint32_t cols, double value = 0.) {
   uint32_t shape[2] = {rows, cols};
   matrice_c::matrix *ptr =
       matrice_c::matrix_arange(data_type::kdouble, shape, &value);
   return matrix(ptr);
}

} // namespace matrice

#endif
