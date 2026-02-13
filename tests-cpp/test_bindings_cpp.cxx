#include "../bindings/cpp/matrice.hxx"

int main() {
   using namespace matrice;

   {
      vector x = arange(10, 1.0f);
      vector y = arange(10, 1.0f);
      auto z = x + y;
      std::cout << z;
      std::cout << z[0] << std::endl;
      std::cout << z(0) << std::endl;
   }

   {
      matrix x = arange(3, 3, 1.0f);
      matrix y = arange(3, 3, 1.0f);
      std::cout << x;
      std::cout << y;
      auto z = x + y;
      std::cout << z;
   }

   {
      matrix x = arange(3, 3, 1.0f);
      matrix y = arange(3, 3, 1.0f);
      auto z = matmul(x, y);
      std::cout << z;
   }

   {
      vector x({1.0f, 2.0f, 3.0f});
      std::cout << x;
   }

   {
      matrix x({{1.0f, 2.0f}, {3.0f, 4.0f}, {5.0f, 6.0f}});
      std::cout << x;
   }
}

