#include "../bindings/cpp/matrice.hxx"

int main() {
   using namespace matrice;

   {
      vector x = arange(10, 1.0f);
      vector y = arange(10, 1.0f);
      auto z = x + y;
      std::cout << z;
      std::cout << z.at<float>(0) << std::endl;
      std::cout << z.data<float>()[0] << std::endl;
   }

   {
      matrix x = arange(3, 3, 1.0f);
      matrix y = arange(3, 3, 1.0f);
      auto z = x + y;
      std::cout << z;
   }

   {
      matrix x = arange(3, 3, 1.0f);
      matrix y = arange(3, 3, 1.0f);
      auto z = matmul(x, y);
      std::cout << z;
   }

}

