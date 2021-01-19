#include <iostream>
#include <lambda_lanczos.hpp>

const auto I_ = std::complex<double>(0.0, 1.0);

const size_t n = 2;
const size_t N = 1<<n;

const std::complex<double> X[2][2] = 
                    { { 0.0, 1 },
                      { 1, 0.0 } };

const std::complex<double> Y[2][2] = 
                    { { 0.0, -I_ },
                      { I_, 0.0 } };

const std::complex<double> Z[2][2] = 
                    { { 1, 0 },
                      { 0, -1 } };

// Remember to remove this \/\/\/\/\/\/\/\/\/\/\/

using lambda_lanczos::LambdaLanczos;

template<typename T>
using vector = std::vector<T>;

template<typename T>
using complex = std::complex<T>;


template <typename T>
void vector_initializer(vector<T>& v);

template<>
void vector_initializer(vector<double>& v) {
  std::mt19937 mt(1);
  std::uniform_real_distribution<double> rand(-1.0, 1.0);

  size_t n = v.size();
  for(size_t i = 0;i < n;i++) {
    v[i] = rand(mt);
  }
}

template<>
void vector_initializer(vector<complex<double>>& v) {
  std::mt19937 mt(1);
  std::uniform_real_distribution<double> rand(-1.0, 1.0);

  size_t n = v.size();
  for(size_t i = 0;i < n;i++) {
    v[i] = std::complex<double>(rand(mt), rand(mt));
  }
}

// Remember to remove this^^^^^^^^^^^^^^^^^^

void kronecker_prod(std::complex<double> ans[N][N], const std::complex<double> mat[2][2], int a)
{
  int b = n - a - 1;

  if(a < 0 || b < 0)
  {
    std::cout << "ERROR: Negative exponents" << std::endl;
    return;
  }

  if(a==0)
  {
    std::complex<double> B[1<<b][1<<b];
    for(int i = 0; i < 1<<b; i++)
    {
      for(int j = 0; j < 1<<b; j++)
      {
        if(i==j)
          B[i][j] = 1;
        else
          B[i][j] = 0;
      }
    }

    int t1=0, t2=0;
    for (int i = 0; i < 1<<b; i++) { 

      // k loops till rowb 
      for (int k = 0; k < 2*1<<a; k++) { 

        // j loops till cola 
        for (int j = 0; j < 1<<b; j++) { 

          // l loops till colb 
          for (int l = 0; l < 2*1<<a; l++) { 

              // Each element of matrix A is 
              // multiplied by whole Matrix B 
              // resp and stored as Matrix C 
              ans[t1][t2] = B[i][j] * mat[k][l];
              t2++;
          } 
        }
        t1++;
        t2=0;
      } 
    } 
  }

  else if(b==0)
  {
    std::complex<double> A[1<<a][1<<a];
    for(int i = 0; i < 1<<a; i++)
    {
      for(int j = 0; j < 1<<a; j++)
      {
        if(i==j)
          A[i][j] = 1;
        else
          A[i][j] = 0;
      }
    }

    int t1=0, t2=0;
    for (int i = 0; i < 2; i++) { 

      // k loops till rowb 
      for (int k = 0; k < 1<<a; k++) { 

        // j loops till cola 
        for (int j = 0; j < 2; j++) { 

          // l loops till colb 
          for (int l = 0; l < 1<<a; l++) { 

              // Each element of matrix A is 
              // multiplied by whole Matrix B 
              // resp and stored as Matrix C 
              ans[t1][t2] = mat[i][j] * A[k][l];
              t2++;
          } 
        }
        t1++;
        t2 = 0;
      } 
    } 
  }
  else
  {
    std::complex<double> A[1<<a][1<<a];
    std::complex<double> B[1<<b][1<<b];

    for(int i = 0; i < 1<<a; i++)
    {
      for(int j = 0; j < 1<<a; j++)
      {
        if(i==j)
          A[i][j] = 1;
        else
          A[i][j] = 0;
      }
    }

    for(int i = 0; i < 1<<b; i++)
    {
      for(int j = 0; j < 1<<b; j++)
      {
        if(i==j)
          B[i][j] = 1;
        else
          B[i][j] = 0;
      }
    }

    std::complex<double> temp[2*1<<a][2*1<<a];

    int t1=0, t2=0;

    for (int i = 0; i < 2; i++) { 

      // k loops till rowb 
      for (int k = 0; k < 1<<a; k++) { 

        // j loops till cola 
        for (int j = 0; j < 2; j++) { 

          // l loops till colb 
          for (int l = 0; l < 1<<a; l++) { 

              // Each element of matrix A is 
              // multiplied by whole Matrix B 
              // resp and stored as Matrix C 
              temp[t1][t2] = mat[i][j] * A[k][l];
              t2++;
          } 
        }
        t1++;
        t2 = 0;
      } 
    } 

    t1=0;
    for (int i = 0; i < 1<<b; i++) { 

      // k loops till rowb 
      for (int k = 0; k < 2*1<<a; k++) { 

        // j loops till cola 
        for (int j = 0; j < 1<<b; j++) { 

          // l loops till colb 
          for (int l = 0; l < 2*1<<a; l++) { 

              // Each element of matrix A is 
              // multiplied by whole Matrix B 
              // resp and stored as Matrix C 
              ans[t1][t2] = B[i][j] * temp[k][l];
              t2++;
          } 
        }
        t1++;
        t2=0;
      } 
    } 
  }
}

int main()
{
  std::cout << "Enter values of J and h" << std::endl;

  double J, h;

  std::cin >> J;
  std::cin >> h;

  // std::complex<double> matrix[n][n] = { { 0.0, -I_ },
  //                 { I_, 0.0 } };
  // /* Its eigenvalues are {-2, 1, 1} */

  // auto matmul = [&](const vector<complex<double>>& in, vector<complex<double>>& out) {
  //     for(size_t i = 0;i < n;i++) {
  //     for(size_t j = 0;j < n;j++) {
  //     out[i] += matrix[i][j]*in[j];
  //     }
  //     }
  // };

  // // Driver Code Below
  // LambdaLanczos<std::complex<double>> engine(matmul, n);
  // engine.init_vector = vector_initializer<complex<double>>;
  // double eigvalue;
  // std::vector<std::complex<double>> eigvec(n);
  // engine.run(eigvalue, eigvec);

  // std::cout << eigvalue << std::endl;
  
  // for(size_t i = 0;i < n; i++) {
  //     std::cout << eigvec[i] << " ";
  // }
  
  std::complex<double> H[N][N];

  


  // std::cout << N << std::endl;

  // kronecker_prod(ans, X, 1);

  // for(int i = 0; i < N; i++)
  // {
  //   for(int j = 0; j < N; j++)
  //   {
  //     std::cout << std::real(ans[i][j]) << " ";
  //   }
  //   std::cout << std::endl;
  // }

  return 0;
}