#include <fstream>
#include <iostream>
#include <lambda_lanczos.hpp>
#include <eigen3/Eigen/Dense>

const auto I_ = std::complex<double>(0.0, 1.0);

const size_t n = 4;
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

void multiply(std::complex<double> mat1[N][N], std::complex<double> mat2[N][N], std::complex<double> res[N][N])
{
  int i, j, k;
  for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
          res[i][j] = 0;
          for (k = 0; k < N; k++)
              res[i][j] += mat1[i][k] * mat2[k][j];
      }
  }
}

void add(std::complex<double> mat1[N][N], std::complex<double> mat2[N][N], std::complex<double> res[N][N])
{
  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
      res[i][j] = mat1[i][j] + mat2[i][j];
  }
}

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
  std::cout << "Enter values of J>0 and h>0" << std::endl;

  double J, h;

  std::cin >> J;
  std::cin >> h;

  if(J<0 || h<0)
  {
    std::cout << "ERROR: Negative Constants" << std::endl;
    return -1;
  }
  
  std::complex<double> H[N][N];

  std::complex<double> sum1[N][N];

  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
        sum1[i][j] = 0;
  }

  std::complex<double> T[N][N];
  std::complex<double> T1[N][N];
  std::complex<double> T2[N][N];

  for(int i=1; i < n; i++)
  {
    kronecker_prod(T1, Z, n-i);
    kronecker_prod(T2, Z, n-i-1);
    multiply(T1, T2, T);
    add(sum1, T, sum1);
  }

  if(n>2)
  {
    kronecker_prod(T1, Z, 0);
    kronecker_prod(T2, Z, n-1);
    multiply(T1, T2, T);
    add(sum1, T, sum1); 
  }

  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
      sum1[i][j] = -J*sum1[i][j];
  }

  std::complex<double> sum2[N][N];

  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
        sum2[i][j] = 0;
  }

  for(int i = 1; i <= n; i++)
  {
    kronecker_prod(T, X, n-i);
    add(sum2, T, sum2);
  }

  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
      sum2[i][j] = h*sum2[i][j];
  }

  add(sum1, sum2, H);

  typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>  Mat;
  
  Mat H_new;

  H_new.resize(N, N);

  for(int i = 0; i<N; i++)
  {
    for(int j = 0; j<N; j++)
    {
      H_new(i, j) = H[i][j];
    }
  }

  Eigen::SelfAdjointEigenSolver<Mat> es(N);

  es.compute(H_new);

  // for(int i = 0; i < N; i++)
  // {
  //   for(int j = 0; j < N; j++)
  //   {
  //     std::cout << H[i][j] << " ";
  //   }
  //   std::cout << std::endl;
  // }

  // Driver Code Below
  auto matmul = [&](const vector<complex<double>>& in, vector<complex<double>>& out) {
      for(size_t i = 0;i < N;i++) {
      for(size_t j = 0;j < N;j++) {
      out[i] += H[i][j]*in[j];
      }
      }
  };

  LambdaLanczos<std::complex<double>> engine(matmul, N, true);
  engine.init_vector = vector_initializer<complex<double>>;
  double eigvalue;
  std::vector<std::complex<double>> eigvec(N);
  engine.run(eigvalue, eigvec);

  std::ofstream outfile;
  outfile.open("Energies.txt");

  outfile << es.eigenvalues();

  std::cout << es.eigenvectors() << std::endl;

  std::cout << "Eigenvalue:" << eigvalue << std::endl;

  return 0;
}