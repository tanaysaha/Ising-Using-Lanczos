#include <fstream>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <bits/stdc++.h>

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>  Mat;

const auto I_ = std::complex<double>(0.0, 1.0);

const size_t n = 8;
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

  // Shows the Hamiltonian
  // for(int i = 0; i < N; i++)
  // {
  //   for(int j = 0; j < N; j++)
  //   {
  //     std::cout << H[i][j] << " ";
  //   }
  //   std::cout << std::endl;
  // }

  std::ofstream outfile1;
  outfile1.open("Energies.txt");

  outfile1 << es.eigenvalues();

  std::cout << "Enter values of L1, given L = "<< n << std::endl;

  int L1, L2;

  std::cin >> L1;
  L2 = n - L1;

  std::ofstream outfile2;
  outfile2.open("Entropies.txt");

  double EntropyList[N];

  for(int sn = 0; sn < N; sn++)
  {
    L2 = n - L1;
    Mat psi = Mat::Zero(1<<L2, 1<<L1);

    for(int i=0; i < N; i++)
    {
      int r = i%(1<<L1);
      int l = i/(1<<L2);
      psi(l, r) = es.eigenvectors()(sn, i); //First Index is the State Number
    }

    Eigen::JacobiSVD<Mat> svd(psi, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::VectorXd S = svd.singularValues();

    double Entropy = 0;

    for(int i=0; i < svd.singularValues().size(); i++)
    {
      if(S(i)!=0)
        Entropy += -2*(S(i))*(log(S(i)))*(S(i));
    }
    EntropyList[sn] = Entropy;
    outfile2 << Entropy << std::endl;
  }

  std::ofstream outfile3;
  outfile3.open("Entropies_gs.txt");

  for(L1 = 0; L1 <= n/2; L1++)
  {
    L2 = n - L1;
    Mat psi = Mat::Zero(1<<L2, 1<<L1);

    for(int i=0; i < N; i++)
    {
      int r = i%(1<<L1);
      int l = i/(1<<L2);
      psi(l, r) = es.eigenvectors()(0, i); //First Index is the State Number
    }

    Eigen::BDCSVD<Mat> svd(psi, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::VectorXd S = svd.singularValues();

    double Entropy = 0;

    for(int i=0; i < svd.singularValues().size(); i++)
    {
      if(S(i)!=0)
        Entropy += -2*(S(i))*(log(S(i)))*(S(i));
    }
    outfile3 << Entropy << std::endl;
  }

  double min = EntropyList[0];
  int min_loc_temp = 0;

  std::vector<int> min_loc;

  for(int i=0; i < N; i++)
  {
    if(EntropyList[i] < min)
    {
      min = EntropyList[i];
      min_loc_temp = i;
    }
  }

  for(int i=0; i < N; i++)
  {
    if(EntropyList[i] > 0.99*min && EntropyList[i] < 1.01*min)
      min_loc.push_back(i);
  }

  std::cout << min << "locs:" << std::endl;

  for(int i=0; i<min_loc.size(); i++)
    std::cout << min_loc.at(i) << std::endl;

  Mat Z_t;
  Z_t.resize(N, N);

  Mat T_new;
  T_new.resize(N, N);

  for(int i = 0; i<N; i++)
  {
    for(int j = 0; j<N; j++)
    {
      T_new(i, j) = T1[i][j];
    }
  }

  Z_t = (Mat::Identity(N, N) + I_*H_new)*T_new*(Mat::Identity(N, N) - I_*H_new);
  Mat state;
  state.resize(N, 1);
  for(int i = 0; i < N; i++)
    state(i) = es.eigenvectors()(0, i);

  auto avg_Z = state.adjoint()*(Z_t*state);

  std::cout << avg_Z << std::endl;

  return 0;
}