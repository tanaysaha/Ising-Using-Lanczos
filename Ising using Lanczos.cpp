#include <eigen3/Eigen/Dense>
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

Mat kronecker_prod(const std::complex<double> mat[2][2], int a)
{
  int b = n - a - 1;

  Mat ans;
  ans.resize(N, N);

  if(a < 0 || b < 0)
    throw "ERROR: Negative exponents";

  if(a==0)
  {
    Mat B = Mat::Identity(1<<b, 1<<b);

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
              ans(t1, t2) = B(i, j) * mat[k][l];
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
    Mat A = Mat::Identity(1<<a, 1<<a);

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
              ans(t1, t2) = mat[i][j] * A(k, l);
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
    Mat B = Mat::Identity(1<<b, 1<<b);
    Mat A = Mat::Identity(1<<a, 1<<a);

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
              temp[t1][t2] = mat[i][j] * A(k, l);
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
              ans(t1, t2) = B(i, j) * temp[k][l];
              t2++;
          } 
        }
        t1++;
        t2=0;
      } 
    } 
  }
  return ans;
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

  Mat T1;
  Mat T2;
  Mat T;
  Mat sum1 = Mat::Zero(N, N);

  for(int i=1; i < n; i++)
  {
    T1 = kronecker_prod(Z, n-i);
    T2 = kronecker_prod(Z, n-i-1);
    T = T1*T2;
    sum1 += T;
  }

  if(n>2)
  {
    T1 = kronecker_prod(Z, 0);
    T2 = kronecker_prod(Z, n-1);
    T = T1*T2;
    sum1 += T;
  }

  sum1 = -J*sum1;

  Mat sum2 = Mat::Zero(N, N);

  for(int i = 1; i <= n; i++)
  {
    T = kronecker_prod(X, n-i);
    sum2 += T;
  }

  sum2 = h*sum2;
  
  Mat H = sum1 + sum2;

  Eigen::SelfAdjointEigenSolver<Mat> es(N);

  es.compute(H);

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
  } //Ground State Length Calculations

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

  // std::cout << min << "locs:" << std::endl;

  // for(int i=0; i<min_loc.size(); i++)
  //   std::cout << min_loc.at(i) << std::endl;

  // std::cout << "Enter time" << std::endl;

  // double t;
  // std::cin >> t;

  Mat Z_t;
  Z_t.resize(N, N);

  Mat state;
  state.resize(N, 1);
  for(int i = 0; i < N; i++)
    state(i) = es.eigenvectors()(0, i);

  std::ofstream outfile4;
  outfile4.open("Avg_Z.txt");

  for(double t = 0; t < std::min( 1/(100*J), 1/(100*h)); t += 0.001)
  {
    Z_t = (Mat::Identity(N, N) + t*I_*H)*T1*(Mat::Identity(N, N) - t*I_*H); // Time Evolved Operator

    auto avg_Z = state.adjoint()*(Z_t*state);
    outfile4 << avg_Z.real() << std::endl;
  }

  return 0;
}