#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <bits/stdc++.h>

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>  Mat;

const auto I_ = std::complex<double>(0.0, 1.0);

const std::complex<double> X[2][2] = 
                    { { 0.0, 1 },
                      { 1, 0.0 } };

const std::complex<double> Y[2][2] = 
                    { { 0.0, -I_ },
                      { I_, 0.0 } };

const std::complex<double> Z[2][2] = 
                    { { 1, 0 },
                      { 0, -1 } };

Mat kronecker_prod(const std::complex<double> mat[2][2], int a, int n)
{
  int b = n - a - 1;

  Mat ans;
  ans.resize(1<<n, 1<<n);

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
  std::ofstream outfile0;
  outfile0.open("params.txt");

  std::cout << "Enter value of n, the lattice site number" << std::endl;

  int n;

  std::cin >> n;

  std::cout << "Enter values of J>0 and h1, h2>0" << std::endl;

  double J, h1, h2;

  std::cin >> J;
  std::cin >> h1;
  std::cin >> h2;

  if(J<0 || h1<0 || h2<0)
  {
    std::cout << "ERROR: Negative Constants" << std::endl;
    return -1;
  }

  Mat T1;
  Mat T2;
  Mat T;
  Mat sum1 = Mat::Zero(1<<n, 1<<n);

  for(int i=1; i < n; i++)
  {
    T1 = kronecker_prod(Z, n-i, n);
    T2 = kronecker_prod(Z, n-i-1, n);
    T = T1*T2;
    sum1 += T;
  }

  if(n>2)
  {
    T1 = kronecker_prod(Z, 0, n);
    T2 = kronecker_prod(Z, n-1, n);
    T = T1*T2;
    sum1 += T;
  }

  sum1 = -J*sum1;

  Mat sum2 = Mat::Zero(1<<n, 1<<n);

  for(int i = 1; i <= n; i++)
  {
    T = kronecker_prod(X, n-i, n);
    sum2 += T;
  }

  sum2 = h1*sum2;
  
  Mat sum3 = Mat::Zero(1<<n, 1<<n);

  for(int i = 1; i<= n; i++)
  {
    T = kronecker_prod(Z, n-i, n);
    sum3 += T;
  }

  sum3 = h2*sum3;

  Mat H = sum1 + sum2 + sum3;

  Eigen::SelfAdjointEigenSolver<Mat> es(1<<n);

  es.compute(H);

  std::ofstream outfile1;
  outfile1.open("Energies.txt");

  outfile1 << es.eigenvalues() << std::endl;

  std::cout << "Enter values of L1, given L = "<< n << std::endl;

  int L1, L2;

  std::cin >> L1;
  L2 = n - L1;

  outfile0 << "n = " << n << std::endl;
  outfile0 << "J = " << J << std::endl;
  outfile0 << "h1 = " << h1 << std::endl;
  outfile0 << "h2 = " << h2 << std::endl;
  outfile0 << "L1 = " << L1 << std::endl;

  std::ofstream outfile2;
  outfile2.open("Entropies.txt");

  double EntropyList[1<<n];

  for(int sn = 0; sn < 1<<n; sn++)
  {
    L2 = n - L1;
    Mat psi = Mat::Zero(1<<L2, 1<<L1);

    for(int i=0; i < 1<<n; i++)
    {
      int r = i%(1<<L1);
      int l = i/(1<<L2);
      psi(l, r) = es.eigenvectors()(sn, i); //First Index is the State Number
    }

    Eigen::BDCSVD<Mat> svd(psi);

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

    for(int i=0; i < 1<<n; i++)
    {
      int r = i%(1<<L1);
      int l = i/(1<<L2);
      psi(l, r) = es.eigenvectors()(0, i); //First Index is the State Number
    }

    Eigen::BDCSVD<Mat> svd(psi);

    Eigen::VectorXd S = svd.singularValues();

    double Entropy = 0;

    for(int i=0; i < svd.singularValues().size(); i++)
    {
      if(S(i)!=0)
        Entropy += -2*(S(i))*(log(S(i)))*(S(i));
    }
    outfile3 << Entropy << std::endl;
  } //Ground State Length Calculations

  // double min = EntropyList[0];
  // double max = EntropyList[0];

  std::vector<int> min_loc;

  double List_avg = 0;
  double List_sq_avg = 0;

  for(int i=0; i < 1<<n; i++)
  {
    // if(EntropyList[i] < min)
    //   min = EntropyList[i];
    // if(EntropyList[i] > max)
    //   max = EntropyList[i];
    List_avg += EntropyList[i];
    List_sq_avg += EntropyList[i]*EntropyList[i];
  }

  List_avg = List_avg/(1<<n);
  List_sq_avg = List_sq_avg/(1<<n);

  double List_std_dev = std::sqrt(List_sq_avg - List_avg*List_avg);

  for(int i=0; i < 1<<n; i++)
  {
    if(EntropyList[i] < List_avg - 2*List_std_dev) //Scar Filter
      min_loc.push_back(i);
  }

  std::cout << "No. of Scars:" << min_loc.size() << std::endl;

  Mat Z_t;
  Z_t.resize(1<<n, 1<<n);

  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> state = Mat::Zero(1<<n, 1);

  for(int i = 0; i < 1<<n; i++)
  {
    for(int scar_n = 0; scar_n < min_loc.size(); scar_n++)
    {
      state(i) += std::exp(std::complex<double>(scar_n, 0)*I_)*es.eigenvectors()(min_loc.at(scar_n), i); //sum of Scar States
    }
  }

  auto normalization = state.dot(state);

  if(normalization.imag() > 0.00001)
    throw "normalization imaginary";

  state = (1/std::sqrt(normalization.real()))*state;

  if(std::abs((state.dot(state)).real() - 1) > 0.0001)
    throw "State not normalized";

  std::cout << "Enter indices of Z(t) and Z(0)" << std::endl;
  int indext, index0;
  std::cin >> indext;
  std::cin >> index0;

  outfile0 << "indext = " << indext << std::endl;
  outfile0 << "index0 = " << index0 << std::endl;

  T1 = kronecker_prod(Z, n-indext, n); //Z_t(0)
  T2 = kronecker_prod(Z, n-index0, n); //Z_0(0)

  std::ofstream outfile4, outfile5;
  outfile4.open("Avg_Z.txt");
  outfile5.open("Avg_Z0.txt");

  int itr = 0;

  for(double t = 0; itr < 100; t += 0.001)
  {
    Mat U = (t*I_*H).exp(); //Careful with unsupported function
    Z_t = U*T1*U.adjoint(); // Time Evolved Operator

    auto avg_Z = state.dot(Z_t*state);
    auto avg_Z0 = state.dot(Z_t*T2*state);
    if(avg_Z.imag() > 0.00001)
      throw "Avg_Z imaginary";

    if(avg_Z0.imag() > 0.00001)
      throw "Avg_Z imaginary";

    outfile4 << avg_Z.real() << std::endl;
    outfile5 << avg_Z0.real() << std::endl;
    itr++;
  }

  return 0;
}