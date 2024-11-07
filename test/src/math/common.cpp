#include"../../../include/fmma/math.hpp"
#include"../test_common.hpp"
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<vector>
#include<array>

//{{{ test_Chebyshev
template<typename TYPE>
bool test_Chebyshev(int n, int size, TYPE tol){
  TYPE diff = 0.0;
  for(int i=0; i<size; ++i){
    TYPE x = (TYPE)rand()*2.0/RAND_MAX-1.0;
    TYPE ans = fmma::Chebyshev(n, x);
    TYPE exact = 1;
    if(n==0){
      exact = 1;
    }else if(n==1){
      exact = x;
    }else if(n==2){
      exact = 2*x*x-1;
    }else if(n==3){
      exact = 4*x*x*x-3*x;
    }else{
      fprintf(stderr, "%s:%d ERROR : Chebyshev test is for n<=3 but %d\n", __FILE__, __LINE__, n);
      exit(EXIT_FAILURE);
    }
    TYPE d = ans - exact;
    diff += d*d;
  }
  diff = sqrt(diff / size);

  if(diff < tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    fprintf(stderr, "diff = %e > tol = %e\n", diff, tol);
    failed(__FILE__, __func__);
    return false;
  }
}
//}}}

//{{{ test_axpy
template<typename TYPE>
bool test_axpy(int n, TYPE tol){
  TYPE diff = 0.0;
  std::vector<TYPE> x(n), y_ans(n), y_exact(n);
  TYPE alpha = (TYPE)rand()/RAND_MAX-0.5;
  for(int i=0; i<n; ++i){
    x[i] = (TYPE)rand()/RAND_MAX-0.5;
    y_ans[i] = (TYPE)rand()/RAND_MAX-0.5;
    y_exact[i] = y_ans[i];
  }

  fmma::axpy(alpha, x, y_ans);
  for(int i=0; i<n; ++i){
    y_exact[i] += alpha*x[i];
  }

  for(int i=0; i<n; ++i){
    TYPE d = y_ans[i] - y_exact[i];
    diff += d*d;
  }
  diff = sqrt(diff / n);

  if(diff < tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    fprintf(stderr, "diff = %e > tol = %e\n", diff, tol);
    failed(__FILE__, __func__);
    return false;
  }
}
//}}}

//{{{ test_dot
template<typename TYPE>
bool test_dot(int n, TYPE tol){
  std::vector<TYPE> x(n), y(n);
  for(int i=0; i<n; ++i){
    x[i] = (TYPE)rand()/RAND_MAX-0.5;
    y[i] = (TYPE)rand()/RAND_MAX-0.5;
  }

  TYPE ans = fmma::dot(x, y);
  TYPE exact = 0.0;
  for(int i=0; i<n; ++i){
    exact += x[i]*y[i];
  }

  TYPE diff = std::fabs(exact-ans);

  if(diff < tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    fprintf(stderr, "diff = %e > tol = %e\n", diff, tol);
    failed(__FILE__, __func__);
    return false;
  }
}
//}}}

//{{{ test_matvec
template<typename TYPE>
bool test_matvec(int n, int m, TYPE tol){
  TYPE diff = 0.0;
  std::vector<TYPE> A(n*m), x(m), y_ans(n), y_exact(n);
  for(int i=0; i<n; ++i){
    for(int j=0; j<m; ++j){
      A[i*m+j] = (TYPE)rand()/RAND_MAX-0.5;
    }
  }
  for(int j=0; j<m; ++j){
    x[j] = (TYPE)rand()/RAND_MAX-0.5;
  }

  fmma::matvec(A, x, y_ans);
  for(int i=0; i<n; ++i){
    y_exact[i] = 0.0;
    for(int j=0; j<m; ++j){
      y_exact[i] += A[i*m+j] * x[j];
    }
  }

  for(int i=0; i<n; ++i){
    TYPE d = y_ans[i] - y_exact[i];
    diff += d*d;
  }
  diff = sqrt(diff / n);

  if(diff < tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    fprintf(stderr, "diff = %e > tol = %e\n", diff, tol);
    std::string comment = "(" + std::to_string(n) + "," + std::to_string(m) + ")";
    failed(__FILE__, __func__, comment);
    return false;
  }
}
//}}}

//{{{ test_solve_gcr
template<typename TYPE>
bool test_solve_gcr(int n, TYPE tol){
  TYPE diff = 0.0;
  std::vector<TYPE> A(n*n), x_ans(n), x_exact(n), b(n);
  for(int i=0; i<n; ++i){
    for(int j=0; j<n; ++j){
      A[i*n+j] = (TYPE)rand()/RAND_MAX-0.5;
    }
  }
  for(int j=0; j<n; ++j){
    x_exact[j] = (TYPE)rand()/RAND_MAX-0.5;
  }

  fmma::matvec(A, x_exact, b);
  fmma::solve_gcr(A, x_ans, b);

  for(int i=0; i<n; ++i){
    TYPE d = x_ans[i] - x_exact[i];
    diff += d*d;
  }
  diff = sqrt(diff / n);

  if(diff < tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    fprintf(stderr, "diff = %e > tol = %e\n", diff, tol);
    failed(__FILE__, __func__);
    return false;
  }
}
//}}}

int main(void){
  srand(0);

  if(!test_Chebyshev<double>(0, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_Chebyshev<double>(1, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_Chebyshev<double>(2, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_Chebyshev<double>(3, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_Chebyshev<float>(0, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_Chebyshev<float>(1, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_Chebyshev<float>(2, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_Chebyshev<float>(3, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_axpy<double>(100, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_axpy<float>(100, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_dot<double>(100, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_dot<float>(100, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_matvec<double>(10, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_matvec<double>(20, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_matvec<double>(10, 20, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_matvec<float>(10, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_matvec<float>(20, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_matvec<float>(10, 20, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_solve_gcr<double>(10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_solve_gcr<float>(10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }

  return 0;
}

