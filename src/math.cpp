#include"../include/fmma/math.hpp"
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<complex>
#include<vector>
#include<array>
#ifdef FMMA_USE_BLAS
#include<cblas.h>
#endif

namespace fmma {

template<typename TYPE>
TYPE Chebyshev(int n, TYPE x){
  if(x < -1.0 || 1.0 < x){
    fprintf(stderr, "%s:%d ERROR : Chebyshev input should be in [-1.0, 1.0] but %lf\n", __FILE__, __LINE__, x);
    exit(EXIT_FAILURE);
  }
  if(n < 0){
    fprintf(stderr, "%s:%d ERROR : Chebyshev dim should be >= 0 but %d\n", __FILE__, __LINE__, n);
    exit(EXIT_FAILURE);
  }
  TYPE t = acos(x);
  return cos(n*t);
}

template double Chebyshev(int n, double x);
template float Chebyshev(int n, float x);

template<typename TYPE>
TYPE SChebyshev(int n, TYPE x, TYPE y){
  if(x < -1.0 || 1.0 < x){
    fprintf(stderr, "%s:%d ERROR : SChebyshev input x should be in [-1.0, 1.0] but %lf\n", __FILE__, __LINE__, x);
    exit(EXIT_FAILURE);
  }
  if(y < -1.0 || 1.0 < y){
    fprintf(stderr, "%s:%d ERROR : SChebyshev input y should be in [-1.0, 1.0] but %lf\n", __FILE__, __LINE__, y);
    exit(EXIT_FAILURE);
  }
  if(n <= 0){
    fprintf(stderr, "%s:%d ERROR : SChebyshev dim should be > 0 but %d\n", __FILE__, __LINE__, n);
    exit(EXIT_FAILURE);
  }
  TYPE ans = (TYPE)1.0/n;
  for(int k=1; k<n; ++k){
    ans += (TYPE)2.0/n*Chebyshev(k, x)*Chebyshev(k, y);
  }
  return ans;
}

template double SChebyshev(int n, double x, double y);
template float SChebyshev(int n, float x, float y);

template<typename TYPE>
void copy(const std::size_t N, const TYPE* a, TYPE* y){
  for(std::size_t n=0; n<N; ++n){
    y[n] = a[n];
  }
  return;
}

template<>
void copy<double>(const std::size_t N, const double* a, double* y){
#if FMMA_USE_BLAS
  cblas_dcopy(N, a, 1, y, 1);
#else
  for(std::size_t n=0; n<N; ++n){
    y[n] = a[n];
  }
#endif
  return;
}

template<>
void copy<float>(const std::size_t N, const float* a, float* y){
#if FMMA_USE_BLAS
  cblas_scopy(N, a, 1, y, 1);
#else
  for(std::size_t n=0; n<N; ++n){
    y[n] = a[n];
  }
#endif
  return;
}

template void copy(const std::size_t N, const double* a, double* y);
template void copy(const std::size_t N, const float* a, float* y);

template<typename TYPE>
TYPE dot(const std::size_t N, const TYPE* x, const TYPE* y){
  TYPE ans = 0.0;
  for(std::size_t n=0; n<N; ++n){
    ans += x[n]*y[n];
  }
  return ans;
}

template<>
double dot<double>(const std::size_t N, const double* x, const double* y){
  double ans = 0.0;
#if FMMA_USE_BLAS
  ans = cblas_ddot(N, x, 1, y, 1);
#else
  for(std::size_t n=0; n<N; ++n){
    ans += x[n]*y[n];
  }
#endif
  return ans;
}

template<>
float dot<float>(const std::size_t N, const float* x, const float* y){
  float ans = 0.0;
#if FMMA_USE_BLAS
  ans = cblas_sdot(N, x, 1, y, 1);
#else
  for(std::size_t n=0; n<N; ++n){
    ans += x[n]*y[n];
  }
#endif
  return ans;
}

template double dot(const std::size_t N, const double* x, const double* y);
template float dot(const std::size_t N, const float* x, const float* y);

template<typename TYPE>
void matmul(const std::size_t M, const std::size_t N, const std::size_t K, const TYPE alpha, const std::vector<TYPE>& A, const std::vector<TYPE>& B, const TYPE beta, std::vector<TYPE>& C){
  // matrix-matrix multiplication (AB = C // MK * KN = MN)
  if(A.size() != M*K){
    fprintf(stderr, "%s:%d ERROR : matrix size error (A:%zu != M(%zu)*K(%zu))\n", __FILE__, __LINE__, A.size(), M, K);
    exit(EXIT_FAILURE);
  }
  if(B.size() != K*N){
    fprintf(stderr, "%s:%d ERROR : matrix size error (B:%zu != K(%zu)*N(%zu))\n", __FILE__, __LINE__, B.size(), K, N);
    exit(EXIT_FAILURE);
  }
  C.resize(M*N);

  for(std::size_t m=0; m<M; ++m){
    for(std::size_t n=0; n<N; ++n){
      C[m*N+n] = beta*C[m*N+n];
    }
    for(std::size_t k=0; k<K; ++k){
      for(std::size_t n=0; n<N; ++n){
        C[m*N+n] += alpha*A[m*K+k]*B[k*N+n];
      }
    }
  }

  return;
}


template<>
void matmul<double>(const std::size_t M, const std::size_t N, const std::size_t K, const double alpha, const std::vector<double>& A, const std::vector<double>& B, const double beta, std::vector<double>& C){
  // matrix-matrix multiplication (AB = C // MK * KN = MN)
  if(A.size() != M*K){
    fprintf(stderr, "%s:%d ERROR : matrix size error (A:%zu != M(%zu)*K(%zu))\n", __FILE__, __LINE__, A.size(), M, K);
    exit(EXIT_FAILURE);
  }
  if(B.size() != K*N){
    fprintf(stderr, "%s:%d ERROR : matrix size error (B:%zu != K(%zu)*N(%zu))\n", __FILE__, __LINE__, B.size(), K, N);
    exit(EXIT_FAILURE);
  }
  C.resize(M*N);

#if FMMA_USE_BLAS
  const auto *Ad = A.data();
  const auto *Bd = B.data();
  auto *Cd = C.data();

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, Ad, K, Bd, N, beta, Cd, N);
#else
  for(std::size_t m=0; m<M; ++m){
    for(std::size_t n=0; n<N; ++n){
      C[m*N+n] = beta*C[m*N+n];
    }
    for(std::size_t k=0; k<K; ++k){
      for(std::size_t n=0; n<N; ++n){
        C[m*N+n] += alpha*A[m*K+k]*B[k*N+n];
      }
    }
  }

#endif

  return;
}


template<>
void matmul<float>(const std::size_t M, const std::size_t N, const std::size_t K, const float alpha, const std::vector<float>& A, const std::vector<float>& B, const float beta, std::vector<float>& C){
  // matrix-matrix multiplication (AB = C // MK * KN = MN)
  if(A.size() != M*K){
    fprintf(stderr, "%s:%d ERROR : matrix size error (A:%zu != M(%zu)*K(%zu))\n", __FILE__, __LINE__, A.size(), M, K);
    exit(EXIT_FAILURE);
  }
  if(B.size() != K*N){
    fprintf(stderr, "%s:%d ERROR : matrix size error (B:%zu != K(%zu)*N(%zu))\n", __FILE__, __LINE__, B.size(), K, N);
    exit(EXIT_FAILURE);
  }
  C.resize(M*N);

#if FMMA_USE_BLAS
  const auto *Ad = A.data();
  const auto *Bd = B.data();
  auto *Cd = C.data();

  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, Ad, K, Bd, N, beta, Cd, N);
#else
  for(std::size_t m=0; m<M; ++m){
    for(std::size_t n=0; n<N; ++n){
      C[m*N+n] = beta*C[m*N+n];
    }
    for(std::size_t k=0; k<K; ++k){
      for(std::size_t n=0; n<N; ++n){
        C[m*N+n] += alpha*A[m*K+k]*B[k*N+n];
      }
    }
  }

#endif

  return;
}

template void matmul(const std::size_t M, const std::size_t N, const std::size_t K, const double alpha, const std::vector<double>& A, const std::vector<double>& B, const double beta, std::vector<double>& C);
template void matmul(const std::size_t M, const std::size_t N, const std::size_t K, const float alpha, const std::vector<float>& A, const std::vector<float>& B, const float beta, std::vector<float>& C);

template<typename TYPE>
void matmul(const std::size_t M, const std::size_t N, const std::size_t K, const std::vector<TYPE>& A, const std::vector<TYPE>& B, std::vector<TYPE>& C){
  matmul(M, N, K, (TYPE)1.0, A, B, (TYPE)0.0, C);
  return;
}

template void matmul(const std::size_t M, const std::size_t N, const std::size_t K, const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C);
template void matmul(const std::size_t M, const std::size_t N, const std::size_t K, const std::vector<float>& A, const std::vector<float>& B, std::vector<float>& C);

template<typename TYPE>
void matvec(const TYPE alpha, const std::vector<TYPE>& A, const std::vector<TYPE>& x, const TYPE beta, std::vector<TYPE>& ans){
  // matrix-vector multiplication (ans = Ax)
  std::size_t N = x.size();
  std::size_t M = A.size()/N;
  if(A.size() != N*M){
    fprintf(stderr, "%s:%d ERROR : matvec size error (A:%zu, x:%zu)\n", __FILE__, __LINE__, A.size(), x.size());
    exit(EXIT_FAILURE);
  }
  ans.resize(M);

  for(std::size_t m=0; m<M; ++m){
    ans[m] = beta*ans[m];
    for(std::size_t n=0; n<N; ++n){
      ans[m] += alpha*A[m*N+n]*x[n];
    }
  }

  return;
}


template<>
void matvec<double>(const double alpha, const std::vector<double>& A, const std::vector<double>& x, const double beta, std::vector<double>& ans){
  // matrix-vector multiplication (ans = Ax)
  std::size_t N = x.size();
  std::size_t M = A.size()/N;
  if(A.size() != N*M){
    fprintf(stderr, "%s:%d ERROR : matvec size error (A:%zu, x:%zu)\n", __FILE__, __LINE__, A.size(), x.size());
    exit(EXIT_FAILURE);
  }
  ans.resize(M);

#if FMMA_USE_BLAS
  const auto *xd = x.data();
  auto *yd = ans.data();
  const auto *vald = A.data();
  const auto m = M;
  const auto n = N;

  cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, alpha, vald, n, xd, 1, beta, yd, 1);

#else
  for(std::size_t m=0; m<M; ++m){
    ans[m] = beta*ans[m];
    for(std::size_t n=0; n<N; ++n){
      ans[m] += alpha*A[m*N+n]*x[n];
    }
  }

#endif

  return;
}

template<>
void matvec<float>(const float alpha, const std::vector<float>& A, const std::vector<float>& x, const float beta, std::vector<float>& ans){
  // matrix-vector multiplication (ans = Ax)
  std::size_t N = x.size();
  std::size_t M = A.size()/N;
  if(A.size() != N*M){
    fprintf(stderr, "%s:%d ERROR : matvec size error (A:%zu, x:%zu)\n", __FILE__, __LINE__, A.size(), x.size());
    exit(EXIT_FAILURE);
  }
  ans.resize(M);

#if FMMA_USE_BLAS
  const auto *xd = x.data();
  auto *yd = ans.data();
  const auto *vald = A.data();
  const auto m = M;
  const auto n = N;

  cblas_sgemv(CblasRowMajor, CblasNoTrans, m, n, alpha, vald, n, xd, 1, beta, yd, 1);

#else
  for(std::size_t m=0; m<M; ++m){
    ans[m] = beta*ans[m];
    for(std::size_t n=0; n<N; ++n){
      ans[m] += alpha*A[m*N+n]*x[n];
    }
  }

#endif

  return;
}

template void matvec(const double alpha, const std::vector<double>& A, const std::vector<double>& x, const double beta, std::vector<double>& ans);
template void matvec(const float alpha, const std::vector<float>& A, const std::vector<float>& x, const float beta, std::vector<float>& ans);

template<typename TYPE>
void matvec(const std::vector<TYPE>& A, const std::vector<TYPE>& x, std::vector<TYPE>& ans){
  matvec((TYPE)1.0, A, x, (TYPE)0.0, ans);
  return;
}

template void matvec(const std::vector<double>& A, const std::vector<double>& x, std::vector<double>& ans);
template void matvec(const std::vector<float>& A, const std::vector<float>& x, std::vector<float>& ans);

template<typename TYPE>
TYPE dot(const std::vector<TYPE>& x, const std::vector<TYPE>& y){
  // inner production
  std::size_t N = x.size();
  if(y.size() != N){
    fprintf(stderr, "%s:%d ERROR : dot size error x(%zu)!=y(%zu)\n", __FILE__, __LINE__, x.size(), y.size());
    exit(EXIT_FAILURE);
  }
  return dot(N, x.data(), y.data());
}

template double dot(const std::vector<double>& x, const std::vector<double>& y);
template float dot(const std::vector<float>& x, const std::vector<float>& y);

template<typename TYPE>
void axpy(TYPE a, const std::vector<TYPE>& x, std::vector<TYPE>& y){
  std::size_t N = x.size();
  if(y.size() != N){
    fprintf(stderr, "%s:%d ERROR : dot size error x(%zu)!=y(%zu)\n", __FILE__, __LINE__, x.size(), y.size());
    exit(EXIT_FAILURE);
  }
  for(std::size_t i=0; i<N; ++i){
    y[i] += a*x[i];
  }
  return;
}

template void axpy(double a, const std::vector<double>& x, std::vector<double>& y);
template void axpy(float a, const std::vector<float>& x, std::vector<float>& y);

template<typename TYPE>
void solve_gcr(const std::vector<TYPE>& A, std::vector<TYPE>& x, const std::vector<TYPE>& b, std::size_t ITR){
  std::size_t N = b.size();
  if(A.size() != N*N){
    fprintf(stderr, "%s:%d ERROR : A size error A(%zu,%zu), b(%zu)\n", __FILE__, __LINE__, A.size()/N, N, N);
    exit(EXIT_FAILURE);
  }

  std::vector<TYPE> r = b;
  std::vector<TYPE> Ax;
  x.resize(N);
  for(std::size_t i=0; i<N; ++i) x[i] = ((TYPE)rand()/RAND_MAX-0.5)/N;
  matvec(A, x, Ax);
  for(std::size_t i=0; i<N; ++i) r[i] -= Ax[i];
  std::vector<TYPE> p_new = r;

  std::vector<std::vector<TYPE>> p;
  std::vector<std::vector<TYPE>> Ap;
  std::vector<TYPE> Ap_dot;

  ITR = std::min(ITR, N);

  {
    TYPE resid = dot(r, r);
    fprintf(stderr, "resid = %lf\n", resid);
    if(resid < 1.0e-6){
      return;
    }
  }

  for(std::size_t itr=0; itr<ITR; ++itr){
    fprintf(stderr, "itr : %zu\n", itr);
    std::vector<TYPE> Ap_new;
    matvec(A, p_new, Ap_new);
    p.push_back(p_new);
    Ap.push_back(Ap_new);
    TYPE tmp = dot(Ap_new, Ap_new);
    Ap_dot.push_back(tmp);

    TYPE alpha = dot(r, Ap_new)/tmp;
    axpy(alpha, p_new, x);
    axpy(-alpha, Ap_new, r);

    TYPE resid = dot(r, r);
    fprintf(stderr, "resid = %lf   (alpha = %lf)\n", resid, alpha);
    if(resid < 1.0e-6){
      break;
    }

    std::vector<TYPE> Ar_new;
    matvec(A, r, Ar_new);
    std::vector<TYPE> beta(itr+1);
    for(std::size_t k=0; k<N; ++k){
      p_new[k] = r[k];
    }
    for(std::size_t i=0; i<=itr; ++i){
      beta[i] = -dot(Ar_new, Ap[i])/Ap_dot[i];
      axpy(beta[i], p[i], p_new);
    }
  }

  return;
}

template void solve_gcr(const std::vector<double>& A, std::vector<double>& x, const std::vector<double>& b, std::size_t ITR);
template void solve_gcr(const std::vector<float>& A, std::vector<float>& x, const std::vector<float>& b, std::size_t ITR);

template<typename TYPE>
void solve(const std::vector<TYPE>& A, std::vector<TYPE>& x, const std::vector<TYPE>& b){
  // solve linear equation (Ax = b)
  solve_gcr(A, x, b);
  return;
}

template void solve(const std::vector<double>& A, std::vector<double>& x, const std::vector<double>& b);
template void solve(const std::vector<float>& A, std::vector<float>& x, const std::vector<float>& b);

}; // namespace fmma
