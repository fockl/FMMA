#pragma once

#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<complex>
#include<vector>
#include<array>

namespace fmma {

template<typename TYPE>
TYPE Chebyshev(int n, TYPE x);

template<typename TYPE>
TYPE SChebyshev(int n, TYPE x, TYPE y);

template<typename TYPE>
void copy(const std::size_t N, const TYPE* a, TYPE* y);

template<typename TYPE>
TYPE dot(const std::size_t N, const TYPE* x, const TYPE* y);

template<typename TYPE>
void matmul(const std::size_t M, const std::size_t N, const std::size_t K, const TYPE alpha, const std::vector<TYPE>& A, const std::vector<TYPE>& B, const TYPE beta, std::vector<TYPE>& C);

template<typename TYPE>
void matmul(const std::size_t M, const std::size_t N, const std::size_t K, const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C);

template<typename TYPE>
void matvec(const std::vector<TYPE>& A, const std::vector<TYPE>& x, std::vector<TYPE>& ans);

template<typename TYPE>
void matvec(const TYPE alpha, const std::vector<TYPE>& A, const std::vector<TYPE>& x, const TYPE beta, std::vector<TYPE>& ans);

template<typename TYPE>
TYPE dot(const std::vector<TYPE>& x, const std::vector<TYPE>& y);

template<typename TYPE>
void axpy(TYPE a, const std::vector<TYPE>& x, std::vector<TYPE>& y);

template<typename TYPE>
void solve_gcr(const std::vector<TYPE>& A, std::vector<TYPE>& x, const std::vector<TYPE>& b, std::size_t ITR = 100);

template<typename TYPE>
  void solve(const std::vector<TYPE>& A, std::vector<TYPE>& x, const std::vector<TYPE>& b);

}; // nemaspace fmma
