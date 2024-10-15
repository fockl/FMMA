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
void matvec(const std::vector<TYPE>& A, const std::vector<TYPE>& x, std::vector<TYPE>& ans);

template<typename TYPE>
TYPE dot(const std::vector<TYPE>& x, const std::vector<TYPE>& y);

template<typename TYPE>
void axpy(TYPE a, const std::vector<TYPE>& x, std::vector<TYPE>& y);

template<typename TYPE>
void solve_gcr(const std::vector<TYPE>& A, std::vector<TYPE>& x, const std::vector<TYPE>& b, std::size_t ITR = 100);

template<typename TYPE>
  void solve(const std::vector<TYPE>& A, std::vector<TYPE>& x, const std::vector<TYPE>& b);

}; // nemaspace fmma
