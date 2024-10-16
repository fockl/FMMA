#include<cstdlib>
#include<cstdio>
#include<vector>
#include<array>
#include<functional>
#include<cmath>
#include"../include/fmma/fmma.hpp"

namespace fmma {

template<typename TYPE, std::size_t DIM>
FMMA<TYPE, DIM>::FMMA(void){
  return;
};

template FMMA<double, 1>::FMMA(void);
template FMMA<double, 2>::FMMA(void);

template<typename TYPE, std::size_t DIM>
FMMA<TYPE, DIM>::~FMMA(void){
  return;
};

template FMMA<double, 1>::~FMMA(void);
template FMMA<double, 2>::~FMMA(void);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::exact(const std::vector<std::array<TYPE, DIM>>& source, const std::vector<std::array<TYPE, DIM>>& target, std::vector<TYPE>& ans){
  std::size_t N = target.size();
  ans.resize(N);
  for(std::size_t i=0; i<N; ++i){
    ans[i] = 0.0;
    for(std::size_t j=0; j<source.size(); ++j){
      ans[i] += fn(source[j], target[i]);
    }
  }
  return;
};

template void FMMA<double, 1>::exact(const std::vector<std::array<double, 1>>& source, const std::vector<std::array<double, 1>>& target, std::vector<double>& ans);
template void FMMA<double, 2>::exact(const std::vector<std::array<double, 2>>& source, const std::vector<std::array<double, 2>>& target, std::vector<double>& ans);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::solve(const std::vector<std::array<TYPE, DIM>>& source, const std::vector<std::array<TYPE, DIM>>& target, std::vector<TYPE>& ans){
  if(this->solve_type == "exact"){
    exact(source, target, ans);
  }else if(this->solve_type == "nrnmm"){
    nrnmm(source, target, ans);
  }else{
    fprintf(stderr, "%s:%d ERROR : solve type %s not undefined\n", __FILE__, __LINE__, this->solve_type.c_str());
    exit(EXIT_FAILURE);
  }
  return;
};

template void FMMA<double, 1>::solve(const std::vector<std::array<double, 1>>& source, const std::vector<std::array<double, 1>>& target, std::vector<double>& ans);
template void FMMA<double, 2>::solve(const std::vector<std::array<double, 2>>& source, const std::vector<std::array<double, 2>>& target, std::vector<double>& ans);

} // namespace fmma
