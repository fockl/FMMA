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

template<typename TYPE, std::size_t DIM>
FMMA<TYPE, DIM>::~FMMA(void){
  return;
};

template FMMA<double, 1>::~FMMA(void);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::exact(std::vector<std::array<TYPE, DIM>>& source, std::vector<std::array<TYPE, DIM>>& target, std::vector<TYPE>& ans){
  if(ans.size() != target.size()){
    fprintf(stderr, "%s:%d ERROR : target size %zu != ans size %zu\n", __FILE__, __LINE__, target.size(), ans.size());
    exit(EXIT_FAILURE);
  }
  for(std::size_t i=0; i<ans.size(); ++i){
    ans[i] = 0.0;
    for(std::size_t j=0; j<source.size(); ++j){
      ans[i] += fn(source[j], target[i]);
    }
  }
  return;
};

template void FMMA<double, 1>::exact(std::vector<std::array<double, 1>>& source, std::vector<std::array<double, 1>>& target, std::vector<double>& ans);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::solve(std::vector<std::array<TYPE, DIM>>& source, std::vector<std::array<TYPE, DIM>>& target, std::vector<TYPE>& ans){
  if(this->solve_type == "exact"){
    exact(source, target, ans);
  }else{
    fprintf(stderr, "%s:%d ERROR : solve type %s not undefined\n", __FILE__, __LINE__, this->solve_type.c_str());
    exit(EXIT_FAILURE);
  }
  return;
};

template void FMMA<double, 1>::solve(std::vector<std::array<double, 1>>& source, std::vector<std::array<double, 1>>& target, std::vector<double>& ans);

} // namespace fmma
