#include"../include/fmma/fmma.hpp"
#include"../include/fmma/math.hpp"
#include<cstdlib>
#include<cstdio>
#include<vector>
#include<array>
#include<functional>
#include<cmath>
#include<chrono>
#include<algorithm>

namespace fmma {

template<typename TYPE, std::size_t DIM>
FMMA<TYPE, DIM>::FMMA(void){
  return;
};

template FMMA<double, 1>::FMMA(void);
template FMMA<double, 2>::FMMA(void);
template FMMA<double, 3>::FMMA(void);

template<typename TYPE, std::size_t DIM>
FMMA<TYPE, DIM>::~FMMA(void){
  return;
};

template FMMA<double, 1>::~FMMA(void);
template FMMA<double, 2>::~FMMA(void);
template FMMA<double, 3>::~FMMA(void);

template<typename TYPE, std::size_t DIM>
bool FMMA<TYPE, DIM>::check_blas(void){
#if FMMA_USE_BLAS
  return true;
#else
  return false;
#endif
};

template bool FMMA<double, 1>::check_blas(void);
template bool FMMA<double, 2>::check_blas(void);
template bool FMMA<double, 3>::check_blas(void);


template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::exact(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){
  std::size_t N = target.size();
  ans.resize(N);
  for(std::size_t i=0; i<N; ++i){
    ans[i] = 0.0;
    for(std::size_t j=0; j<source.size(); ++j){
      ans[i] += fn(target[i]-source[j]);
    }
  }
  return;
};

template void FMMA<double, 1>::exact(const std::vector<std::array<double, 1>>& target, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::exact(const std::vector<std::array<double, 2>>& target, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);
template void FMMA<double, 3>::exact(const std::vector<std::array<double, 3>>& target, const std::vector<std::array<double, 3>>& source, std::vector<double>& ans);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::exact(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){
  std::size_t N = target.size();
  ans.resize(N);
  for(std::size_t i=0; i<N; ++i){
    ans[i] = 0.0;
    for(std::size_t j=0; j<source.size(); ++j){
      ans[i] += source_weight[j] * fn(target[i]-source[j]);
    }
  }
  return;
};

template void FMMA<double, 1>::exact(const std::vector<std::array<double, 1>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::exact(const std::vector<std::array<double, 2>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);
template void FMMA<double, 3>::exact(const std::vector<std::array<double, 3>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 3>>& source, std::vector<double>& ans);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::exact_matvec(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){
  std::size_t N = target.size();
  std::size_t M = source.size();
  ans.resize(N);
  std::vector<TYPE> Mat(N*M);
  for(std::size_t i=0; i<N; ++i){
    for(std::size_t j=0; j<M; ++j){
      Mat[i*M+j] = fn(target[i]-source[j]);
    }
  }
  matvec(Mat, source_weight, ans);
  return;
};

template void FMMA<double, 1>::exact_matvec(const std::vector<std::array<double, 1>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::exact_matvec(const std::vector<std::array<double, 2>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);
template void FMMA<double, 3>::exact_matvec(const std::vector<std::array<double, 3>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 3>>& source, std::vector<double>& ans);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::exact_matvec(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){
  std::size_t M = source.size();
  std::vector<TYPE> ones(M);
  for(std::size_t j=0; j<M; ++j){
    ones[j] = 1.0;
  }
  exact_matvec(target, ones, source, ans);
  return;
};

template void FMMA<double, 1>::exact_matvec(const std::vector<std::array<double, 1>>& target, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::exact_matvec(const std::vector<std::array<double, 2>>& target, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);
template void FMMA<double, 3>::exact_matvec(const std::vector<std::array<double, 3>>& target, const std::vector<std::array<double, 3>>& source, std::vector<double>& ans);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::solve(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){
  if(this->solve_type == "exact"){
    exact(target, source, ans);
  }else if(this->solve_type == "exact_matvec"){
    exact_matvec(target, source, ans);
  }else if(this->solve_type == "nrnmm"){
    nrnmm(target, source, ans);
  }else if(this->solve_type == "tree"){
    tree(target, source, ans);
  }else if(this->solve_type == "fmm"){
    fmm(target, source, ans);
  }else{
    fprintf(stderr, "%s:%d ERROR : solve type %s not undefined\n", __FILE__, __LINE__, this->solve_type.c_str());
    exit(EXIT_FAILURE);
  }
  return;
};

template void FMMA<double, 1>::solve(const std::vector<std::array<double, 1>>& target, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::solve(const std::vector<std::array<double, 2>>& target, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);
template void FMMA<double, 3>::solve(const std::vector<std::array<double, 3>>& target, const std::vector<std::array<double, 3>>& source, std::vector<double>& ans);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::solve(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){
  std::chrono::system_clock::time_point start, end;
  start = std::chrono::system_clock::now();
  if(this->solve_type == "exact"){
    exact(target, source_weight, source, ans);
  }else if(this->solve_type == "exact_matvec"){
    exact_matvec(target, source_weight, source, ans);
  }else if(this->solve_type == "nrnmm"){
    nrnmm(target, source_weight, source, ans);
  }else if(this->solve_type == "tree"){
    tree(target, source_weight, source, ans);
  }else if(this->solve_type == "fmm"){
    fmm(target, source_weight, source, ans);
  }else{
    fprintf(stderr, "%s:%d ERROR : solve type %s not undefined\n", __FILE__, __LINE__, this->solve_type.c_str());
    exit(EXIT_FAILURE);
  }
  end = std::chrono::system_clock::now();

  time_log["total time(" + this->solve_type + ")"] = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

  std::multimap<double, std::pair<std::string, double>, std::greater<double>> mmap;
  for(auto itr : time_log){
    mmap.insert(std::make_pair(itr.second, itr));
  }
  FILE *fp;
  fp = fopen("fmma.log", "w");
  fprintf(fp, "name time[ms]\n");
  for(auto itr : mmap){
    fprintf(stderr, "%s %lf\n", itr.second.first.c_str(), itr.second.second);
    fprintf(fp, "%s %lf\n", itr.second.first.c_str(), itr.second.second);
  }
  fclose(fp);

  return;
};

template void FMMA<double, 1>::solve(const std::vector<std::array<double, 1>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::solve(const std::vector<std::array<double, 2>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);
template void FMMA<double, 3>::solve(const std::vector<std::array<double, 3>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 3>>& source, std::vector<double>& ans);

} // namespace fmma
