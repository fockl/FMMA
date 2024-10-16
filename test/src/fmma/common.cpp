#include"../../../include/fmma/fmma.hpp"
#include"../test_common.hpp"
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<vector>
#include<array>

//{{{ test_exact
template<typename TYPE, std::size_t DIM>
bool test_exact(int ssize, int tsize, TYPE tol){
  std::vector<std::array<TYPE, DIM>> source(ssize), target(tsize);
  for(int i=0; i<ssize; ++i){
    for(std::size_t dim=0; dim<DIM; ++dim){
      source[i][dim] = (TYPE)rand()/RAND_MAX-0.5;
    }
  }
  for(int i=0; i<tsize; ++i){
    for(std::size_t dim=0; dim<DIM; ++dim){
      target[i][dim] = (TYPE)rand()/RAND_MAX-0.5;
    }
  }

  fmma::FMMA<TYPE, DIM> fmma;
  std::vector<TYPE> ans(tsize);
  fmma.exact(source, target, ans);

  std::vector<TYPE> exact(tsize);
  for(int i=0; i<ssize; ++i){
    for(int j=0; j<tsize; ++j){
      TYPE len = 0.0;
      for(std::size_t dim=0; dim<DIM; ++dim){
        TYPE d = source[i][dim] - target[j][dim];
        len += d * d;
      }
      len = std::sqrt(len);
      exact[j] += 1.0/len;
    }
  }

  TYPE diff = 0.0;
  for(int i=0; i<tsize; ++i){
    TYPE d = ans[i]-exact[i];
    diff += d*d;
  }
  diff = sqrt(diff / tsize);

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

//{{{ test_exact_1_r2
template<typename TYPE, std::size_t DIM>
bool test_exact_1_r2(int ssize, int tsize, TYPE tol){
  std::vector<std::array<TYPE, DIM>> source(ssize), target(tsize);
  for(int i=0; i<ssize; ++i){
    for(std::size_t dim=0; dim<DIM; ++dim){
      source[i][dim] = (TYPE)rand()/RAND_MAX-0.5;
    }
  }
  for(int i=0; i<tsize; ++i){
    for(std::size_t dim=0; dim<DIM; ++dim){
      target[i][dim] = (TYPE)rand()/RAND_MAX-0.5;
    }
  }

  fmma::FMMA<TYPE, DIM> fmma;
  std::vector<TYPE> ans(tsize);
  auto fn = [](const std::array<TYPE, DIM>& source, const std::array<TYPE, DIM>& target){
    TYPE len2 = 0.0;
    for(std::size_t dim=0; dim<DIM; ++dim){
      TYPE diff = source[dim]-target[dim];
      len2 += diff*diff;
    }
    return 1.0/len2;
  };
  fmma.fn = fn;
  fmma.exact(source, target, ans);

  std::vector<TYPE> exact(tsize);
  for(int i=0; i<ssize; ++i){
    for(int j=0; j<tsize; ++j){
      TYPE len = 0.0;
      for(std::size_t dim=0; dim<DIM; ++dim){
        TYPE d = source[i][dim] - target[j][dim];
        len += d * d;
      }
      exact[j] += 1.0/len;
    }
  }

  TYPE diff = 0.0;
  for(int i=0; i<tsize; ++i){
    TYPE d = ans[i]-exact[i];
    diff += d*d;
  }
  diff = sqrt(diff / tsize);

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

  if(!test_exact<double, 1>(10, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_exact<double, 1>(10, 20, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_exact<double, 1>(20, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_exact<double, 2>(10, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_exact<double, 2>(10, 20, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_exact<double, 2>(20, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_exact_1_r2<double, 1>(10, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_exact_1_r2<double, 1>(10, 20, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_exact_1_r2<double, 1>(20, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_exact_1_r2<double, 2>(10, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_exact_1_r2<double, 2>(10, 20, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_exact_1_r2<double, 2>(20, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  return 0;
}
