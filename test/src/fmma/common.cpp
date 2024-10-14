#include"../../../include/fmma.hpp"
#include"../test_common.hpp"
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<vector>
#include<array>

template<typename TYPE, std::size_t DIM>
bool test_solve(int ssize, int tsize, TYPE tol){
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
  fmma.solve(source, target, ans);

  std::vector<TYPE> exact(tsize);
  for(int i=0; i<ssize; ++i){
    for(int j=0; j<tsize; ++j){
      double len = 0.0;
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
    fprintf(stderr, "diff = %e > rol = %e\n", diff, tol);
    failed(__FILE__, __func__);
    return false;
  }
}

int main(void){
  srand(0);

  if(!test_solve<double, 1>(10, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  return 0;
}
