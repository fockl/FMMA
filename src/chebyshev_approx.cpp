#include<cstdlib>
#include<cstdio>
#include<vector>
#include<array>
#include<functional>
#include<cmath>
#include"../include/fmma/math.hpp"
#include"../include/fmma/chebyshev_approx.hpp"

namespace fmma {

template<typename TYPE, std::size_t DIM>
CHEBYSHEV_APPROX<TYPE, DIM>::CHEBYSHEV_APPROX(void){
  return;
};

template CHEBYSHEV_APPROX<double, 1>::CHEBYSHEV_APPROX(void);
template CHEBYSHEV_APPROX<double, 2>::CHEBYSHEV_APPROX(void);
template CHEBYSHEV_APPROX<double, 3>::CHEBYSHEV_APPROX(void);

template<typename TYPE, std::size_t DIM>
CHEBYSHEV_APPROX<TYPE, DIM>::~CHEBYSHEV_APPROX(void){
  return;
};

template CHEBYSHEV_APPROX<double, 1>::~CHEBYSHEV_APPROX(void);
template CHEBYSHEV_APPROX<double, 2>::~CHEBYSHEV_APPROX(void);
template CHEBYSHEV_APPROX<double, 3>::~CHEBYSHEV_APPROX(void);

template<typename TYPE, std::size_t DIM>
void CHEBYSHEV_APPROX<TYPE, DIM>::initialize(void){
  size=1;
  for(std::size_t i=0; i<DIM; ++i){
    size *= (Nth+1);
  }
  coef.resize(size);

  std::vector<TYPE> chebyshev_node(Nth+1);
  for(int k=0; k<=Nth; ++k){
    chebyshev_node[k] = cos((2.0*k+1.0)/(2*Nth+2)*M_PI);
  }

  std::array<TYPE, DIM> pos;
  std::vector<TYPE> Mat(size*size);
  std::vector<TYPE> b(size);
  for(std::size_t i=0; i<size; ++i){
    std::size_t i_copy = i;
    for(std::size_t dim1=0; dim1<DIM; ++dim1){
      pos[dim1] = chebyshev_node[i_copy%(Nth+1)];
      i_copy /= (Nth+1);
    }
    b[i] = fn(pos);

    for(std::size_t j=0; j<size; ++j){
      std::size_t j_copy = j;
      TYPE val = 1.0;
      for(std::size_t dim2=0; dim2<DIM; ++dim2){
        val *= Chebyshev(j_copy%(Nth+1), pos[dim2]);
        j_copy /= (Nth+1);
      }
      Mat[i*size+j] = val;
    }
  }

  solve(Mat, coef, b);

  return;
}

template void CHEBYSHEV_APPROX<double, 1>::initialize(void);
template void CHEBYSHEV_APPROX<double, 2>::initialize(void);
template void CHEBYSHEV_APPROX<double, 3>::initialize(void);

template<typename TYPE, std::size_t DIM>
TYPE CHEBYSHEV_APPROX<TYPE, DIM>::predict(const std::array<TYPE, DIM>& pos){
  TYPE ans = 0.0;
  for(std::size_t i=0; i<size; ++i){
    std::size_t i_copy = i;
    TYPE val = 1.0;
    for(std::size_t dim1=0; dim1<DIM; ++dim1){
      val *= Chebyshev(i_copy%(Nth+1), pos[dim1]);
      i_copy /= (Nth+1);
    }
    ans += coef[i]*val;
  }

  return ans;
}

template double CHEBYSHEV_APPROX<double, 1>::predict(const std::array<double, 1>& pos);
template double CHEBYSHEV_APPROX<double, 2>::predict(const std::array<double, 2>& pos);
template double CHEBYSHEV_APPROX<double, 3>::predict(const std::array<double, 3>& pos);

}; // namespace fmma

