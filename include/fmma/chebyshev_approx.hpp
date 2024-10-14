#pragma once

#include<cstdlib>
#include<cstdio>
#include<vector>
#include<array>
#include<functional>
#include<cmath>

namespace fmma {

template<typename TYPE, std::size_t DIM>
class CHEBYSHEV_APPROX{
  public:
    int Nth = 0; // 近似次数
    std::function<TYPE(const std::array<TYPE, DIM>& pos)> fn = 
      [](const std::array<TYPE, DIM>& pos){
        double len2 = 0.0;
        for(std::size_t dim=0; dim<DIM; ++dim){
          len2 += pos[dim]*pos[dim];
        }
        return len2;
      };

  public:
    std::vector<TYPE> coef; // 係数
    std::size_t size; // 要素数

  public:
    CHEBYSHEV_APPROX(void);
    ~CHEBYSHEV_APPROX(void);
    void initialize(void);
    TYPE predict(const std::array<TYPE, DIM>& pos);

};


}; // namespace fmma

