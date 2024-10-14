#pragma once

#include<cstdlib>
#include<cstdio>
#include<vector>
#include<array>
#include<functional>
#include<cmath>

namespace fmma {

template<typename TYPE, std::size_t DIM>
class FMMA{
  public:
    std::function<TYPE(std::array<TYPE, DIM>& source, std::array<TYPE, DIM>& target)> fn = 
      [](std::array<TYPE, DIM>& source, std::array<TYPE, DIM>& target){
        double len = 0.0;
        for(std::size_t dim=0; dim<DIM; ++dim){
          double diff = source[dim]-target[dim];
          len += diff*diff;
        }
        return 1.0/std::sqrt(len);
      };

    std::string solve_type = "exact";

  public:
    FMMA(void);
    ~FMMA(void);
    void solve(std::vector<std::array<TYPE, DIM>>& source, std::vector<std::array<TYPE, DIM>>& target, std::vector<TYPE>& ans);
    void exact(std::vector<std::array<TYPE, DIM>>& source, std::vector<std::array<TYPE, DIM>>& target, std::vector<TYPE>& ans);
};

} // namespace fmma
