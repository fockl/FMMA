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
    std::function<TYPE(const std::array<TYPE, DIM>& source, const std::array<TYPE, DIM>& target)> fn = 
      [](const std::array<TYPE, DIM>& source, const std::array<TYPE, DIM>& target){
        double len = 0.0;
        for(std::size_t dim=0; dim<DIM; ++dim){
          double diff = source[dim]-target[dim];
          len += diff*diff;
        }
        return 1.0/std::sqrt(len);
      };

    std::string solve_type = "exact";

    int nrn_N = -1; // 1辺辺りの分割数 nrnで使用
    int poly_ord = 1;

  public:
    FMMA(void);
    ~FMMA(void);
    void solve(const std::vector<std::array<TYPE, DIM>>& source, const std::vector<std::array<TYPE, DIM>>& target, std::vector<TYPE>& ans);
    void exact(const std::vector<std::array<TYPE, DIM>>& source, const std::vector<std::array<TYPE, DIM>>& target, std::vector<TYPE>& ans);
    void nrnmm(const std::vector<std::array<TYPE, DIM>>& source, const std::vector<std::array<TYPE, DIM>>& target, std::vector<TYPE>& ans);
};

} // namespace fmma
