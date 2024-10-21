#pragma once

#include<cstdlib>
#include<cstdio>
#include<vector>
#include<array>
#include<functional>
#include<cmath>
#include<string>

namespace fmma {

template<typename TYPE, std::size_t DIM>
std::array<TYPE, DIM> operator+(const std::array<TYPE, DIM>& lhs, const std::array<TYPE, DIM>& rhs){
  std::array<TYPE, DIM> ans;
  for(std::size_t dim=0; dim<DIM; ++dim){
    ans[dim] = lhs[dim]+rhs[dim];
  }
  return ans;
};

template<typename TYPE, std::size_t DIM>
std::array<TYPE, DIM> operator-(const std::array<TYPE, DIM>& lhs, const std::array<TYPE, DIM>& rhs){
  std::array<TYPE, DIM> ans;
  for(std::size_t dim=0; dim<DIM; ++dim){
    ans[dim] = lhs[dim]-rhs[dim];
  }
  return ans;
};

template<typename TYPE, std::size_t DIM>
std::array<TYPE, DIM> operator*(const std::array<TYPE, DIM>& lhs, const std::array<TYPE, DIM>& rhs){
  std::array<TYPE, DIM> ans;
  for(std::size_t dim=0; dim<DIM; ++dim){
    ans[dim] = lhs[dim]*rhs[dim];
  }
  return ans;
};

template<typename TYPE, std::size_t DIM>
std::array<TYPE, DIM> operator/(const std::array<TYPE, DIM>& lhs, const std::array<TYPE, DIM>& rhs){
  std::array<TYPE, DIM> ans;
  for(std::size_t dim=0; dim<DIM; ++dim){
    ans[dim] = lhs[dim]/rhs[dim];
  }
  return ans;
};


template<typename TYPE, std::size_t DIM>
class FMMA{
  public:
    std::function<TYPE(const std::array<TYPE, DIM>& target_source)> fn = 
      [](const std::array<TYPE, DIM>& target_source){
        double len = 0.0;
        for(std::size_t dim=0; dim<DIM; ++dim){
          double diff = target_source[dim];
          len += diff*diff;
        }
        return 1.0/std::sqrt(len);
      };

    void set_fn(const std::function<TYPE(const std::array<TYPE, DIM>& target_source)>& fn){
      this->fn = fn;
      return;
    }

    void set_fn(const std::function<TYPE(const std::array<TYPE, DIM>& target, const std::array<TYPE, DIM>& source)>& fn){
      this->fn = [fn](const std::array<TYPE, DIM>& target_source){
        std::array<TYPE, DIM> zero;
        for(std::size_t dim=0; dim<DIM; ++dim){
          zero[dim] = 0.0;
        }
        return fn(target_source, zero);
      };
      return;
    }
    /*
    std::function<TYPE(const std::array<TYPE, DIM>& target, const std::array<TYPE, DIM>& source)> fn = 
      [](const std::array<TYPE, DIM>& target, const std::array<TYPE, DIM>& source){
        double len = 0.0;
        for(std::size_t dim=0; dim<DIM; ++dim){
          double diff = target[dim] - source[dim];
          len += diff*diff;
        }
        return 1.0/std::sqrt(len);
      };
      */

    std::string solve_type = "exact";

    int nrn_N = -1; // 1辺辺りの分割数 nrnで使用
    int poly_ord = 1;

  public:
    FMMA(void);
    ~FMMA(void);
    void solve(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void solve(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void exact(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void exact(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void nrnmm(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void nrnmm(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
};

} // namespace fmma
