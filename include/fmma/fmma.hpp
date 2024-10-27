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

    std::string solve_type = "exact";

    int nrn_N = -1; // 1辺辺りの分割数 nrnで指定
    int poly_ord = 1;
    int Depth = -1; // 深さ、2^Depth = 最深部の一辺当たりの分割数

  public:
    FMMA(void);
    ~FMMA(void);
    void solve(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void solve(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void exact(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void exact(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void nrnmm(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void nrnmm(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void tree(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void tree(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);

  private:
    void get_minmax(const std::vector<std::array<TYPE, DIM>>& array1, const std::vector<std::array<TYPE, DIM>>& array2, std::array<TYPE, DIM>& min_array, std::array<TYPE, DIM>& max_array);
  protected:
    void get_origin_and_length(const std::vector<std::array<TYPE, DIM>>& array1, const std::vector<std::array<TYPE, DIM>>& array2, std::array<TYPE, DIM>& origin, TYPE& length);
    void P2M(const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, int N, const std::array<TYPE, DIM>& min_pos, const TYPE len, std::vector<std::vector<std::size_t>>& source_ind_in_box, std::vector<std::vector<TYPE>>& Wm, std::vector<std::array<TYPE, DIM>>& chebyshev_node_all);
    void M2M(const std::size_t N, const std::vector<std::array<TYPE, DIM>>& chebyshev_node_all, const std::vector<std::vector<TYPE>>& Wm_in, std::vector<std::vector<TYPE>>& Wm_out);
  private:
    std::array<std::size_t, DIM> get_box_ind_of_ind(const std::size_t ind, int N);
    std::size_t get_ind_of_box_ind(const std::array<int, DIM>& box_ind, int N);
    std::vector<std::size_t> multipole_calc_box_indices(const std::array<int, DIM>& box_ind, int N);
    std::vector<std::size_t> exact_calc_box_indices(const std::array<int, DIM>& box_ind, int N);

};

} // namespace fmma
