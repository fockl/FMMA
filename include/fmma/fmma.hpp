#pragma once

#include<cstdlib>
#include<cstdio>
#include<vector>
#include<array>
#include<functional>
#include<cmath>
#include<string>
#include<chrono>
#include<map>

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
    std::function<TYPE(const std::array<TYPE, DIM>& target, const std::array<TYPE, DIM>& source)> fn = 
      [](const std::array<TYPE, DIM>& target, const std::array<TYPE, DIM>& source){
        double len = 0.0;
        for(std::size_t dim=0; dim<DIM; ++dim){
          double diff = target[dim]-source[dim];
          len += diff*diff;
        }
        return 1.0/std::sqrt(len);
      };

    void set_fn(const std::function<TYPE(const std::array<TYPE, DIM>& target_source)>& fn){
      // fn = f(x) と1変数関数で与えた場合、
      // f(x, y) = f(x-y)と解釈する
      this->fn = [fn](const std::array<TYPE, DIM>& target, const std::array<TYPE, DIM>& source){
        return fn(target-source);
      };
      return;
    }

    void set_fn1(const std::function<TYPE(const std::array<TYPE, DIM>& target_source)>& fn){
      set_fn(fn);
      return;
    }

    void set_fn(const std::function<TYPE(const std::array<TYPE, DIM>& target, const std::array<TYPE, DIM>& source)>& fn){
      this->fn = fn;
      return;
    }
    
    void set_fn2(const std::function<TYPE(const std::array<TYPE, DIM>& target, const std::array<TYPE, DIM>& source)>& fn){
      set_fn(fn);
      return;
    }

  private:
    std::string solver_type = "exact";
  public:
    void set_solver_type(const std::string new_type){
      solver_type = new_type;
      return;
    }

  private:
    int nrn_N = -1; // 1辺辺りの分割数 nrnで指定
  public:
    void set_nrn_N(int new_nrn_N){
      nrn_N = new_nrn_N;
      return;
    }

  private:
    int poly_ord = 1;
  public:
    void set_poly_ord(int new_poly_ord){
      poly_ord = new_poly_ord;
      return;
    }

  private:
    int Depth = -1; // 深さ、2^Depth = 最深部の一辺当たりの分割数
  public:
    void set_Depth(int new_Depth){
      Depth = new_Depth;
      return;
    }

  private:
    int trans_sym_flag = 0; // 並進対称性flag, 0: 対称性なし、1: f(x, y) = f(x-y)
  public:
    void set_trans_sym_flag(int new_flag){
      trans_sym_flag = new_flag;
      return;
    }

  public:
    FMMA(void);
    ~FMMA(void);

    bool check_blas(void);

    void solve(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void solve(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);

    void exact(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void exact(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void exact_matvec(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void exact_matvec(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void nrnmm(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void nrnmm(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void tree(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void tree(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void fmm(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);
    void fmm(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans);


  private:
    void get_minmax(const std::vector<std::array<TYPE, DIM>>& array1, const std::vector<std::array<TYPE, DIM>>& array2, std::array<TYPE, DIM>& min_array, std::array<TYPE, DIM>& max_array);
  protected:
    void get_origin_and_length(const std::vector<std::array<TYPE, DIM>>& array1, const std::vector<std::array<TYPE, DIM>>& array2, std::array<TYPE, DIM>& origin, TYPE& length);
    void P2M(const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, int N, const std::array<TYPE, DIM>& min_pos, const TYPE len, std::vector<std::vector<std::size_t>>& source_ind_in_box, std::vector<std::vector<TYPE>>& Wm, std::vector<std::array<TYPE, DIM>>& chebyshev_node_all);
    void M2M(const std::size_t N, const std::vector<std::array<TYPE, DIM>>& chebyshev_node_all, const std::vector<std::vector<TYPE>>& Wm_in, std::vector<std::vector<TYPE>>& Wm_out);
    void M2P(const std::vector<std::array<TYPE, DIM>>& target, const std::size_t N, const std::array<TYPE, DIM>& origin, const TYPE Len, const std::vector<std::array<TYPE, DIM>>& chebyshev_node_all, const std::vector<std::vector<TYPE>>& Wm, std::vector<TYPE>& ans);
    void M2L(const std::size_t N, const std::array<TYPE, DIM>& origin, const TYPE Len, const std::vector<std::array<TYPE, DIM>>& chebyshev_node_all, const std::vector<std::vector<TYPE>>& Wm, std::vector<std::vector<TYPE>>& Wl);
    void L2L(const std::size_t N, const std::vector<std::array<TYPE, DIM>>& chebyshev_node_all, const std::vector<std::vector<TYPE>>& Wl_in, std::vector<std::vector<TYPE>>& Wl_out);
    void L2P(const std::vector<std::array<TYPE, DIM>>& target, const std::array<TYPE, DIM>& origin, const TYPE Len, const std::vector<std::array<TYPE, DIM>>& chebyshev_node_all, const std::vector<std::vector<TYPE>>& Wl, std::vector<TYPE>& ans);
  private:
    std::array<std::size_t, DIM> get_box_ind_of_ind(const std::size_t ind, int N);
    std::size_t get_ind_of_box_ind(const std::array<int, DIM>& box_ind, int N);
    template<typename INT>
    std::vector<std::size_t> multipole_calc_box_indices(const std::array<INT, DIM>& box_ind, int N);
    std::vector<std::size_t> exact_calc_box_indices(const std::array<int, DIM>& box_ind, int N);

  private:
    std::map<std::string, double> time_log;

};

} // namespace fmma
