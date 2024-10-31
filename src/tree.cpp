#include"../include/fmma/math.hpp"
#include"../include/fmma/fmma.hpp"
#include<cstdlib>
#include<cstdio>
#include<vector>
#include<array>
#include<functional>
#include<cmath>

namespace fmma {

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::tree(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){
  std::chrono::system_clock::time_point start, end;

  {
    std::size_t M = source.size();
    std::size_t N = target.size();
    ans.resize(N);
    if(Depth <= 0){
      Depth = (int)(log(N)/DIM)+1;
      if(Depth <= 0){
        Depth = 1;
      }
    }
    if(M == 0 || N == 0){
      return;
    }
    if(source_weight.size() != M){
      fprintf(stderr, "%s:%d ERROR : source_weight size (%zu) != source size (%zu)\n", __FILE__, __LINE__, source_weight.size(), source.size());
      exit(EXIT_FAILURE);
    }
  }

  std::array<TYPE, DIM> origin;
  TYPE Len;
  get_origin_and_length(target, source, origin, Len);

  std::vector<std::vector<std::vector<TYPE>>> Wm(Depth);
  std::vector<std::vector<std::size_t>> source_ind_in_box;
  std::vector<std::array<TYPE, DIM>> chebyshev_node_all;

  {
    // P2M
    std::size_t tmp_N = 1<<(Depth-1);
    start = std::chrono::system_clock::now();
    P2M(source_weight, source, tmp_N, origin, Len/tmp_N, source_ind_in_box, Wm[Depth-1], chebyshev_node_all);
    end = std::chrono::system_clock::now();
    time_log["P2M"] += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    // M2M
    for(int depth=0; depth+1<Depth; ++depth){
      start = std::chrono::system_clock::now();
      M2M(tmp_N, chebyshev_node_all, Wm[Depth-depth-1], Wm[Depth-depth-2]);
      end = std::chrono::system_clock::now();
      time_log["M2M"] += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
      tmp_N /= 2;
    }
  }

  // M2P
  start = std::chrono::system_clock::now();
  for(std::size_t t=0; t<target.size(); ++t){
    ans[t] = 0.0;
  }
  int tmp_N = 1;
  for(int depth=0; depth<Depth; ++depth){
    M2P(target, tmp_N, origin, Len, chebyshev_node_all, Wm[depth], ans);
    tmp_N *= 2;
  }
  end = std::chrono::system_clock::now();
  time_log["M2P"] += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

  start = std::chrono::system_clock::now();
  // P2P
#pragma omp parallel for
  for(std::size_t t=0; t<target.size(); ++t){
    std::size_t tmp_N = 1<<(Depth-1);
    TYPE len = Len/tmp_N;

    std::array<int, DIM> target_ind_of_box;
    for(std::size_t dim=0; dim<DIM; ++dim){
      target_ind_of_box[dim] = std::min((int)((target[t][dim]-origin[dim])/len), (int)tmp_N-1);
    }
    std::vector<std::size_t> indices = exact_calc_box_indices(target_ind_of_box, tmp_N);
    for(std::size_t i=0; i<indices.size(); ++i){
      for(std::size_t j=0; j<source_ind_in_box[indices[i]].size(); ++j){
        ans[t] += source_weight[source_ind_in_box[indices[i]][j]]*fn(target[t], source[source_ind_in_box[indices[i]][j]]);
      }
    }
  }
  end = std::chrono::system_clock::now();
  time_log["P2P"] += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

  return;
};

template void FMMA<double, 1>::tree(const std::vector<std::array<double, 1>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::tree(const std::vector<std::array<double, 2>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);
template void FMMA<double, 3>::tree(const std::vector<std::array<double, 3>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 3>>& source, std::vector<double>& ans);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::tree(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){
  std::vector<TYPE> ones(source.size());
  for(std::size_t i=0; i<source.size(); ++i){
    ones[i] = 1.0;
  }

  tree(target, ones, source, ans);

  return;
};

template void FMMA<double, 1>::tree(const std::vector<std::array<double, 1>>& target, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::tree(const std::vector<std::array<double, 2>>& target, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);
template void FMMA<double, 3>::tree(const std::vector<std::array<double, 3>>& target, const std::vector<std::array<double, 3>>& source, std::vector<double>& ans);

}; // namespace fmma
