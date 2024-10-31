#include"../include/fmma/math.hpp"
#include"../include/fmma/fmma.hpp"
#include<cstdlib>
#include<cstdio>
#include<vector>
#include<array>
#include<functional>
#include<cmath>
#if FMMA_USE_OPENMP
#include<omp.h>
#endif

namespace fmma {

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::fmm(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){
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

  fprintf(stderr, "fmm poly order = %d\n", poly_ord);
  fprintf(stderr, "fmm Depth = %d\n", Depth);
#if FMMA_USE_OPENMP
  fprintf(stderr, "omp parallel num = %d\n", omp_get_max_threads());
#endif

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
    // depth=0, 1のWmは使わない
    start = std::chrono::system_clock::now();
    for(int depth=0; depth+3<Depth; ++depth){
      M2M(tmp_N, chebyshev_node_all, Wm[Depth-depth-1], Wm[Depth-depth-2]);
      tmp_N /= 2;
    }
    end = std::chrono::system_clock::now();
    time_log["M2M"] += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
  }

  std::vector<std::vector<std::vector<TYPE>>> Wl(Depth);
    // depth=0, 1のWlは使わない
  for(int depth=2; depth<Depth; ++depth){
    Wl[depth].resize(Wm[depth].size());
    for(std::size_t ind=0; ind<Wl[depth].size(); ++ind){
      Wl[depth][ind].resize(chebyshev_node_all.size());
      for(std::size_t k=0; k<chebyshev_node_all.size(); ++k){
        Wl[depth][ind][k] = 0.0;
      }
    }
  }

  {
    // M2L
    std::size_t tmp_N = 1;
    start = std::chrono::system_clock::now();
    for(int depth=0; depth<Depth; ++depth){
      if(depth>=2){
        // depth=0, 1のWmは使わない
        M2L(tmp_N, origin, Len, chebyshev_node_all, Wm[depth], Wl[depth]);
      }
      tmp_N *= 2;
    }
    end = std::chrono::system_clock::now();
    time_log["M2L"] += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    // L2L
    tmp_N = 1;
    start = std::chrono::system_clock::now();
    for(int depth=0; depth+1<Depth; ++depth){
      if(depth>=2){
        // depth=0, 1のWlは使わない
        L2L(tmp_N*2, chebyshev_node_all, Wl[depth], Wl[depth+1]);
      }
      tmp_N *= 2;
    }
    end = std::chrono::system_clock::now();
    time_log["L2L"] += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
  }

  start = std::chrono::system_clock::now();
  if(Depth-1>=2){
    // depth=0, 1のWlは使わない
    L2P(target, origin, Len, chebyshev_node_all, Wl[Depth-1], ans);
  }
  end = std::chrono::system_clock::now();
  time_log["L2P"] += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();


  start = std::chrono::system_clock::now();
  // P2P
#pragma omp parallel for
  for(std::size_t t=0; t<target.size(); ++t){
    std::size_t tmp_N = 1<<(Depth-1);
    std::array<int, DIM> target_ind_of_box;
    for(std::size_t dim=0; dim<DIM; ++dim){
      TYPE len = Len/tmp_N;
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

template void FMMA<double, 1>::fmm(const std::vector<std::array<double, 1>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::fmm(const std::vector<std::array<double, 2>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);
template void FMMA<double, 3>::fmm(const std::vector<std::array<double, 3>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 3>>& source, std::vector<double>& ans);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::fmm(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){
  std::vector<TYPE> ones(source.size());
  for(std::size_t i=0; i<source.size(); ++i){
    ones[i] = 1.0;
  }

  fmm(target, ones, source, ans);

  return;
};

template void FMMA<double, 1>::fmm(const std::vector<std::array<double, 1>>& target, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::fmm(const std::vector<std::array<double, 2>>& target, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);
template void FMMA<double, 3>::fmm(const std::vector<std::array<double, 3>>& target, const std::vector<std::array<double, 3>>& source, std::vector<double>& ans);

}; // namespace fmma
