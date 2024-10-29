#include<cstdlib>
#include<cstdio>
#include<vector>
#include<array>
#include<functional>
#include<cmath>
#include"../include/fmma/fmma.hpp"
#include"../include/fmma/math.hpp"

namespace fmma {

 
template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::nrnmm(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){

  std::size_t M = source.size();
  std::size_t N = target.size();
  ans.resize(N);
  if(nrn_N <= 0){
    nrn_N = (int)(sqrt(3.0)*pow((TYPE)M, 1.0/(2.0*DIM)))+1.0;
  }
  if(M == 0 || N == 0){
    return;
  }
  if(source_weight.size() != M){
    fprintf(stderr, "%s:%d ERROR : source_weight size (%zu) != source size (%zu)\n", __FILE__, __LINE__, source_weight.size(), source.size());
    exit(EXIT_FAILURE);
  }

  std::size_t SIZE = 1;
  for(std::size_t i=0; i<DIM; ++i){
    SIZE *= nrn_N;
  }

  std::array<TYPE, DIM> origin;
  TYPE Len;
  get_origin_and_length(target, source, origin, Len);
  TYPE len = Len/nrn_N;

  std::vector<std::vector<TYPE>> Wm;
  std::vector<std::vector<std::size_t>> source_ind_in_box;
  std::vector<std::array<TYPE, DIM>> chebyshev_node_all;

  P2M(source_weight, source, nrn_N, origin, len, source_ind_in_box, Wm, chebyshev_node_all);

  std::size_t poly_ord_all = chebyshev_node_all.size();

  std::array<TYPE, DIM> chebyshev_real_pos;
  std::array<TYPE, DIM> relative_orig_pos;
  std::array<int, DIM> target_ind_of_box;
  std::array<int, DIM> ind_of_box;
  for(std::size_t t=0; t<target.size(); ++t){
    ans[t] = 0.0;

    for(std::size_t dim=0; dim<DIM; ++dim){
      target_ind_of_box[dim] = std::min((int)((target[t][dim]-origin[dim])/len), nrn_N-1);
    }

    for(std::size_t s=0; s<SIZE; ++s){
      std::size_t s_copy = s;
      int max_dist_from_t = 0;
      for(std::size_t dim=0; dim<DIM; ++dim){
        relative_orig_pos[DIM-1-dim] = len*(s_copy%nrn_N)+origin[DIM-1-dim];
        ind_of_box[DIM-1-dim] = s_copy%nrn_N;
        s_copy /= nrn_N;

        max_dist_from_t = std::max(max_dist_from_t, std::abs(ind_of_box[DIM-1-dim]-target_ind_of_box[DIM-1-dim]));
      }

      if(max_dist_from_t <= 1){
        for(std::size_t i=0; i<source_ind_in_box[s].size(); ++i){
          ans[t] += source_weight[source_ind_in_box[s][i]]*fn(target[t]-source[source_ind_in_box[s][i]]);
        }
      }else{
        for(std::size_t k=0; k<poly_ord_all; ++k){
          for(std::size_t dim=0; dim<DIM; ++dim){
            chebyshev_real_pos[dim] = (chebyshev_node_all[k][dim]+1.0)/2.0*len+relative_orig_pos[dim];
          }
          ans[t] += fn(target[t]-chebyshev_real_pos)*Wm[s][k];
        }
      }
    }
  }

  return;
};

template void FMMA<double, 1>::nrnmm(const std::vector<std::array<double, 1>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::nrnmm(const std::vector<std::array<double, 2>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::nrnmm(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){
  std::vector<TYPE> ones(source.size());
  for(std::size_t i=0; i<source.size(); ++i){
    ones[i] = 1.0;
  }

  nrnmm(target, ones, source, ans);

  return;
};

template void FMMA<double, 1>::nrnmm(const std::vector<std::array<double, 1>>& target, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::nrnmm(const std::vector<std::array<double, 2>>& target, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);

}; // namespace fmma
