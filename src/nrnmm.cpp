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

  std::size_t SIZE = 1;
  for(std::size_t i=0; i<DIM; ++i){
    SIZE *= nrn_N;
  }

  std::array<TYPE, DIM> min_pos, max_pos;
  min_pos = source[0];
  max_pos = source[0];
  for(std::size_t s=0; s<source.size(); ++s){
    for(std::size_t dim=0; dim<DIM; ++dim){
      min_pos[dim] = std::min(min_pos[dim], source[s][dim]);
      max_pos[dim] = std::max(max_pos[dim], source[s][dim]);
    }
  }
  for(std::size_t t=0; t<target.size(); ++t){
    for(std::size_t dim=0; dim<DIM; ++dim){
      min_pos[dim] = std::min(min_pos[dim], target[t][dim]);
      max_pos[dim] = std::max(max_pos[dim], target[t][dim]);
    }
  }
  TYPE Len = 0.0;
  for(std::size_t dim=0; dim<DIM; ++dim){
    Len = std::max(Len, max_pos[dim] - min_pos[dim]);
  }
  TYPE len = Len/nrn_N;

  for(std::size_t dim=0; dim<DIM; ++dim){
    min_pos[dim] = (max_pos[dim] + min_pos[dim])/2 - Len/2;
    max_pos[dim] = min_pos[dim] + Len;
  }

  std::size_t poly_ord_all = 1;
  for(std::size_t dim=0; dim<DIM; ++dim){
    poly_ord_all *= (poly_ord+1);
  }
  std::vector<std::vector<TYPE>> Wm(SIZE, std::vector<TYPE>(poly_ord_all));
  std::vector<std::vector<std::size_t>> source_ind_in_box(SIZE);
  std::vector<TYPE> chebyshev_node(poly_ord+1);
  for(int k=0; k<=poly_ord; ++k){
    chebyshev_node[k] = cos((2.0*k+1.0)/(2*poly_ord+2)*M_PI);
  }
  std::vector<std::array<TYPE, DIM>> chebyshev_node_all(poly_ord_all);
  for(std::size_t k=0; k<poly_ord_all; ++k){
    std::size_t k_copy = k;
    for(std::size_t dim=0; dim<DIM; ++dim){
      std::size_t ord = k_copy%(poly_ord+1);
      chebyshev_node_all[k][DIM-1-dim] = chebyshev_node[ord];
      k_copy /= (poly_ord+1);
    }
  }

  std::array<TYPE, DIM> relative_pos;
  for(std::size_t s=0; s<source.size(); ++s){
    std::array<TYPE, DIM> source_pos = source[s];
    std::size_t pos_node = 0;
    for(std::size_t dim=0; dim<DIM; ++dim){
      int tmp = std::min((int)((source_pos[dim]-min_pos[dim])/len), nrn_N-1);
      pos_node *= nrn_N;
      pos_node += tmp;
      relative_pos[dim] = std::max(std::min(2.0*((source_pos[dim]-min_pos[dim])/len-tmp)-1.0, 1.0), -1.0);
    }

    source_ind_in_box[pos_node].push_back(s);

    for(std::size_t k=0; k<poly_ord_all; ++k){
      TYPE val = 1.0;
      for(std::size_t dim=0; dim<DIM; ++dim){
        val *= SChebyshev(poly_ord+1, chebyshev_node_all[k][dim], relative_pos[dim]);
      }
      Wm[pos_node][k] += source_weight[s]*val;
    }
  }

  std::array<TYPE, DIM> chebyshev_real_pos;
  std::array<TYPE, DIM> relative_orig_pos;
  std::array<int, DIM> target_ind_of_box;
  std::array<int, DIM> ind_of_box;
  for(std::size_t t=0; t<target.size(); ++t){
    ans[t] = 0.0;

    for(std::size_t dim=0; dim<DIM; ++dim){
      target_ind_of_box[DIM-1-dim] = std::min((int)((target[t][dim]-min_pos[dim])/len), nrn_N-1);
    }

    for(std::size_t s=0; s<SIZE; ++s){
      std::size_t s_copy = s;
      int max_dist_from_t = 0;
      for(std::size_t dim=0; dim<DIM; ++dim){
        relative_orig_pos[DIM-1-dim] = len*(s_copy%nrn_N)+min_pos[DIM-1-dim];
        ind_of_box[dim] = s_copy%nrn_N;
        s_copy /= nrn_N;

        max_dist_from_t = std::max(max_dist_from_t, std::abs(ind_of_box[dim]-target_ind_of_box[dim]));
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

};
