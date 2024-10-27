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
  void FMMA<TYPE, DIM>::get_minmax(const std::vector<std::array<TYPE, DIM>>& array1, const std::vector<std::array<TYPE, DIM>>& array2, std::array<TYPE, DIM>& min_array, std::array<TYPE, DIM>& max_array){
    if(array1.size()>0){
      min_array = array1[0];
      max_array = array1[0];
    }else if(array2.size()>0){
      min_array = array2[0];
      max_array = array2[0];
    }
    for(std::size_t i=0; i<array1.size(); ++i){
      for(std::size_t dim=0; dim<DIM; ++dim){
        min_array[dim] = std::min(min_array[dim], array1[i][dim]);
        max_array[dim] = std::max(max_array[dim], array1[i][dim]);
      }
    }
    for(std::size_t i=0; i<array2.size(); ++i){
      for(std::size_t dim=0; dim<DIM; ++dim){
        min_array[dim] = std::min(min_array[dim], array2[i][dim]);
        max_array[dim] = std::max(max_array[dim], array2[i][dim]);
      }
    }
    return;
  }

template void FMMA<double, 1>::get_minmax(const std::vector<std::array<double, 1>>& array1, const std::vector<std::array<double, 1>>& array2, std::array<double, 1>& min_array, std::array<double, 1>& max_array);
template void FMMA<double, 2>::get_minmax(const std::vector<std::array<double, 2>>& array1, const std::vector<std::array<double, 2>>& array2, std::array<double, 2>& min_array, std::array<double, 2>& max_array);

template<typename TYPE, std::size_t DIM>
  void FMMA<TYPE, DIM>::get_origin_and_length(const std::vector<std::array<TYPE, DIM>>& array1, const std::vector<std::array<TYPE, DIM>>& array2, std::array<TYPE, DIM>& origin, TYPE& Len){
    std::array<TYPE, DIM> min_pos, max_pos;
    get_minmax(array1, array2, min_pos, max_pos);

    Len = 0.0;
    for(std::size_t dim=0; dim<DIM; ++dim){
      Len = std::max(Len, max_pos[dim] - min_pos[dim]);
    }

    for(std::size_t dim=0; dim<DIM; ++dim){
      origin[dim] = (max_pos[dim] + min_pos[dim])/2 - Len/2;
    }
    return;
  }

template void FMMA<double, 1>::get_origin_and_length(const std::vector<std::array<double, 1>>& array1, const std::vector<std::array<double, 1>>& array2, std::array<double, 1>& origin, double& Len);
template void FMMA<double, 2>::get_origin_and_length(const std::vector<std::array<double, 2>>& array1, const std::vector<std::array<double, 2>>& array2, std::array<double, 2>& origin, double& Len);

template<typename TYPE, std::size_t DIM>
  void FMMA<TYPE, DIM>::P2M(const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, int N, const std::array<TYPE, DIM>& min_pos, const TYPE len, std::vector<std::vector<std::size_t>>& source_ind_in_box, std::vector<std::vector<TYPE>>& Wm, std::vector<std::array<TYPE, DIM>>& chebyshev_node_all){
    std::size_t SIZE = 1;
    for(std::size_t dim=0; dim<DIM; ++dim){
      SIZE *= N;
    }

    std::size_t poly_ord_all = 1;
    for(std::size_t dim=0; dim<DIM; ++dim){
      poly_ord_all *= (poly_ord+1);
    }

    Wm.resize(SIZE);
    for(std::size_t i=0; i<SIZE; ++i){
      Wm[i].resize(poly_ord_all);
    }
    source_ind_in_box.resize(SIZE);
    std::vector<TYPE> chebyshev_node(poly_ord+1);
    for(int k=0; k<=poly_ord; ++k){
      chebyshev_node[k] = cos((2.0*k+1.0)/(2*poly_ord+2)*M_PI);
    }
    chebyshev_node_all.resize(poly_ord_all);
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
        int tmp = std::min((int)((source_pos[dim]-min_pos[dim])/len), N-1);
        pos_node *= N;
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

    return;
  }

template void FMMA<double, 1>::P2M(const std::vector<double>& source_weight, const std::vector<std::array<double, 1>>& source, int N, const std::array<double, 1>& min_pos, const double len, std::vector<std::vector<std::size_t>>& source_ind_in_box, std::vector<std::vector<double>>& Wm, std::vector<std::array<double, 1>>& chebyshev_node_all);
template void FMMA<double, 2>::P2M(const std::vector<double>& source_weight, const std::vector<std::array<double, 2>>& source, int N, const std::array<double, 2>& min_pos, const double len, std::vector<std::vector<std::size_t>>& source_ind_in_box, std::vector<std::vector<double>>& Wm, std::vector<std::array<double, 2>>& chebyshev_node_all);

  
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

  std::array<TYPE, DIM> min_pos, max_pos;
  get_minmax(target, source, min_pos, max_pos);

  TYPE Len = 0.0;
  for(std::size_t dim=0; dim<DIM; ++dim){
    Len = std::max(Len, max_pos[dim] - min_pos[dim]);
  }
  TYPE len = Len/nrn_N;

  for(std::size_t dim=0; dim<DIM; ++dim){
    min_pos[dim] = (max_pos[dim] + min_pos[dim])/2 - Len/2;
    max_pos[dim] = min_pos[dim] + Len;
  }

  std::vector<std::vector<TYPE>> Wm;
  std::vector<std::vector<std::size_t>> source_ind_in_box;
  std::vector<std::array<TYPE, DIM>> chebyshev_node_all;

  P2M(source_weight, source, nrn_N, min_pos, len, source_ind_in_box, Wm, chebyshev_node_all);

  std::size_t poly_ord_all = chebyshev_node_all.size();

  std::array<TYPE, DIM> chebyshev_real_pos;
  std::array<TYPE, DIM> relative_orig_pos;
  std::array<int, DIM> target_ind_of_box;
  std::array<int, DIM> ind_of_box;
  for(std::size_t t=0; t<target.size(); ++t){
    ans[t] = 0.0;

    for(std::size_t dim=0; dim<DIM; ++dim){
      target_ind_of_box[dim] = std::min((int)((target[t][dim]-min_pos[dim])/len), nrn_N-1);
    }

    for(std::size_t s=0; s<SIZE; ++s){
      std::size_t s_copy = s;
      int max_dist_from_t = 0;
      for(std::size_t dim=0; dim<DIM; ++dim){
        relative_orig_pos[DIM-1-dim] = len*(s_copy%nrn_N)+min_pos[DIM-1-dim];
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
