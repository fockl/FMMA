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
std::array<std::size_t, DIM> FMMA<TYPE, DIM>::get_box_ind_of_ind(const std::size_t ind, int N){
  std::array<std::size_t, DIM> ans;
  std::size_t ind_copy = ind;
  for(std::size_t dim=0; dim<DIM; ++dim){
    ans[DIM-1-dim] = ind_copy%N;
    ind_copy /= N;
  }
  return ans;
}

template std::array<std::size_t, 1> FMMA<double, 1>::get_box_ind_of_ind(std::size_t ind, int N);
template std::array<std::size_t, 2> FMMA<double, 2>::get_box_ind_of_ind(std::size_t ind, int N);

template<typename TYPE, std::size_t DIM>
std::size_t FMMA<TYPE, DIM>::get_ind_of_box_ind(const std::array<int, DIM>& box_ind, int N){
  std::size_t ans = 0;
  for(std::size_t dim=0; dim<DIM; ++dim){
    ans *= N;
    ans += box_ind[dim];
  }
  return ans;
}

template std::size_t FMMA<double, 1>::get_ind_of_box_ind(const std::array<int, 1>& box_ind, int N);
template std::size_t FMMA<double, 2>::get_ind_of_box_ind(const std::array<int, 2>& box_ind, int N);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::M2M(const std::size_t N, const std::vector<std::array<TYPE, DIM>>& chebyshev_node_all, const std::vector<std::vector<TYPE>>& Wm_in, std::vector<std::vector<TYPE>>& Wm_out){
  {
    std::size_t Nd = 1;
    for(std::size_t dim=0; dim<DIM; ++dim){
      Nd *= N;
    }
    if(Wm_in.size() != Nd){
      fprintf(stderr, "%s:%d ERROR : size error %zu != %zu\n", __FILE__, __LINE__, Wm_in.size(), Nd);
      exit(EXIT_FAILURE);
    }
    Wm_out.resize(Wm_in.size()>>DIM);
    if(Wm_in.size() != Wm_out.size()<<DIM){
      fprintf(stderr, "%s:%d ERROR : size error %zu != %zu\n", __FILE__, __LINE__, Wm_in.size(), Wm_out.size()<<DIM);
      exit(EXIT_FAILURE);
    }
  }
  std::size_t poly_ord_all = chebyshev_node_all.size();
  for(std::size_t i=0; i<Wm_out.size(); ++i){
    Wm_out[i].resize(poly_ord_all);
    for(std::size_t k=0; k<poly_ord_all; ++k){
      Wm_out[i][k] = 0.0;
    }
  }

  std::array<int, DIM> box_ind_out;
  std::array<TYPE, DIM> shift;
  for(std::size_t ind_in=0; ind_in<Wm_in.size(); ++ind_in){
    std::array<std::size_t, DIM> box_ind_in = get_box_ind_of_ind(ind_in, N);
    for(std::size_t dim=0; dim<DIM; ++dim){
      box_ind_out[dim] = box_ind_in[dim]/2;
      shift[dim] = 2.0*(box_ind_in[dim]%2);
      shift[dim] -= 1.0;
    }

    std::size_t ind_out = get_ind_of_box_ind(box_ind_out, N/2);

    for(std::size_t kout=0; kout<poly_ord_all; ++kout){
      for(std::size_t kin=0; kin<poly_ord_all; ++kin){
        TYPE tmp = Wm_in[ind_in][kin];
        for(std::size_t dim=0; dim<DIM; ++dim){
          tmp *= SChebyshev(poly_ord+1, chebyshev_node_all[kout][dim], (chebyshev_node_all[kin][dim]+shift[dim])/2);
        }
        Wm_out[ind_out][kout] += tmp;
      }
    }
  }

  return;
}

template void FMMA<double, 1>::M2M(const std::size_t N, const std::vector<std::array<double, 1>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wm_in, std::vector<std::vector<double>>& Wm_out);
template void FMMA<double, 2>::M2M(const std::size_t N, const std::vector<std::array<double, 2>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wm_in, std::vector<std::vector<double>>& Wm_out);


template<typename TYPE, std::size_t DIM>
template<typename UINT>
std::vector<std::size_t> FMMA<TYPE, DIM>::multipole_calc_box_indices(const std::array<UINT, DIM>& box_ind, int N){
  std::vector<std::size_t> ans;
  std::array<int, DIM> lower, upper;
  for(std::size_t dim=0; dim<DIM; ++dim){
    if(box_ind[dim]/2>0){
      lower[dim] = box_ind[dim]/2;
      lower[dim] -= 1;
    }else{
      lower[dim] = box_ind[dim]/2;
    }
    if(box_ind[dim]/2+1<(UINT)N/2){
      upper[dim] = box_ind[dim]/2;
      upper[dim] += 1;
    }else{
      upper[dim] = box_ind[dim]/2;
    }
    lower[dim] *= 2;
    upper[dim] *= 2;
    upper[dim] += 1;
  }
  std::array<std::size_t, DIM> shape;
  std::size_t SIZE = 1;
  for(std::size_t dim=0; dim<DIM; ++dim){
    shape[dim] = upper[dim]-lower[dim]+1;
    SIZE *= shape[dim];
  }

  for(std::size_t s=0; s<SIZE; ++s){
    std::size_t s_copy = s;
    int max_dist_from_ind = 0;
    std::array<int, DIM> pos = lower;
    for(std::size_t dim=0; dim<DIM; ++dim){
      pos[DIM-1-dim] += s_copy%shape[DIM-1-dim];
      max_dist_from_ind = std::max(max_dist_from_ind, std::abs(pos[DIM-1-dim]-(int)box_ind[DIM-1-dim]));
      s_copy /= shape[DIM-1-dim];
    }

    if(max_dist_from_ind <= 1){
      continue;
    }

    ans.push_back(get_ind_of_box_ind(pos, N));
  }

  return ans;
}

template std::vector<std::size_t> FMMA<double, 1>::multipole_calc_box_indices(const std::array<int, 1>& box_ind, int N);
template std::vector<std::size_t> FMMA<double, 2>::multipole_calc_box_indices(const std::array<int, 2>& box_ind, int N);
template std::vector<std::size_t> FMMA<double, 1>::multipole_calc_box_indices(const std::array<std::size_t, 1>& box_ind, int N);
template std::vector<std::size_t> FMMA<double, 2>::multipole_calc_box_indices(const std::array<std::size_t, 2>& box_ind, int N);


template<typename TYPE, std::size_t DIM>
std::vector<std::size_t> FMMA<TYPE, DIM>::exact_calc_box_indices(const std::array<int, DIM>& box_ind, int N){
  std::vector<std::size_t> ans;
  std::array<int, DIM> lower, upper;
  for(std::size_t dim=0; dim<DIM; ++dim){
    if(box_ind[dim]>0){
      lower[dim] = box_ind[dim]-1;
    }else{
      lower[dim] = box_ind[dim];
    }
    if(box_ind[dim]+1<N){
      upper[dim] = box_ind[dim]+1;
    }else{
      upper[dim] = box_ind[dim];
    }
  }
  std::array<std::size_t, DIM> shape;
  std::size_t SIZE = 1;
  for(std::size_t dim=0; dim<DIM; ++dim){
    shape[dim] = upper[dim]-lower[dim]+1;
    SIZE *= shape[dim];
  }

  for(std::size_t s=0; s<SIZE; ++s){
    std::size_t s_copy = s;
    std::array<int, DIM> pos = lower;
    for(std::size_t dim=0; dim<DIM; ++dim){
      pos[DIM-1-dim] += s_copy%shape[DIM-1-dim];
      s_copy /= shape[DIM-1-dim];
    }

    ans.push_back(get_ind_of_box_ind(pos, N));
  }

  return ans;
}

template std::vector<std::size_t> FMMA<double, 1>::exact_calc_box_indices(const std::array<int, 1>& box_ind, int N);
template std::vector<std::size_t> FMMA<double, 2>::exact_calc_box_indices(const std::array<int, 2>& box_ind, int N);


template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::tree(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){

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
    P2M(source_weight, source, tmp_N, origin, Len/tmp_N, source_ind_in_box, Wm[Depth-1], chebyshev_node_all);
    // M2M
    for(int depth=0; depth+1<Depth; ++depth){
      M2M(tmp_N, chebyshev_node_all, Wm[Depth-depth-1], Wm[Depth-depth-2]);
      tmp_N /= 2;
    }
  }


  std::array<TYPE, DIM> chebyshev_real_pos;
  std::array<TYPE, DIM> relative_orig_pos;
  std::array<int, DIM> target_ind_of_box;
  for(std::size_t t=0; t<target.size(); ++t){
    ans[t] = 0.0;

    // M2P
    int tmp_N = 1;
    for(int depth=0; depth<Depth; ++depth){
      std::size_t poly_ord_all = chebyshev_node_all.size();
      TYPE len = Len/tmp_N;

      for(std::size_t dim=0; dim<DIM; ++dim){
        target_ind_of_box[dim] = std::min((int)((target[t][dim]-origin[dim])/len), tmp_N-1);
      }

      std::vector<std::size_t> indices = multipole_calc_box_indices(target_ind_of_box, tmp_N);

      for(std::size_t i=0; i<indices.size(); ++i){
        std::size_t s = indices[i];
        std::size_t s_copy = s;
        for(std::size_t dim=0; dim<DIM; ++dim){
          relative_orig_pos[DIM-1-dim] = len*(s_copy%tmp_N)+origin[DIM-1-dim];
          s_copy /= tmp_N;
        }

        for(std::size_t k=0; k<poly_ord_all; ++k){
          for(std::size_t dim=0; dim<DIM; ++dim){
            chebyshev_real_pos[dim] = (chebyshev_node_all[k][dim]+1.0)/2.0*len+relative_orig_pos[dim];
          }
          ans[t] += fn(target[t]-chebyshev_real_pos)*Wm[depth][s][k];
        }
      }
      tmp_N *= 2;
    } // depth
    tmp_N /= 2;

    std::vector<std::size_t> indices = exact_calc_box_indices(target_ind_of_box, tmp_N);
    // P2P
    for(std::size_t i=0; i<indices.size(); ++i){
      for(std::size_t j=0; j<source_ind_in_box[indices[i]].size(); ++j){
        ans[t] += source_weight[source_ind_in_box[indices[i]][j]]*fn(target[t]-source[source_ind_in_box[indices[i]][j]]);
      }
    }
  }

  return;
};

template void FMMA<double, 1>::tree(const std::vector<std::array<double, 1>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::tree(const std::vector<std::array<double, 2>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);

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

}; // namespace fmma
