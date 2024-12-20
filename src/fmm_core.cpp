#include"../include/fmma/fmma.hpp"
#include"../include/fmma/math.hpp"
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

//{{{ multipole_calc_box_indices
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
template std::vector<std::size_t> FMMA<double, 3>::multipole_calc_box_indices(const std::array<int, 3>& box_ind, int N);
template std::vector<std::size_t> FMMA<double, 1>::multipole_calc_box_indices(const std::array<std::size_t, 1>& box_ind, int N);
template std::vector<std::size_t> FMMA<double, 2>::multipole_calc_box_indices(const std::array<std::size_t, 2>& box_ind, int N);
template std::vector<std::size_t> FMMA<double, 3>::multipole_calc_box_indices(const std::array<std::size_t, 3>& box_ind, int N);
template std::vector<std::size_t> FMMA<float, 1>::multipole_calc_box_indices(const std::array<int, 1>& box_ind, int N);
template std::vector<std::size_t> FMMA<float, 2>::multipole_calc_box_indices(const std::array<int, 2>& box_ind, int N);
template std::vector<std::size_t> FMMA<float, 3>::multipole_calc_box_indices(const std::array<int, 3>& box_ind, int N);
template std::vector<std::size_t> FMMA<float, 1>::multipole_calc_box_indices(const std::array<std::size_t, 1>& box_ind, int N);
template std::vector<std::size_t> FMMA<float, 2>::multipole_calc_box_indices(const std::array<std::size_t, 2>& box_ind, int N);
template std::vector<std::size_t> FMMA<float, 3>::multipole_calc_box_indices(const std::array<std::size_t, 3>& box_ind, int N);
//}}}

//{{{ exact_calc_box_indices
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
template std::vector<std::size_t> FMMA<double, 3>::exact_calc_box_indices(const std::array<int, 3>& box_ind, int N);
template std::vector<std::size_t> FMMA<float, 1>::exact_calc_box_indices(const std::array<int, 1>& box_ind, int N);
template std::vector<std::size_t> FMMA<float, 2>::exact_calc_box_indices(const std::array<int, 2>& box_ind, int N);
template std::vector<std::size_t> FMMA<float, 3>::exact_calc_box_indices(const std::array<int, 3>& box_ind, int N);
//}}}

//{{{ get_minmax
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
template void FMMA<double, 3>::get_minmax(const std::vector<std::array<double, 3>>& array1, const std::vector<std::array<double, 3>>& array2, std::array<double, 3>& min_array, std::array<double, 3>& max_array);
template void FMMA<float, 1>::get_minmax(const std::vector<std::array<float, 1>>& array1, const std::vector<std::array<float, 1>>& array2, std::array<float, 1>& min_array, std::array<float, 1>& max_array);
template void FMMA<float, 2>::get_minmax(const std::vector<std::array<float, 2>>& array1, const std::vector<std::array<float, 2>>& array2, std::array<float, 2>& min_array, std::array<float, 2>& max_array);
template void FMMA<float, 3>::get_minmax(const std::vector<std::array<float, 3>>& array1, const std::vector<std::array<float, 3>>& array2, std::array<float, 3>& min_array, std::array<float, 3>& max_array);
//}}}

//{{{ get_origin_and_length
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
template void FMMA<double, 3>::get_origin_and_length(const std::vector<std::array<double, 3>>& array1, const std::vector<std::array<double, 3>>& array2, std::array<double, 3>& origin, double& Len);
template void FMMA<float, 1>::get_origin_and_length(const std::vector<std::array<float, 1>>& array1, const std::vector<std::array<float, 1>>& array2, std::array<float, 1>& origin, float& Len);
template void FMMA<float, 2>::get_origin_and_length(const std::vector<std::array<float, 2>>& array1, const std::vector<std::array<float, 2>>& array2, std::array<float, 2>& origin, float& Len);
template void FMMA<float, 3>::get_origin_and_length(const std::vector<std::array<float, 3>>& array1, const std::vector<std::array<float, 3>>& array2, std::array<float, 3>& origin, float& Len);
//}}} 

//{{{ get_box_ind_of_ind
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
template std::array<std::size_t, 3> FMMA<double, 3>::get_box_ind_of_ind(std::size_t ind, int N);
template std::array<std::size_t, 1> FMMA<float, 1>::get_box_ind_of_ind(std::size_t ind, int N);
template std::array<std::size_t, 2> FMMA<float, 2>::get_box_ind_of_ind(std::size_t ind, int N);
template std::array<std::size_t, 3> FMMA<float, 3>::get_box_ind_of_ind(std::size_t ind, int N);
//}}}

//{{{ get_ind_of_box_ind
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
template std::size_t FMMA<double, 3>::get_ind_of_box_ind(const std::array<int, 3>& box_ind, int N);
template std::size_t FMMA<float, 1>::get_ind_of_box_ind(const std::array<int, 1>& box_ind, int N);
template std::size_t FMMA<float, 2>::get_ind_of_box_ind(const std::array<int, 2>& box_ind, int N);
template std::size_t FMMA<float, 3>::get_ind_of_box_ind(const std::array<int, 3>& box_ind, int N);
//}}}

//{{{ P2M
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
#pragma omp parallel for
    for(std::size_t k=0; k<poly_ord_all; ++k){
      std::size_t k_copy = k;
      for(std::size_t dim=0; dim<DIM; ++dim){
        std::size_t ord = k_copy%(poly_ord+1);
        chebyshev_node_all[k][DIM-1-dim] = chebyshev_node[ord];
        k_copy /= (poly_ord+1);
      }
    }

    std::vector<std::array<TYPE, DIM>> relative_pos(source.size());
    std::vector<std::size_t> pos_node(source.size(), 0);

    for(std::size_t s=0; s<source.size(); ++s){
      std::array<TYPE, DIM> source_pos = source[s];
      for(std::size_t dim=0; dim<DIM; ++dim){
        int tmp = std::min((int)((source_pos[dim]-min_pos[dim])/len), N-1);
        pos_node[s] *= N;
        pos_node[s] += tmp;
        relative_pos[s][dim] = std::max(std::min(2.0*((source_pos[dim]-min_pos[dim])/len-tmp)-1.0, 1.0), -1.0);
      }
      source_ind_in_box[pos_node[s]].push_back(s);
    }

#pragma omp parallel for
    for(std::size_t siz=0; siz<SIZE; ++siz){
      for(const auto s : source_ind_in_box[siz]){
        for(std::size_t k=0; k<poly_ord_all; ++k){
          TYPE val = 1.0;
          for(std::size_t dim=0; dim<DIM; ++dim){
            val *= SChebyshev(poly_ord+1, chebyshev_node_all[k][dim], relative_pos[s][dim]);
          }
          Wm[siz][k] += source_weight[s]*val;
        }
      }
    }

    return;
  }

template void FMMA<double, 1>::P2M(const std::vector<double>& source_weight, const std::vector<std::array<double, 1>>& source, int N, const std::array<double, 1>& min_pos, const double len, std::vector<std::vector<std::size_t>>& source_ind_in_box, std::vector<std::vector<double>>& Wm, std::vector<std::array<double, 1>>& chebyshev_node_all);
template void FMMA<double, 2>::P2M(const std::vector<double>& source_weight, const std::vector<std::array<double, 2>>& source, int N, const std::array<double, 2>& min_pos, const double len, std::vector<std::vector<std::size_t>>& source_ind_in_box, std::vector<std::vector<double>>& Wm, std::vector<std::array<double, 2>>& chebyshev_node_all);
template void FMMA<double, 3>::P2M(const std::vector<double>& source_weight, const std::vector<std::array<double, 3>>& source, int N, const std::array<double, 3>& min_pos, const double len, std::vector<std::vector<std::size_t>>& source_ind_in_box, std::vector<std::vector<double>>& Wm, std::vector<std::array<double, 3>>& chebyshev_node_all);
template void FMMA<float, 1>::P2M(const std::vector<float>& source_weight, const std::vector<std::array<float, 1>>& source, int N, const std::array<float, 1>& min_pos, const float len, std::vector<std::vector<std::size_t>>& source_ind_in_box, std::vector<std::vector<float>>& Wm, std::vector<std::array<float, 1>>& chebyshev_node_all);
template void FMMA<float, 2>::P2M(const std::vector<float>& source_weight, const std::vector<std::array<float, 2>>& source, int N, const std::array<float, 2>& min_pos, const float len, std::vector<std::vector<std::size_t>>& source_ind_in_box, std::vector<std::vector<float>>& Wm, std::vector<std::array<float, 2>>& chebyshev_node_all);
template void FMMA<float, 3>::P2M(const std::vector<float>& source_weight, const std::vector<std::array<float, 3>>& source, int N, const std::array<float, 3>& min_pos, const float len, std::vector<std::vector<std::size_t>>& source_ind_in_box, std::vector<std::vector<float>>& Wm, std::vector<std::array<float, 3>>& chebyshev_node_all);
//}}}

//{{{ M2M
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

  std::vector<std::vector<std::vector<TYPE>>> vals(poly_ord_all);
  for(std::size_t kin=0; kin<poly_ord_all; ++kin){
    vals[kin].resize(poly_ord_all);
    for(std::size_t kout=0; kout<poly_ord_all; ++kout){
      vals[kin][kout].resize(1<<DIM);
      for(std::size_t sind=0; sind<(1<<DIM); ++sind){
        vals[kin][kout][sind] = 1.0;
        std::size_t s_copy = sind;
        for(std::size_t dim=0; dim<DIM; ++dim){
          int shift = 2.0*(s_copy%2);
          shift -= 1.0;
          s_copy /= 2;

          vals[kin][kout][sind] *= SChebyshev(poly_ord+1, chebyshev_node_all[kout][DIM-1-dim], (chebyshev_node_all[kin][DIM-1-dim]+shift)/2);
        }
      }
    }
  }

  std::array<int, DIM> box_ind_out;
  std::array<TYPE, DIM> shift;
  for(std::size_t ind_in=0; ind_in<Wm_in.size(); ++ind_in){
    std::array<std::size_t, DIM> box_ind_in = get_box_ind_of_ind(ind_in, N);
    std::size_t sind = 0;
    for(std::size_t dim=0; dim<DIM; ++dim){
      box_ind_out[dim] = box_ind_in[dim]/2;
      shift[dim] = 2.0*(box_ind_in[dim]%2);
      shift[dim] -= 1.0;
      sind *= 2;
      sind += box_ind_in[dim]%2;
    }

    std::size_t ind_out = get_ind_of_box_ind(box_ind_out, N/2);

    for(std::size_t kout=0; kout<poly_ord_all; ++kout){
      for(std::size_t kin=0; kin<poly_ord_all; ++kin){
        TYPE tmp = Wm_in[ind_in][kin];
        tmp *= vals[kin][kout][sind];
        Wm_out[ind_out][kout] += tmp;
      }
    }
  }

  return;
}

template void FMMA<double, 1>::M2M(const std::size_t N, const std::vector<std::array<double, 1>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wm_in, std::vector<std::vector<double>>& Wm_out);
template void FMMA<double, 2>::M2M(const std::size_t N, const std::vector<std::array<double, 2>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wm_in, std::vector<std::vector<double>>& Wm_out);
template void FMMA<double, 3>::M2M(const std::size_t N, const std::vector<std::array<double, 3>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wm_in, std::vector<std::vector<double>>& Wm_out);
template void FMMA<float, 1>::M2M(const std::size_t N, const std::vector<std::array<float, 1>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wm_in, std::vector<std::vector<float>>& Wm_out);
template void FMMA<float, 2>::M2M(const std::size_t N, const std::vector<std::array<float, 2>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wm_in, std::vector<std::vector<float>>& Wm_out);
template void FMMA<float, 3>::M2M(const std::size_t N, const std::vector<std::array<float, 3>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wm_in, std::vector<std::vector<float>>& Wm_out);
//}}}

  //{{{ M2L
template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::M2L(const std::size_t N, const std::array<TYPE, DIM>& origin, const TYPE Len, const std::vector<std::array<TYPE, DIM>>& chebyshev_node_all, const std::vector<std::vector<TYPE>>& Wm, std::vector<std::vector<TYPE>>& Wl){
  TYPE len = Len/N;
  std::size_t poly_ord_all = chebyshev_node_all.size();
  Wl.resize(Wm.size());
  for(std::size_t ind=0; ind<Wl.size(); ++ind){
    Wl[ind].resize(poly_ord_all);
  }

  if(this->trans_sym_flag==1){
    std::size_t SIZE = 1;
    for(std::size_t dim=0; dim<DIM; ++dim){
      SIZE *= 7;
    }
    std::vector<std::vector<TYPE>> fns(SIZE, std::vector<TYPE>(poly_ord_all*poly_ord_all));
    //std::vector<TYPE> fns_matrix(poly_ord_all*SIZE*poly_ord_all); // kl x (in x km)
#pragma omp parallel
    for(std::size_t sind=0; sind<SIZE; ++sind){
      std::array<int, DIM> shift;
      std::size_t s_copy = sind;
      int max_dist_from_ind = 0;
      for(std::size_t dim=0; dim<DIM; ++dim){
        shift[DIM-1-dim] = s_copy%7-3;
        max_dist_from_ind = std::max(max_dist_from_ind, std::abs(shift[DIM-1-dim]));
        s_copy /= 7;
      }
      if(max_dist_from_ind <= 1){
        continue;
      }
      std::vector<std::array<TYPE, DIM>> relative_vector_l(poly_ord_all);
      std::vector<std::array<TYPE, DIM>> relative_vector_m(poly_ord_all);
      for(std::size_t k=0; k<poly_ord_all; ++k){
        for(std::size_t dim=0; dim<DIM; ++dim){
          relative_vector_l[k][dim] = len*(chebyshev_node_all[k][dim]/2);
          relative_vector_m[k][dim] = len*((TYPE)shift[dim]+chebyshev_node_all[k][dim]/2);
        }
      }
      for(std::size_t kl=0; kl<poly_ord_all; ++kl){
        for(std::size_t km=0; km<poly_ord_all; ++km){
          TYPE tmp = fn(relative_vector_l[kl], relative_vector_m[km]);
          fns[sind][kl*poly_ord_all+km] = tmp;
          //fns_matrix[kl*SIZE*poly_ord_all + sind*poly_ord_all + km] = tmp;
        }
      }
    }
#pragma omp parallel
    {
      std::size_t ind_begin = 0;
      std::size_t ind_end = Wl.size();
#if FMMA_USE_OPENMP
      ind_begin = Wl.size()*omp_get_thread_num()/omp_get_num_threads();
      ind_end = Wl.size()*(omp_get_thread_num()+1)/omp_get_num_threads();
#endif
      //std::vector<TYPE> Wm_vector(SIZE*poly_ord_all); // (ind x km)
      //std::vector<TYPE> Wm_zero(SIZE*poly_ord_all, 0.0);
      for(std::size_t ind=ind_begin; ind<ind_end; ++ind){
        std::array<std::size_t, DIM> l_box_ind = get_box_ind_of_ind(ind, N);
        std::vector<std::size_t> indices = multipole_calc_box_indices(l_box_ind, N);

        //for(std::size_t i=0; i<SIZE*poly_ord_all; ++i){
        //  Wm_vector[i] = 0.0;
        //}
        //copy(SIZE*poly_ord_all, Wm_zero.data(), Wm_vector.data());

        for(std::size_t i=0; i<indices.size(); ++i){
          std::array<std::size_t, DIM> m_box_ind = get_box_ind_of_ind(indices[i], N);
          int sind = 0;
          for(std::size_t dim=0; dim<DIM; ++dim){
            sind *= 7;
            sind += (int)m_box_ind[dim]+3-(int)l_box_ind[dim];
          }

          //for(std::size_t k=0; k<poly_ord_all; ++k){
          //  Wm_vector[sind*poly_ord_all + k] = Wm[indices[i]][k];
          //}
          //copy(poly_ord_all, Wm[indices[i]].data(), Wm_vector.data()+sind*poly_ord_all);

          matvec((TYPE)1.0, fns[sind], Wm[indices[i]], (TYPE)1.0, Wl[ind]);
        }
        //matvec(fns_matrix, Wm_vector, Wl[ind]);
      }
    }
  }else{
#pragma omp parallel
    {
      std::size_t ind_begin = 0;
      std::size_t ind_end = Wl.size();
#if FMMA_USE_OPENMP
      ind_begin = Wl.size()*omp_get_thread_num()/omp_get_num_threads();
      ind_end = Wl.size()*(omp_get_thread_num()+1)/omp_get_num_threads();
#endif
      for(std::size_t ind=ind_begin; ind<ind_end; ++ind){
        std::array<std::size_t, DIM> l_box_ind = get_box_ind_of_ind(ind, N);
        std::vector<std::size_t> indices = multipole_calc_box_indices(l_box_ind, N);
        std::vector<std::array<TYPE, DIM>> relative_vector_l(poly_ord_all);
        for(std::size_t kl=0; kl<poly_ord_all; ++kl){
          for(std::size_t dim=0; dim<DIM; ++dim){
            relative_vector_l[kl][dim] = len*((TYPE)l_box_ind[dim]+(1.0+chebyshev_node_all[kl][dim])/2)+origin[dim];
          }
        }
        std::vector<std::array<TYPE, DIM>> relative_vector_m(poly_ord_all);
        std::vector<TYPE> fns(poly_ord_all*poly_ord_all);
        for(std::size_t i=0; i<indices.size(); ++i){
          std::array<std::size_t, DIM> m_box_ind = get_box_ind_of_ind(indices[i], N);
          for(std::size_t km=0; km<poly_ord_all; ++km){
            for(std::size_t dim=0; dim<DIM; ++dim){
              relative_vector_m[km][dim] = len*((TYPE)m_box_ind[dim]+(1.0+chebyshev_node_all[km][dim])/2)+origin[dim];
            }
          }
          for(std::size_t kl=0; kl<poly_ord_all; ++kl){
            for(std::size_t km=0; km<poly_ord_all; ++km){
              fns[kl*poly_ord_all+km] = fn(relative_vector_l[kl], relative_vector_m[km]);
            }
          }
          matvec((TYPE)1.0, fns, Wm[indices[i]], (TYPE)1.0, Wl[ind]);
        }
      }
    }
  }

  return;
}

template void FMMA<double, 1>::M2L(const std::size_t N, const std::array<double, 1>& origin, const double Len, const std::vector<std::array<double, 1>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wm, std::vector<std::vector<double>>& Wl);
template void FMMA<double, 2>::M2L(const std::size_t N, const std::array<double, 2>& origin, const double Len, const std::vector<std::array<double, 2>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wm, std::vector<std::vector<double>>& Wl);
template void FMMA<double, 3>::M2L(const std::size_t N, const std::array<double, 3>& origin, const double Len, const std::vector<std::array<double, 3>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wm, std::vector<std::vector<double>>& Wl);
template void FMMA<float, 1>::M2L(const std::size_t N, const std::array<float, 1>& origin, const float Len, const std::vector<std::array<float, 1>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wm, std::vector<std::vector<float>>& Wl);
template void FMMA<float, 2>::M2L(const std::size_t N, const std::array<float, 2>& origin, const float Len, const std::vector<std::array<float, 2>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wm, std::vector<std::vector<float>>& Wl);
template void FMMA<float, 3>::M2L(const std::size_t N, const std::array<float, 3>& origin, const float Len, const std::vector<std::array<float, 3>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wm, std::vector<std::vector<float>>& Wl);
//}}}

//{{{ L2L
template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::L2L(const std::size_t N, const std::vector<std::array<TYPE, DIM>>& chebyshev_node_all, const std::vector<std::vector<TYPE>>& Wl_in, std::vector<std::vector<TYPE>>& Wl_out){
  {
    std::size_t Nd = 1;
    for(std::size_t dim=0; dim<DIM; ++dim){
      Nd *= (N/2);
    }
    if(Wl_in.size() != Nd){
      fprintf(stderr, "%s:%d ERROR : size error %zu != %zu\n", __FILE__, __LINE__, Wl_in.size(), Nd);
      exit(EXIT_FAILURE);
    }
    Wl_out.resize(Wl_in.size()<<DIM);
  }
  std::size_t poly_ord_all = chebyshev_node_all.size();
  for(std::size_t i=0; i<Wl_out.size(); ++i){
    Wl_out[i].resize(poly_ord_all);
  }

  std::vector<std::vector<std::vector<TYPE>>> vals(poly_ord_all);
  for(std::size_t kin=0; kin<poly_ord_all; ++kin){
    vals[kin].resize(poly_ord_all);
    for(std::size_t kout=0; kout<poly_ord_all; ++kout){
      vals[kin][kout].resize(1<<DIM);
      for(std::size_t sind=0; sind<(1<<DIM); ++sind){
        vals[kin][kout][sind] = 1.0;
        std::size_t s_copy = sind;
        for(std::size_t dim=0; dim<DIM; ++dim){
          int shift = 2.0*(s_copy%2);
          shift -= 1.0;
          s_copy /= 2;

          vals[kin][kout][sind] *= SChebyshev(poly_ord+1, (chebyshev_node_all[kout][DIM-1-dim]+shift)/2, chebyshev_node_all[kin][DIM-1-dim]);
        }
      }
    }
  }

  std::array<int, DIM> box_ind_in;
  std::array<TYPE, DIM> shift;
  for(std::size_t ind_out=0; ind_out<Wl_out.size(); ++ind_out){
    std::array<std::size_t, DIM> box_ind_out = get_box_ind_of_ind(ind_out, N);
    std::size_t sind = 0;
    for(std::size_t dim=0; dim<DIM; ++dim){
      box_ind_in[dim] = box_ind_out[dim]/2;
      shift[dim] = 2.0*(box_ind_out[dim]%2);
      shift[dim] -= 1.0;
      sind *= 2;
      sind += box_ind_out[dim]%2;
    }

    std::size_t ind_in = get_ind_of_box_ind(box_ind_in, N/2);

    for(std::size_t kout=0; kout<poly_ord_all; ++kout){
      for(std::size_t kin=0; kin<poly_ord_all; ++kin){
        TYPE tmp = Wl_in[ind_in][kin];
        tmp *= vals[kin][kout][sind];
        Wl_out[ind_out][kout] += tmp;
      }
    }
  }

  return;
}

template void FMMA<double, 1>::L2L(const std::size_t N, const std::vector<std::array<double, 1>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wl_in, std::vector<std::vector<double>>& Wl_out);
template void FMMA<double, 2>::L2L(const std::size_t N, const std::vector<std::array<double, 2>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wl_in, std::vector<std::vector<double>>& Wl_out);
template void FMMA<double, 3>::L2L(const std::size_t N, const std::vector<std::array<double, 3>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wl_in, std::vector<std::vector<double>>& Wl_out);
template void FMMA<float, 1>::L2L(const std::size_t N, const std::vector<std::array<float, 1>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wl_in, std::vector<std::vector<float>>& Wl_out);
template void FMMA<float, 2>::L2L(const std::size_t N, const std::vector<std::array<float, 2>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wl_in, std::vector<std::vector<float>>& Wl_out);
template void FMMA<float, 3>::L2L(const std::size_t N, const std::vector<std::array<float, 3>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wl_in, std::vector<std::vector<float>>& Wl_out);
//}}}

//{{{ L2P
template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::L2P(const std::vector<std::array<TYPE, DIM>>& target, const std::array<TYPE, DIM>& origin, const TYPE Len, const std::vector<std::array<TYPE, DIM>>& chebyshev_node_all, const std::vector<std::vector<TYPE>>& Wl, std::vector<TYPE>& ans){

  std::vector<std::array<TYPE, DIM>> relative_pos(target.size());
  std::vector<std::size_t> ind(target.size());
  std::size_t poly_ord_all = chebyshev_node_all.size();

#pragma omp parallel for
  for(std::size_t t=0; t<target.size(); ++t){
    ans[t] = 0.0;

    std::size_t N = 1<<(Depth-1);
    std::array<int, DIM> target_ind_of_box;

    TYPE len = Len/N;
    for(std::size_t dim=0; dim<DIM; ++dim){
      target_ind_of_box[dim] = std::min((int)((target[t][dim]-origin[dim])/len), (int)N-1);
      relative_pos[t][dim] = std::max(std::min(2.0*((target[t][dim]-origin[dim])/len-target_ind_of_box[dim])-1.0, 1.0), -1.0);
    }

    ind[t] = get_ind_of_box_ind(target_ind_of_box, N);
  }

#pragma omp parallel for
  for(std::size_t t=0; t<target.size(); ++t){
    std::vector<TYPE> vals(poly_ord_all);
    for(std::size_t k=0; k<poly_ord_all; ++k){
      vals[k] = 1.0;
      for(std::size_t dim=0; dim<DIM; ++dim){
        vals[k] *= SChebyshev(poly_ord+1, relative_pos[t][dim], chebyshev_node_all[k][dim]);
      }
    }
    ans[t] += dot(poly_ord_all, Wl[ind[t]].data(), vals.data());
  }

  return;
}

template void FMMA<double, 1>::L2P(const std::vector<std::array<double, 1>>& target, const std::array<double, 1>& origin, const double Len, const std::vector<std::array<double, 1>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wl, std::vector<double>& ans);
template void FMMA<double, 2>::L2P(const std::vector<std::array<double, 2>>& target, const std::array<double, 2>& origin, const double Len, const std::vector<std::array<double, 2>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wl, std::vector<double>& ans);
template void FMMA<double, 3>::L2P(const std::vector<std::array<double, 3>>& target, const std::array<double, 3>& origin, const double Len, const std::vector<std::array<double, 3>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wl, std::vector<double>& ans);
template void FMMA<float, 1>::L2P(const std::vector<std::array<float, 1>>& target, const std::array<float, 1>& origin, const float Len, const std::vector<std::array<float, 1>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wl, std::vector<float>& ans);
template void FMMA<float, 2>::L2P(const std::vector<std::array<float, 2>>& target, const std::array<float, 2>& origin, const float Len, const std::vector<std::array<float, 2>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wl, std::vector<float>& ans);
template void FMMA<float, 3>::L2P(const std::vector<std::array<float, 3>>& target, const std::array<float, 3>& origin, const float Len, const std::vector<std::array<float, 3>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wl, std::vector<float>& ans);
//}}}

//{{{ M2P
template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::M2P(const std::vector<std::array<TYPE, DIM>>& target, const std::size_t N, const std::array<TYPE, DIM>& origin, const TYPE Len, const std::vector<std::array<TYPE, DIM>>& chebyshev_node_all, const std::vector<std::vector<TYPE>>& Wm, std::vector<TYPE>& ans){
  std::size_t poly_ord_all = chebyshev_node_all.size();
  TYPE len = Len/N;

  auto ansp = ans.data();

#pragma omp parallel for
  for(std::size_t t=0; t<target.size(); ++t){
    std::array<int, DIM> target_ind_of_box;
    for(std::size_t dim=0; dim<DIM; ++dim){
      target_ind_of_box[dim] = std::min((int)((target[t][dim]-origin[dim])/len), (int)N-1);
    }

    std::vector<std::size_t> indices = multipole_calc_box_indices(target_ind_of_box, N);
    std::array<TYPE, DIM> relative_orig_pos;
    std::array<TYPE, DIM> chebyshev_real_pos;

    for(std::size_t i=0; i<indices.size(); ++i){
      std::size_t s = indices[i];
      std::size_t s_copy = s;
      for(std::size_t dim=0; dim<DIM; ++dim){
        relative_orig_pos[DIM-1-dim] = len*(s_copy%N)+origin[DIM-1-dim];
        s_copy /= N;
      }

      for(std::size_t k=0; k<poly_ord_all; ++k){
        for(std::size_t dim=0; dim<DIM; ++dim){
          chebyshev_real_pos[dim] = (chebyshev_node_all[k][dim]+1.0)/2.0*len+relative_orig_pos[dim];
        }
        ansp[t] += fn(target[t], chebyshev_real_pos)*Wm[s][k];
      }
    }
  }
  return;
};

template void FMMA<double, 1>::M2P(const std::vector<std::array<double, 1>>& target, const std::size_t N, const std::array<double, 1>& origin, const double Len, const std::vector<std::array<double, 1>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wm, std::vector<double>& ans);
template void FMMA<double, 2>::M2P(const std::vector<std::array<double, 2>>& target, const std::size_t N, const std::array<double, 2>& origin, const double Len, const std::vector<std::array<double, 2>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wm, std::vector<double>& ans);
template void FMMA<double, 3>::M2P(const std::vector<std::array<double, 3>>& target, const std::size_t N, const std::array<double, 3>& origin, const double Len, const std::vector<std::array<double, 3>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wm, std::vector<double>& ans);
template void FMMA<float, 1>::M2P(const std::vector<std::array<float, 1>>& target, const std::size_t N, const std::array<float, 1>& origin, const float Len, const std::vector<std::array<float, 1>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wm, std::vector<float>& ans);
template void FMMA<float, 2>::M2P(const std::vector<std::array<float, 2>>& target, const std::size_t N, const std::array<float, 2>& origin, const float Len, const std::vector<std::array<float, 2>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wm, std::vector<float>& ans);
template void FMMA<float, 3>::M2P(const std::vector<std::array<float, 3>>& target, const std::size_t N, const std::array<float, 3>& origin, const float Len, const std::vector<std::array<float, 3>>& chebyshev_node_all, const std::vector<std::vector<float>>& Wm, std::vector<float>& ans);
//}}}

}; // namespace fmma
