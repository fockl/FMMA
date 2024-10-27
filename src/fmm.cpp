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
void FMMA<TYPE, DIM>::M2L(const std::size_t N, const TYPE Len, const std::vector<std::array<TYPE, DIM>>& chebyshev_node_all, const std::vector<std::vector<TYPE>>& Wm, std::vector<std::vector<TYPE>>& Wl){
  TYPE len = Len/N;
  std::array<TYPE, DIM> relative_vector;
  std::size_t poly_ord_all = chebyshev_node_all.size();
  Wl.resize(Wm.size());
  for(std::size_t ind=0; ind<Wl.size(); ++ind){
    Wl[ind].resize(poly_ord_all);
    for(std::size_t kl=0; kl<poly_ord_all; ++kl){
      Wl[ind][kl] = 0.0;
    }
    std::array<std::size_t, DIM> l_box_ind = get_box_ind_of_ind(ind, N);
    std::vector<std::size_t> indices = multipole_calc_box_indices(l_box_ind, N);
    for(std::size_t i=0; i<indices.size(); ++i){
      std::array<std::size_t, DIM> m_box_ind = get_box_ind_of_ind(indices[i], N);
      for(std::size_t kl=0; kl<poly_ord_all; ++kl){
        for(std::size_t km=0; km<poly_ord_all; ++km){
          for(std::size_t dim=0; dim<DIM; ++dim){
            relative_vector[dim] = len*(((TYPE)l_box_ind[dim]+chebyshev_node_all[kl][dim]/2)-((TYPE)m_box_ind[dim]+chebyshev_node_all[km][dim]/2));
          }
          Wl[ind][kl] += fn(relative_vector)*Wm[indices[i]][km];
        }
      }
    }
  }

  return;
}

template void FMMA<double, 1>::M2L(const std::size_t N, const double Len, const std::vector<std::array<double, 1>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wm, std::vector<std::vector<double>>& Wl);
template void FMMA<double, 2>::M2L(const std::size_t N, const double Len, const std::vector<std::array<double, 2>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wm, std::vector<std::vector<double>>& Wl);

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
    for(std::size_t k=0; k<poly_ord_all; ++k){
      //Wl_out[i][k] = 0.0;
    }
  }

  std::array<int, DIM> box_ind_in;
  std::array<TYPE, DIM> shift;
  for(std::size_t ind_out=0; ind_out<Wl_out.size(); ++ind_out){
    std::array<std::size_t, DIM> box_ind_out = get_box_ind_of_ind(ind_out, N);
    for(std::size_t dim=0; dim<DIM; ++dim){
      box_ind_in[dim] = box_ind_out[dim]/2;
      shift[dim] = 2.0*(box_ind_out[dim]%2);
      shift[dim] -= 1.0;
    }

    std::size_t ind_in = get_ind_of_box_ind(box_ind_in, N/2);

    for(std::size_t kout=0; kout<poly_ord_all; ++kout){
      for(std::size_t kin=0; kin<poly_ord_all; ++kin){
        TYPE tmp = Wl_in[ind_in][kin];
        for(std::size_t dim=0; dim<DIM; ++dim){
          tmp *= SChebyshev(poly_ord+1, (chebyshev_node_all[kout][dim]+shift[dim])/2, chebyshev_node_all[kin][dim]);
        }
        Wl_out[ind_out][kout] += tmp;
      }
    }
  }

  return;
}

template void FMMA<double, 1>::L2L(const std::size_t N, const std::vector<std::array<double, 1>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wl_in, std::vector<std::vector<double>>& Wl_out);
template void FMMA<double, 2>::L2L(const std::size_t N, const std::vector<std::array<double, 2>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wl_in, std::vector<std::vector<double>>& Wl_out);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::L2P(const std::array<TYPE, DIM>& target, const std::size_t N, const std::array<TYPE, DIM>& origin, const TYPE Len, const std::vector<std::array<TYPE, DIM>>& chebyshev_node_all, const std::vector<std::vector<TYPE>>& Wl, std::array<int, DIM>& target_ind_of_box, TYPE& ans){
  TYPE len = Len/N;
  std::size_t poly_ord_all = chebyshev_node_all.size();
  std::array<TYPE, DIM> relative_pos;
  for(std::size_t dim=0; dim<DIM; ++dim){
    target_ind_of_box[dim] = std::min((int)((target[dim]-origin[dim])/len), (int)N-1);
    relative_pos[dim] = std::max(std::min(2.0*((target[dim]-origin[dim])/len-target_ind_of_box[dim])-1.0, 1.0), -1.0);
  }

  std::size_t ind = get_ind_of_box_ind(target_ind_of_box, N);

  for(std::size_t k=0; k<poly_ord_all; ++k){
    TYPE val = 1.0;
    for(std::size_t dim=0; dim<DIM; ++dim){
      val *= SChebyshev(poly_ord+1, relative_pos[dim], chebyshev_node_all[k][dim]);
    }
    ans += Wl[ind][k]*val;
  }

  return;
}

template void FMMA<double, 1>::L2P(const std::array<double, 1>& target, const std::size_t N, const std::array<double, 1>& origin, const double Len, const std::vector<std::array<double, 1>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wl, std::array<int, 1>& target_ind_of_box, double& ans);
template void FMMA<double, 2>::L2P(const std::array<double, 2>& target, const std::size_t N, const std::array<double, 2>& origin, const double Len, const std::vector<std::array<double, 2>>& chebyshev_node_all, const std::vector<std::vector<double>>& Wl, std::array<int, 2>& target_ind_of_box, double& ans);

template<typename TYPE, std::size_t DIM>
void FMMA<TYPE, DIM>::fmm(const std::vector<std::array<TYPE, DIM>>& target, const std::vector<TYPE>& source_weight, const std::vector<std::array<TYPE, DIM>>& source, std::vector<TYPE>& ans){

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

  std::vector<std::vector<std::vector<TYPE>>> Wl(Depth);
  {
    // M2L
    std::size_t tmp_N = 1;
    for(int depth=0; depth<Depth; ++depth){
      M2L(tmp_N, Len, chebyshev_node_all, Wm[depth], Wl[depth]);
      tmp_N *= 2;
    }

    // L2L
    tmp_N = 1;
    for(int depth=0; depth+1<Depth; ++depth){
      L2L(tmp_N*2, chebyshev_node_all, Wl[depth], Wl[depth+1]);
      tmp_N *= 2;
    }
  }


  for(std::size_t t=0; t<target.size(); ++t){
    ans[t] = 0.0;

    std::size_t tmp_N = 1<<(Depth-1);
    std::array<int, DIM> target_ind_of_box;

    L2P(target[t], tmp_N, origin, Len, chebyshev_node_all, Wl[Depth-1], target_ind_of_box, ans[t]);

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

template void FMMA<double, 1>::fmm(const std::vector<std::array<double, 1>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 1>>& source, std::vector<double>& ans);
template void FMMA<double, 2>::fmm(const std::vector<std::array<double, 2>>& target, const std::vector<double>& source_weight, const std::vector<std::array<double, 2>>& source, std::vector<double>& ans);

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

}; // namespace fmma
