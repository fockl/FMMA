#include"../../../include/fmma/math.hpp"
#include"../../../include/fmma/fmma.hpp"
#include"../test_common.hpp"
#include<cstdlib>
#include<cstdio>
#include<vector>
#include<array>
#include<functional>
#include<cmath>

namespace fmma {

  template<typename TYPE, std::size_t DIM>
class FMMA_TEST : public FMMA<TYPE, DIM>::FMMA{

  public:
bool test_m2m(std::size_t M, std::size_t N, TYPE tol){

  int Depth = 4;
  std::vector<std::array<TYPE, DIM>> target(M);
  std::vector<std::array<TYPE, DIM>> source(N);
  std::vector<TYPE> source_weight(N);

  for(std::size_t m=0; m<M; ++m){
    for(std::size_t dim=0; dim<DIM; ++dim){
      target[m][dim] = (TYPE)rand()/RAND_MAX-0.5;
    }
  }
  for(std::size_t n=0; n<N; ++n){
    for(std::size_t dim=0; dim<DIM; ++dim){
      source[n][dim] = (TYPE)rand()/RAND_MAX-0.5;
    }
    source_weight[n] = (TYPE)rand()/RAND_MAX-0.5;
  }

  std::array<TYPE, DIM> origin;
  TYPE Len;
  this->get_origin_and_length(target, source, origin, Len);

  std::vector<std::vector<std::vector<TYPE>>> Wm_P2M(Depth);
  std::vector<std::vector<std::vector<TYPE>>> Wm_M2M(Depth);
  std::vector<std::vector<std::vector<std::size_t>>> source_ind_in_box(Depth);
  std::vector<std::vector<std::array<TYPE, DIM>>> chebyshev_node_all(Depth);

  {
    std::size_t tmp_N = 1;
    for(int depth=0; depth<Depth; ++depth){
      this->P2M(source_weight, source, tmp_N, origin, Len/tmp_N, source_ind_in_box[depth], Wm_P2M[depth], chebyshev_node_all[depth]);
      tmp_N *= 2;
    }
  }

  {
    Wm_M2M[Depth-1] = Wm_P2M[Depth-1];
    std::size_t tmp_N = 1<<Depth;
    for(int depth=0; depth+1<Depth; ++depth){
      tmp_N /= 2;
      this->M2M(tmp_N, chebyshev_node_all[Depth-depth-1], Wm_M2M[Depth-depth-1], Wm_M2M[Depth-depth-2]);
    }
  }

  std::string comment = "<" + get_type<TYPE>() + "," + std::to_string(DIM) + ">(" + std::to_string(M) + ", " + std::to_string(N) + ", " + std::to_string(tol) + ")";

  for(int depth=0; depth<Depth; ++depth){
    TYPE diff = 0.0;
    for(std::size_t i=0; i<Wm_P2M[depth].size(); ++i){
      for(std::size_t j=0; j<Wm_P2M[depth][i].size(); ++j){
        TYPE tmp = Wm_P2M[depth][i][j] - Wm_M2M[depth][i][j];
        diff += tmp*tmp;
      }
    }

    if(diff >= tol){
      fprintf(stderr, "diff[%d] = %e > tol = %e\n", depth, diff, tol);
      failed(__FILE__, __func__, comment);
      return false;
    }
  }

  pass(__FILE__, __func__, comment);
  return true;
}

};

}; // namespace fmma

int main(void){
  srand(0);

  {
    fmma::FMMA_TEST<double, 1> fmma_test;
    if(!fmma_test.test_m2m(10, 10, 1.0e-6)){
      exit(EXIT_FAILURE);
    }
    if(!fmma_test.test_m2m(10, 20, 1.0e-6)){
      exit(EXIT_FAILURE);
    }
    if(!fmma_test.test_m2m(20, 10, 1.0e-6)){
      exit(EXIT_FAILURE);
    }
  }

  {
    fmma::FMMA_TEST<double, 2> fmma_test;
    if(!fmma_test.test_m2m(10, 10, 1.0e-6)){
      exit(EXIT_FAILURE);
    }
    if(!fmma_test.test_m2m(10, 20, 1.0e-6)){
      exit(EXIT_FAILURE);
    }
    if(!fmma_test.test_m2m(20, 10, 1.0e-6)){
      exit(EXIT_FAILURE);
    }
  }

  return 0;
}
