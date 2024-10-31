#include"../../../include/fmma/fmma.hpp"
#include"../test_common.hpp"
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<vector>
#include<array>

//{{{ test_tree
template<typename TYPE, std::size_t DIM>
bool test_tree(int ssize, int tsize, TYPE tol){
  std::vector<std::array<TYPE, DIM>> source(ssize), target(tsize);
  for(int i=0; i<ssize; ++i){
    for(std::size_t dim=0; dim<DIM; ++dim){
      source[i][dim] = (TYPE)rand()/RAND_MAX-0.5;
    }
  }
  for(int i=0; i<tsize; ++i){
    for(std::size_t dim=0; dim<DIM; ++dim){
      target[i][dim] = (TYPE)rand()/RAND_MAX-0.5;
    }
  }

  fmma::FMMA<TYPE, DIM> fmma;
  std::vector<TYPE> ans(tsize);
  fmma.poly_ord=2;
  fmma.Depth = 3;
  fmma.tree(target, source, ans);

  std::vector<TYPE> exact(tsize);
  for(int i=0; i<ssize; ++i){
    for(int j=0; j<tsize; ++j){
      TYPE len = 0.0;
      for(std::size_t dim=0; dim<DIM; ++dim){
        TYPE d = target[j][dim] - source[i][dim];
        len += d * d;
      }
      len = std::sqrt(len);
      exact[j] += 1.0/len;
    }
  }

  TYPE diff = 0.0;
  for(int i=0; i<tsize; ++i){
    TYPE d = ans[i]-exact[i];
    diff += d*d/exact[i]/exact[i];
  }
  diff = sqrt(diff / tsize);

  std::string comment = "<" + get_type<TYPE>() + "," + std::to_string(DIM) + ">(" + std::to_string(ssize) + ", " + std::to_string(tsize) + ", " + std::to_string(tol) + ")";

  if(diff < tol){
    pass(__FILE__, __func__, comment);
    return true;
  }else{
    fprintf(stderr, "diff = %e > tol = %e\n", diff, tol);
    failed(__FILE__, __func__, comment);
    return false;
  }
}
//}}}

//{{{ test_tree_1_r2
template<typename TYPE, std::size_t DIM>
bool test_tree_1_r2(int ssize, int tsize, TYPE tol){
  std::vector<std::array<TYPE, DIM>> source(ssize), target(tsize);
  std::vector<TYPE> source_weight(ssize);
  for(int i=0; i<ssize; ++i){
    for(std::size_t dim=0; dim<DIM; ++dim){
      source[i][dim] = ((TYPE)rand()/RAND_MAX-0.5)*5;
    }
    source_weight[i] = (TYPE)rand()/RAND_MAX;
  }
  for(int i=0; i<tsize; ++i){
    for(std::size_t dim=0; dim<DIM; ++dim){
      target[i][dim] = ((TYPE)rand()/RAND_MAX-0.5)*5;
    }
  }

  fmma::FMMA<TYPE, DIM> fmma;
  std::vector<TYPE> ans(tsize);
  auto fn = [](const std::array<TYPE, DIM>& target, const std::array<TYPE, DIM>& source){
    TYPE len2 = 0.0;
    TYPE sum = 0.0;
    for(std::size_t dim=0; dim<DIM; ++dim){
      TYPE diff = target[dim]-source[dim];
      len2 += diff*diff;
      sum += diff;
    }
    return sum/len2;
  };
  fmma.set_fn(fn);
  fmma.poly_ord=3;
  fmma.nrn_N = 6;
  fmma.tree(target, source_weight, source, ans);

  std::vector<TYPE> exact(tsize, 0.0);
  for(int i=0; i<ssize; ++i){
    for(int j=0; j<tsize; ++j){
      TYPE len2 = 0.0;
      TYPE sum = 0.0;
      for(std::size_t dim=0; dim<DIM; ++dim){
        TYPE d = target[j][dim] - source[i][dim];
        len2 += d * d;
        sum += d;
      }
      exact[j] += source_weight[i]*sum/len2;
    }
  }

  TYPE diff = 0.0;
  for(int i=0; i<tsize; ++i){
    TYPE d = ans[i]-exact[i];
    diff += d*d/exact[i]/exact[i];
  }
  diff = sqrt(diff / tsize);

  std::string comment = "<" + get_type<TYPE>() + "," + std::to_string(DIM) + ">(" + std::to_string(ssize) + ", " + std::to_string(tsize) + ", " + std::to_string(tol) + ")";

  if(diff < tol){
    pass(__FILE__, __func__, comment);
    return true;
  }else{
    fprintf(stderr, "diff = %e > tol = %e\n", diff, tol);
    failed(__FILE__, __func__, comment);
    return false;
  }
}
//}}}

//{{{ test_tree_sum2_r2
template<typename TYPE, std::size_t DIM>
bool test_tree_sum2_r2(int ssize, int tsize, TYPE tol){
  std::vector<std::array<TYPE, DIM>> source(ssize), target(tsize);
  std::vector<TYPE> source_weight(ssize);
  for(int i=0; i<ssize; ++i){
    for(std::size_t dim=0; dim<DIM; ++dim){
      source[i][dim] = ((TYPE)rand()/RAND_MAX-0.5)*5;
    }
    source_weight[i] = (TYPE)rand()/RAND_MAX;
  }
  for(int i=0; i<tsize; ++i){
    for(std::size_t dim=0; dim<DIM; ++dim){
      target[i][dim] = ((TYPE)rand()/RAND_MAX-0.5)*5;
    }
  }

  fmma::FMMA<TYPE, DIM> fmma;
  std::vector<TYPE> ans(tsize);
  auto fn = [](const std::array<TYPE, DIM>& target, const std::array<TYPE, DIM>& source){
    TYPE len2 = 0.0;
    TYPE sum2 = 0.0;
    for(std::size_t dim=0; dim<DIM; ++dim){
      TYPE diff = target[dim]-source[dim];
      TYPE sum = target[dim]+source[dim];
      len2 += diff*diff;
      sum2 += sum*sum;
    }
    return sum2/len2;
  };
  fmma.set_fn(fn);
  fmma.poly_ord=5;
  fmma.nrn_N = 6;
  fmma.tree(target, source_weight, source, ans);

  std::vector<TYPE> exact(tsize, 0.0);
  for(int i=0; i<ssize; ++i){
    for(int j=0; j<tsize; ++j){
      TYPE len2 = 0.0;
      TYPE sum2 = 0.0;
      for(std::size_t dim=0; dim<DIM; ++dim){
        TYPE d = target[j][dim] - source[i][dim];
        TYPE s = target[j][dim] + source[i][dim];
        len2 += d * d;
        sum2 += s * s;
      }
      exact[j] += source_weight[i]*sum2/len2;
    }
  }

  TYPE diff = 0.0;
  for(int i=0; i<tsize; ++i){
    TYPE d = ans[i]-exact[i];
    diff += d*d;
  }
  diff = sqrt(diff / tsize);

  std::string comment = "<" + get_type<TYPE>() + "," + std::to_string(DIM) + ">(" + std::to_string(ssize) + ", " + std::to_string(tsize) + ", " + std::to_string(tol) + ")";

  if(diff < tol){
    pass(__FILE__, __func__, comment);
    return true;
  }else{
    fprintf(stderr, "diff = %e > tol = %e\n", diff, tol);
    failed(__FILE__, __func__, comment);
    return false;
  }
}
//}}}



int main(void){
  srand(0);

  if(!test_tree<double, 1>(10, 10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_tree<double, 1>(10, 20, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_tree<double, 1>(20, 10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }

  if(!test_tree<double, 2>(10, 10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_tree<double, 2>(10, 20, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_tree<double, 2>(20, 10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }

  if(!test_tree_1_r2<double, 1>(10, 10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_tree_1_r2<double, 1>(10, 20, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_tree_1_r2<double, 1>(20, 10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }

  if(!test_tree_1_r2<double, 2>(10, 10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_tree_1_r2<double, 2>(10, 20, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_tree_1_r2<double, 2>(20, 10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }

  if(!test_tree_sum2_r2<double, 1>(10, 10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_tree_sum2_r2<double, 1>(10, 20, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_tree_sum2_r2<double, 1>(20, 10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }

  if(!test_tree_sum2_r2<double, 2>(10, 10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_tree_sum2_r2<double, 2>(10, 20, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_tree_sum2_r2<double, 2>(20, 10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }

  return 0;
}
