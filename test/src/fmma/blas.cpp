#include"../../../include/fmma/fmma.hpp"
#include"../test_common.hpp"
#include<cstdlib>
#include<string>

int main(int argc, char** argv){
  if(argc==1){
    fprintf(stderr, "Usage: %s blas_flag\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  bool blas_flag = std::stoi(argv[1]);

  fmma::FMMA<double, 1> fmma;

  if(blas_flag != fmma.check_blas()){
    fprintf(stderr, "blas_flag = %d != check_blas() = %d\n", blas_flag, fmma.check_blas());
    failed(__FILE__, __func__);
    exit(EXIT_FAILURE);
  }else{
    pass(__FILE__, __func__);
  }

  return 0;
}
