#include"../../../include/fmma/chebyshev_approx.hpp"
#include"../test_common.hpp"
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<vector>
#include<array>

//{{{ test_coef_1
template<typename TYPE>
bool test_coef_1(int n, int Nth, TYPE tol){
  std::function<TYPE(const std::array<TYPE, 1>& pos)> fn = [n](const std::array<TYPE, 1>& pos){
    if(n==0){
      return 1.0;
    }else if(n==1){
      return pos[0];
    }else if(n==2){
      return 2*pos[0]*pos[0]-1;
    }else if(n==3){
      return 4*pos[0]*pos[0]*pos[0]-3*pos[0];
    }else{
      fprintf(stderr, "%s:%d ERROR : test_coef_1 test is for n<=3 but %d\n", __FILE__, __LINE__, n);
      exit(EXIT_FAILURE);
    }
  };

  fmma::CHEBYSHEV_APPROX<TYPE, 1> che;
  che.fn = fn;
  che.Nth = Nth;
  che.initialize();

  std::vector<TYPE> coef(Nth+1);
  coef[n] = 1.0;

  TYPE diff = 0.0;
  for(int i=0; i<=Nth; ++i){
    TYPE d = coef[i] - che.coef[i];
    diff += d*d;
  }

  if(diff < tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    fprintf(stderr, "diff = %e > rol = %e\n", diff, tol);
    failed(__FILE__, __func__);
    fprintf(stderr, "coef vs che.coef\n");
    for(int i=0; i<=Nth; ++i){
      fprintf(stderr, "%lf %lf\n", coef[i], che.coef[i]);
    }
    return false;
  }
}
//}}}

//{{{ test_coef_2
template<typename TYPE>
bool test_coef_2(int n, int m, int Nth, TYPE tol){
  std::function<TYPE(const std::array<TYPE, 2>& pos)> fn = [n, m](const std::array<TYPE, 2>& pos){
    TYPE val = 1.0;
    if(n==0){
      val *= 1.0;
    }else if(n==1){
      val *= pos[0];
    }else if(n==2){
      val *= 2*pos[0]*pos[0]-1;
    }else if(n==3){
      val *= 4*pos[0]*pos[0]*pos[0]-3*pos[0];
    }else{
      fprintf(stderr, "%s:%d ERROR : test_coef_2 test is for n<=3 but %d\n", __FILE__, __LINE__, n);
      exit(EXIT_FAILURE);
    }
    if(m==0){
      val *= 1.0;
    }else if(m==1){
      val *= pos[1];
    }else if(m==2){
      val *= 2*pos[1]*pos[1]-1;
    }else if(m==3){
      val *= 4*pos[1]*pos[1]*pos[1]-3*pos[1];
    }else{
      fprintf(stderr, "%s:%d ERROR : test_coef_2 test is for n<=3 but %d\n", __FILE__, __LINE__, m);
      exit(EXIT_FAILURE);
    }
    return val;
  };

  fmma::CHEBYSHEV_APPROX<TYPE, 2> che;
  che.fn = fn;
  che.Nth = Nth;
  che.initialize();

  std::vector<TYPE> coef((Nth+1)*(Nth+1));
  coef[m*(Nth+1)+n] = 1.0;

  TYPE diff = 0.0;
  for(int i=0; i<(Nth+1)*(Nth+1); ++i){
    TYPE d = coef[i] - che.coef[i];
    diff += d*d;
  }

  if(diff < tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    fprintf(stderr, "diff = %e > rol = %e\n", diff, tol);
    failed(__FILE__, __func__);
    fprintf(stderr, "coef vs che.coef\n");
    for(int i=0; i<(Nth+1)*(Nth+1); ++i){
      fprintf(stderr, "%lf %lf\n", coef[i], che.coef[i]);
    }
    return false;
  }
}
//}}}

//{{{ test_predict_1
template<typename TYPE>
bool test_predict_1(int n, int Nth, TYPE tol){
  std::function<TYPE(const std::array<TYPE, 1>& pos)> fn = [n](const std::array<TYPE, 1>& pos){
    if(n==0){
      return 1.0;
    }else if(n==1){
      return pos[0];
    }else if(n==2){
      return 2*pos[0]*pos[0]-1;
    }else if(n==3){
      return 4*pos[0]*pos[0]*pos[0]-3*pos[0];
    }else{
      fprintf(stderr, "%s:%d ERROR : test_predict_1 test is for n<=3 but %d\n", __FILE__, __LINE__, n);
      exit(EXIT_FAILURE);
    }
  };

  fmma::CHEBYSHEV_APPROX<TYPE, 1> che;
  che.fn = fn;
  che.Nth = Nth;
  che.initialize();

  TYPE diff = 0.0;
  for(std::size_t d=0; d<10; ++d){
    std::array<TYPE, 1> pos;
    pos[0] = (TYPE)rand()*2.0/RAND_MAX-1.0;

    TYPE ans = che.predict(pos);
    TYPE exact = fn(pos);

    diff += (ans-exact)*(ans-exact);
  }

  if(diff < tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    fprintf(stderr, "diff = %e > rol = %e\n", diff, tol);
    failed(__FILE__, __func__);
    return false;
  }
}
//}}}

//{{{ test_predict_2
template<typename TYPE>
bool test_predict_2(int n, int m, int Nth, TYPE tol){
  std::function<TYPE(const std::array<TYPE, 2>& pos)> fn = [n, m](const std::array<TYPE, 2>& pos){
    TYPE val = 1.0;
    if(n==0){
      val *= 1.0;
    }else if(n==1){
      val *= pos[0];
    }else if(n==2){
      val *= 2*pos[0]*pos[0]-1;
    }else if(n==3){
      val *= 4*pos[0]*pos[0]*pos[0]-3*pos[0];
    }else{
      fprintf(stderr, "%s:%d ERROR : test_coef_2 test is for n<=3 but %d\n", __FILE__, __LINE__, n);
      exit(EXIT_FAILURE);
    }
    if(m==0){
      val *= 1.0;
    }else if(m==1){
      val *= pos[1];
    }else if(m==2){
      val *= 2*pos[1]*pos[1]-1;
    }else if(m==3){
      val *= 4*pos[1]*pos[1]*pos[1]-3*pos[1];
    }else{
      fprintf(stderr, "%s:%d ERROR : test_coef_2 test is for n<=3 but %d\n", __FILE__, __LINE__, m);
      exit(EXIT_FAILURE);
    }
    return val;
  };

  fmma::CHEBYSHEV_APPROX<TYPE, 2> che;
  che.fn = fn;
  che.Nth = Nth;
  che.initialize();

  TYPE diff = 0.0;
  for(std::size_t d=0; d<10; ++d){
    std::array<TYPE, 2> pos;
    pos[0] = (TYPE)rand()*2.0/RAND_MAX-1.0;
    pos[1] = (TYPE)rand()*2.0/RAND_MAX-1.0;

    TYPE ans = che.predict(pos);
    TYPE exact = fn(pos);

    diff += (ans-exact)*(ans-exact);
  }

  if(diff < tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    fprintf(stderr, "diff = %e > rol = %e\n", diff, tol);
    failed(__FILE__, __func__);
    return false;
  }
}
//}}}

//{{{ test_predict_3
template<typename TYPE>
bool test_predict_3(int Nth, TYPE tol){
  std::function<TYPE(const std::array<TYPE, 3>& pos)> fn = [](const std::array<TYPE, 3>& pos){
    TYPE len2 = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];

    return 1.0/(1.0+len2);
  };

  fmma::CHEBYSHEV_APPROX<TYPE, 3> che;
  che.fn = fn;
  che.Nth = Nth;
  che.initialize();

  TYPE diff = 0.0;
  for(std::size_t d=0; d<10; ++d){
    std::array<TYPE, 3> pos;
    pos[0] = (TYPE)rand()*2.0/RAND_MAX-1.0;
    pos[1] = (TYPE)rand()*2.0/RAND_MAX-1.0;
    pos[2] = (TYPE)rand()*2.0/RAND_MAX-1.0;

    TYPE ans = che.predict(pos);
    TYPE exact = fn(pos);

    fprintf(stderr, "%lf vs %lf\n", ans, exact);

    diff += (ans-exact)*(ans-exact);
  }

  if(diff < tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    fprintf(stderr, "diff = %e > rol = %e\n", diff, tol);
    failed(__FILE__, __func__);
    return false;
  }
}
//}}}


int main(void){
  srand(0);

  if(!test_coef_1<double>(0, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_coef_1<double>(1, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_coef_1<double>(2, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_coef_1<double>(3, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_coef_2<double>(0, 0, 5, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_coef_2<double>(0, 1, 5, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_coef_2<double>(2, 0, 5, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_coef_2<double>(1, 3, 5, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_predict_1<double>(0, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_predict_1<double>(1, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_predict_1<double>(2, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_predict_1<double>(3, 10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_predict_2<double>(0, 0, 4, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_predict_2<double>(1, 0, 4, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_predict_2<double>(0, 2, 4, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_predict_2<double>(3, 2, 4, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_predict_3<double>(6, 1.0e-5)){
    exit(EXIT_FAILURE);
  }

  
  return 0;
}
