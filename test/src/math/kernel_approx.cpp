#include"../../../include/fmma/math.hpp"
#include"../test_common.hpp"
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<vector>
#include<array>
#include<functional>

template<typename TYPE>
TYPE predict(std::function<TYPE(TYPE, TYPE)>& fn, TYPE x, TYPE y, int n, TYPE ymin=-1.0, TYPE ymax=1.0){
  std::vector<TYPE> chebyshev_node(n+1);
  for(int k=0; k<=n; ++k){
    chebyshev_node[k] = cos((2.0*k+1.0)/(2*n+2)*M_PI);
  }

  TYPE ans = 0.0;
  for(int k=0; k<=n; ++k){
    TYPE chebyshev_node_original = chebyshev_node[k]*(ymax-ymin)/2.0 + (ymax+ymin)/2.0;
    TYPE y_new = (2.0*y-(ymax+ymin))/(ymax-ymin);
    ans += fn(x, chebyshev_node_original)*fmma::SChebyshev(n+1, chebyshev_node[k], y_new);
  }
  return ans;
}

template<typename TYPE>
TYPE predict(std::function<TYPE(const std::array<TYPE, 2>&, const std::array<TYPE, 2>&)>& fn, const std::array<TYPE, 2>& x, const std::array<TYPE, 2>& y, int n, const std::array<TYPE, 2>& ymin, const std::array<TYPE, 2>& ymax){
  std::vector<TYPE> chebyshev_node(n+1);
  for(int k=0; k<=n; ++k){
    chebyshev_node[k] = cos((2.0*k+1.0)/(2*n+2)*M_PI);
  }

  TYPE ans = 0.0;
  for(int k0=0; k0<=n; ++k0){
    for(int k1=0; k1<=n; ++k1){
      TYPE chebyshev_node_original0 = chebyshev_node[k0]*(ymax[0]-ymin[0])/2.0 + (ymax[0]+ymin[0])/2.0;
      TYPE chebyshev_node_original1 = chebyshev_node[k1]*(ymax[1]-ymin[1])/2.0 + (ymax[1]+ymin[1])/2.0;
      TYPE y0_new = (2.0*y[0]-(ymax[0]+ymin[0]))/(ymax[0]-ymin[0]);
      TYPE y1_new = (2.0*y[1]-(ymax[1]+ymin[1]))/(ymax[1]-ymin[1]);
      std::array<TYPE, 2> chebyshev_node_original = {chebyshev_node_original0, chebyshev_node_original1};
      ans += fn(x, chebyshev_node_original)*fmma::SChebyshev(n+1, chebyshev_node[k0], y0_new)*fmma::SChebyshev(n+1, chebyshev_node[k1], y1_new);
    }
  }
  return ans;
}

//{{{ test_kernel_approx_1
template<typename TYPE>
bool test_kernel_approx_1(int Nmax, TYPE tol){
  std::vector<TYPE> diff(Nmax, 0.0);

  std::function<TYPE(TYPE x, TYPE y)> fn = [](TYPE x, TYPE y){
    return 1.0/((x-y)*(x-y));
  };

  for(std::size_t d=0; d<10; ++d){
    TYPE x = (TYPE)rand()*2.0/RAND_MAX-4.0;
    TYPE y = (TYPE)rand()*2.0/RAND_MAX-1.0;
    for(int n=0; n<Nmax; ++n){
      TYPE ans = predict(fn, x, y, n);
      TYPE exact = fn(x, y);

      diff[n] += (ans-exact)*(ans-exact);
    }
  }

  bool flag = true;
  for(int i=0; i+1<Nmax; ++i){
    if(diff[i] < diff[i+1] && diff[i+1] > tol){
      fprintf(stderr, "diff [%d] %e < diff [%d] %e\n", i, diff[i], i+1, diff[i+1]);
      flag = false;
    }
  }
  if(flag && diff[Nmax-1]<tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    for(int i=0; i<Nmax; ++i){
      fprintf(stderr, "%d-th diff = %e\n", i, diff[i]);
    }
    fprintf(stderr, "%d (tol = %e)\n", flag, tol);
    failed(__FILE__, __func__);
    return false;
  }
}
//}}}

//{{{ test_kernel_approx_2
template<typename TYPE>
bool test_kernel_approx_2(int Nmax, TYPE tol){
  std::vector<TYPE> diff(Nmax, 0.0);

  std::function<TYPE(TYPE x, TYPE y)> fn = [](TYPE x, TYPE y){
    return log(16.0*(x-y)*(x-y)) + exp((x-y)*(x-y)/8.0);
  };

  for(std::size_t d=0; d<10; ++d){
    TYPE x = (TYPE)rand()/RAND_MAX-4.0;
    TYPE y = (TYPE)rand()*4.0/RAND_MAX+1.0;
    TYPE ymin = 1.0, ymax = 5.0;
    TYPE exact = fn(x, y);

    for(int n=0; n<Nmax; ++n){
      TYPE ans = predict(fn, x, y, n, ymin, ymax);

      diff[n] += (ans-exact)*(ans-exact);
    }
  }

  bool flag = true;
  for(int i=0; i+1<Nmax; ++i){
    if(diff[i] < diff[i+1] && diff[+1] > tol){
      fprintf(stderr, "diff [%d] %e < diff [%d] %e\n", i, diff[i], i+1, diff[i+1]);
      flag = false;
    }
  }
  if(flag && diff[Nmax-1]<tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    for(int i=0; i<Nmax; ++i){
      fprintf(stderr, "%d-th diff = %e\n", i, diff[i]);
    }
    fprintf(stderr, "%d (tol = %e)\n", flag, tol);
    failed(__FILE__, __func__);
    return false;
  }
}
//}}}

//{{{ test_kernel_approx_3
template<typename TYPE>
bool test_kernel_approx_3(int Nmax, TYPE tol){
  std::vector<TYPE> diff(Nmax, 0.0);

  std::function<TYPE(const std::array<TYPE, 2>& x, const std::array<TYPE, 2>& y)> fn = [](const std::array<TYPE, 2>& x, const std::array<TYPE, 2>& y){
    return (x[0]-y[0])*(x[0]-y[0])-4*(x[1]-y[1])/(x[0]-y[0]);
  };

  for(std::size_t d=0; d<10; ++d){
    TYPE x0 = (TYPE)rand()*2.0/RAND_MAX-4.0;
    TYPE x1 = (TYPE)rand()*2.0/RAND_MAX-4.0;
    TYPE y0 = (TYPE)rand()*2.0/RAND_MAX-1.0;
    TYPE y1 = (TYPE)rand()*2.0/RAND_MAX-1.0;
    std::array<TYPE, 2> x{x0, x1}, y{y0, y1};
    std::array<TYPE, 2> ymin{-1.0, -1.0}, ymax{1.0, 1.0};
    for(int n=0; n<Nmax; ++n){
      TYPE ans = predict(fn, x, y, n, ymin, ymax);
      TYPE exact = fn(x, y);

      diff[n] += (ans-exact)*(ans-exact);
    }
  }

  bool flag = true;
  for(int i=0; i+1<Nmax; ++i){
    if(diff[i] < diff[i+1] && diff[i+1] > tol){
      fprintf(stderr, "diff [%d] %e < diff [%d] %e\n", i, diff[i], i+1, diff[i+1]);
      flag = false;
    }
  }
  if(flag && diff[Nmax-1]<tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    for(int i=0; i<Nmax; ++i){
      fprintf(stderr, "%d-th diff = %e\n", i, diff[i]);
    }
    fprintf(stderr, "%d (tol = %e)\n", flag, tol);
    failed(__FILE__, __func__);
    return false;
  }
}
//}}}

//{{{ test_kernel_approx_4
template<typename TYPE>
bool test_kernel_approx_4(int Nmax, TYPE tol){
  std::vector<TYPE> diff(Nmax, 0.0);

  std::function<TYPE(const std::array<TYPE, 2>& x, const std::array<TYPE, 2>& y)> fn = [](const std::array<TYPE, 2>& x, const std::array<TYPE, 2>& y){
    TYPE len2 = (x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]);
    return 1.0/len2;
  };

  for(std::size_t d=0; d<10; ++d){
    TYPE x0 = (TYPE)rand()/RAND_MAX-3.0;
    TYPE x1 = (TYPE)rand()/RAND_MAX-3.0;
    TYPE y0 = (TYPE)rand()*4.0/RAND_MAX+1.0;
    TYPE y1 = (TYPE)rand()*4.0/RAND_MAX+1.0;
    std::array<TYPE, 2> x{x0, x1}, y{y0, y1};
    std::array<TYPE, 2> ymin{1.0, 1.0}, ymax{5.0, 5.0};
    for(int n=0; n<Nmax; ++n){
      TYPE ans = predict(fn, x, y, n, ymin, ymax);
      TYPE exact = fn(x, y);

      diff[n] += (ans-exact)*(ans-exact);
    }
  }

  bool flag = true;
  for(int i=0; i+1<Nmax; ++i){
    if(diff[i] < diff[i+1] && diff[i+1] > tol){
      fprintf(stderr, "diff [%d] %e < diff [%d] %e\n", i, diff[i], i+1, diff[i+1]);
      flag = false;
    }
  }
  if(flag && diff[Nmax-1]<tol){
    pass(__FILE__, __func__);
    return true;
  }else{
    for(int i=0; i<Nmax; ++i){
      fprintf(stderr, "%d-th diff = %e\n", i, diff[i]);
    }
    fprintf(stderr, "%d (tol = %e)\n", flag, tol);
    failed(__FILE__, __func__);
    return false;
  }
}
//}}}

int main(void){
  srand(0);

  if(!test_kernel_approx_1<double>(10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_kernel_approx_2<double>(20, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_kernel_approx_3<double>(10, 1.0e-6)){
    exit(EXIT_FAILURE);
  }
  if(!test_kernel_approx_4<double>(20, 1.0e-6)){
    exit(EXIT_FAILURE);
  }

  if(!test_kernel_approx_1<float>(10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_kernel_approx_2<float>(13, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_kernel_approx_3<float>(10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }
  if(!test_kernel_approx_4<float>(10, 1.0e-3)){
    exit(EXIT_FAILURE);
  }

  return 0;
}
