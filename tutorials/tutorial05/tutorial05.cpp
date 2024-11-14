#include<fmma/fmma.hpp>
#include<cstdlib>
#include<cstdio>
#include<array>
#include<vector>

int main(void){

  fmma::FMMA<double, 2> fmma;
  std::size_t N = 10;
  std::vector<std::array<double, 2>> source(N), target(N);
  std::vector<double> source_weight(N);
  std::vector<double> ans(N);
  for(std::size_t i=0; i<N; ++i){
    source[i][0] = (double)i/N;
    source[i][1] = (double)i/N;
    target[i][0] = ((double)i+0.5)/N;
    target[i][1] = ((double)i+0.5)/N;
    source_weight[i] = (double)i/N;
  }

  fmma.set_fn([](const std::array<double, 2>& x, const std::array<double, 2>& y){
      double len = std::sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]));
      return 1.0/len;
      });
  fmma.set_solver_type("fmm");
  fmma.set_Depth(3);
  fmma.set_poly_ord(3);
  fmma.solve(target, source_weight, source, ans);

  fprintf(stderr, "ans :\n");
  for(std::size_t i=0; i<N; ++i){
    fprintf(stderr, "%lf ", ans[i]);
  }
  fprintf(stderr, "\n");

  return 0;
}
