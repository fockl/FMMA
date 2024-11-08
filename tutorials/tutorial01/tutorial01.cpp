#include<fmma/fmma.hpp>
#include<cstdlib>
#include<cstdio>
#include<array>
#include<vector>

int main(void){

  fmma::FMMA<double, 1> fmma;
  std::size_t N = 10;
  std::vector<std::array<double, 1>> source(N), target(N);
  std::vector<double> ans(N);
  for(std::size_t i=0; i<N; ++i){
    source[i][0] = (double)i/N;
    target[i][0] = ((double)i+0.5)/N;
  }

  fmma.set_fn([](const std::array<double, 1>& x, const std::array<double, 1>& y){
      return 1.0/((x[0]-y[0])*(x[0]-y[0]));
      });
  fmma.solve_type = "exact";
  fmma.solve(target, source, ans);

  fprintf(stderr, "ans :\n");
  for(std::size_t i=0; i<N; ++i){
    fprintf(stderr, "%lf ", ans[i]);
  }
  fprintf(stderr, "\n");

  return 0;
}
