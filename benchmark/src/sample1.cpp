#include"../../include/fmma/fmma.hpp"
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<vector>
#include<array>
#include<chrono>

int main(void){

  int N = 1024;
  int REPEAT = 5;
  int ORDER = 3;

  auto fn = [](const std::array<double, 1>& source, const std::array<double, 1>& target){
    double diff = source[0]-target[0];
    double len = diff*diff;
    return 1.0/std::sqrt(len);
  };

  std::vector<double> exact_time(REPEAT);
  std::vector<int> size(REPEAT);
  std::vector<std::vector<double>> nrnmm_time(REPEAT, std::vector<double>(ORDER));
  std::vector<std::vector<double>> nrnmm_error(REPEAT, std::vector<double>(ORDER));

  std::chrono::system_clock::time_point start, end;

  for(int repeat=0; repeat<REPEAT; ++repeat){
    N *= 2;

    size[repeat] = N;
    fmma::FMMA<double, 1> fmma;
    std::vector<std::array<double, 1>> source(N), target(N);
    for(int i=0; i<N; ++i){
      source[i][0] = (double)rand()/RAND_MAX;
      target[i][0] = (double)rand()/RAND_MAX;
    }
    std::vector<double> ans_exact(N);
    fmma.fn = fn;

    start = std::chrono::system_clock::now();
    fmma.exact(source, target, ans_exact);
    end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    exact_time[repeat] = elapsed;

    for(int order=1; order<=ORDER; ++order){
      std::vector<double> ans_nrnmm(N);

      fmma.poly_ord = order;
      start = std::chrono::system_clock::now();
      fmma.nrnmm(source, target, ans_nrnmm);
      end = std::chrono::system_clock::now();
      elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
      nrnmm_time[repeat][order-1] = elapsed;
    }
  } 

  {
    FILE *fp;
    fp = fopen("time.csv", "w");
    fprintf(fp, "N, ");
    fprintf(fp, "exact");
    for(int order=1; order<=ORDER; ++order){
      fprintf(fp, ", nrnmm(%d)", order);
    }
    fprintf(fp, "\n");

    for(int repeat=0; repeat<REPEAT; ++repeat){
      fprintf(fp, "%d", size[repeat]);
      fprintf(fp, ", %lf", exact_time[repeat]);
      for(int order=1; order<=ORDER; ++order){
        fprintf(fp, ", %lf", nrnmm_time[repeat][order-1]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  return 0;
}
