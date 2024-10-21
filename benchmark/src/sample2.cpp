#include"../../include/fmma/fmma.hpp"
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<vector>
#include<array>
#include<chrono>

int main(void){

  int N = 1024;
  int REPEAT = 7;
  int ORDER = 3;

  auto fn = [](const std::array<double, 2>& target, const std::array<double, 2>& source){
    double diff0 = target[0]-source[0];
    double diff1 = target[1]-source[1];
    double len = diff0*diff0 + diff1*diff1;
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
    fmma::FMMA<double, 2> fmma;
    std::vector<std::array<double, 2>> source(N), target(N);
    for(int n=0; n<N; ++n){
      for(int dim=0; dim<2; ++dim){
        source[n][dim] = (double)rand()/RAND_MAX;
        target[n][dim] = (double)rand()/RAND_MAX;
      }
    }
    std::vector<double> ans_exact(N);
    fmma.set_fn(fn);

    start = std::chrono::system_clock::now();
    fprintf(stderr, "exact calculation\n");
    fmma.exact(source, target, ans_exact);
    end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    exact_time[repeat] = elapsed;

    fprintf(stderr, "nrnmm calculation\n");
    for(int order=1; order<=ORDER; ++order){
      std::vector<double> ans_nrnmm(N);

      fmma.poly_ord = order;
      start = std::chrono::system_clock::now();
      fmma.nrnmm(source, target, ans_nrnmm);
      end = std::chrono::system_clock::now();
      elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
      nrnmm_time[repeat][order-1] = elapsed;

      double diff = 0.0;
      for(int n=0; n<N; ++n){
        double tmp = (ans_exact[n]-ans_nrnmm[n])/ans_exact[n];
        diff += tmp*tmp;
      }
      diff = sqrt(diff/N);
      nrnmm_error[repeat][order-1] = diff;
    }
  }

  {
    FILE *fp;
    fp = fopen("time_2.csv", "w");
    fprintf(fp, "N");
    fprintf(fp, ", exact");
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

  {
    FILE *fp;
    fp = fopen("error_2.csv", "w");
    fprintf(fp, "N");
    for(int order=1; order<=ORDER; ++order){
      fprintf(fp, ", nrnmm(%d)", order);
    }
    fprintf(fp, "\n");

    for(int repeat=0; repeat<REPEAT; ++repeat){
      fprintf(fp, "%d", size[repeat]);
      for(int order=1; order<=ORDER; ++order){
        fprintf(fp, ", %e", nrnmm_error[repeat][order-1]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  {
    FILE *fp;
    fp = fopen("comment_2.md", "w");
    fprintf(fp, "2D results\n");
    fprintf(fp, "| type | N | time[ms] | relative error |\n");
    fprintf(fp, "| --- | --- | --- | --- |\n");
    fprintf(fp, "| exact | %d | %e | --- |\n", size[0], exact_time[0]);
    for(int repeat=1; repeat<REPEAT; ++repeat){
      fprintf(fp, "| | %d | %e | --- |\n", size[repeat], exact_time[repeat]);
    }
    for(int order=1; order<=ORDER; ++order){
      fprintf(fp, "| nrnmm(%d) | %d | %e | %e |\n", order, size[0], nrnmm_time[0][order-1], nrnmm_error[0][order-1]);
      for(int repeat=1; repeat<REPEAT; ++repeat){
        fprintf(fp, "| | %d | %e | %e |\n", size[repeat], nrnmm_time[repeat][order-1], nrnmm_error[repeat][order-1]);
      }
    }
    fclose(fp);
  }

  return 0;
}
