#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils/stat.h"

#define F(f, x, y) (f).v[((x) + (f).nx) % (f).nx + (((y) + (f).ny) % (f).ny) * (f).nx]

typedef struct
{
  double *v; // スピンの角度 [0, 1)
  int nx;
  int ny;
} field2d;

/**
 * @brief 2次元XYモデルのシミュレーションを行い、エネルギー、ヘリシティ係数、およびその自己相関時間を計算する
 *
 * @param beta 逆温度 (1/T)
 * @param energy_ptr 計算された平均エネルギーを格納するポインタ
 * @param helicity_ptr 計算された（熱平均の）ヘリシティ係数を格納するポインタ
 * @param autocorr_ptr ヘリシティ係数の自己相関時間を格納するポインタ
 * @param verbose 詳細なログを出力する場合は1、しない場合は0
 * @return int 成功時は0、失敗時は1
 */
int xy2d_mh(double beta, double *energy_ptr, double *helicity_ptr, double *autocorr_ptr, int verbose)
{
  field2d f;
  double energy, *energy_lst, *helicity_inst_lst, *Ex_lst, *Ix_lst;
  int n, mcs, nsamples, mcsequiv, mcsinterval, lag_max;
  FILE *fp, *fpa, *fpm;

  n = 32;
  nsamples = 10000;
  mcsequiv = 10000;
  mcsinterval = 10;
  lag_max = 1000;
  mcs = n * n;
  double d_theta = 0.2;

  if (verbose)
  {
    fp = fopen("xy2d_mh.dat", "w");
    fpa = fopen("xy2d_mh_autocorr.dat", "w");
    fpm = fopen("xy2d_mh_meta.dat", "w");

    fprintf(fpm, "Lattice size: %d x %d\n", n, n);
    fprintf(fpm, "Number of samples: %d\n", nsamples);
    fprintf(fpm, "Equilibration steps: %d\n", mcsequiv);
    fprintf(fpm, "MCS interval: %d\n", mcsinterval);
    fprintf(fpm, "Beta: %f\n", beta);
  }

  srand48(0);

  f.nx = n;
  f.ny = n;
  f.v = (double *)malloc(f.nx * f.ny * sizeof(double));
  energy_lst = (double *)malloc(nsamples * sizeof(double));
  helicity_inst_lst = (double *)malloc(nsamples * sizeof(double));
  Ex_lst = (double *)malloc(nsamples * sizeof(double));
  Ix_lst = (double *)malloc(nsamples * sizeof(double));

  for (int i = 0; i < f.nx * f.ny; i++)
    f.v[i] = drand48();

  energy = 0.0;
  for (int i = 0; i < f.nx; i++)
  {
    for (int j = 0; j < f.ny; j++)
    {
      double theta = F(f, i, j);
      double down = F(f, i, j + 1);
      double right = F(f, i + 1, j);
      energy += cos(2.0 * M_PI * (theta - down)) + cos(2.0 * M_PI * (theta - right));
    }
  }
  energy = -energy / (double)(f.nx * f.ny);

  for (int i = 0; i < mcsequiv + mcsinterval * nsamples; i++)
  {
    for (int j = 0; j < mcs; j++)
    {
      int x = lrand48() % f.nx, y = lrand48() % f.ny;
      double old_theta = F(f, x, y);
      double new_theta = old_theta + (drand48() - 0.5) * d_theta;
      new_theta -= floor(new_theta);

      double up = F(f, x, y + 1), down = F(f, x, y - 1);
      double left = F(f, x - 1, y), right = F(f, x + 1, y);

      double d_cos_sum = (cos(2.0 * M_PI * (new_theta - up)) - cos(2.0 * M_PI * (old_theta - up))) + (cos(2.0 * M_PI * (new_theta - down)) - cos(2.0 * M_PI * (old_theta - down))) + (cos(2.0 * M_PI * (new_theta - left)) - cos(2.0 * M_PI * (old_theta - left))) + (cos(2.0 * M_PI * (new_theta - right)) - cos(2.0 * M_PI * (old_theta - right)));
      double deltaE = -d_cos_sum;

      if (deltaE < 0 || drand48() < exp(-beta * deltaE))
      {
        F(f, x, y) = new_theta;
        energy += deltaE / (double)(f.nx * f.ny);
      }
    }

    if (i >= mcsequiv && (i - mcsequiv) % mcsinterval == 0)
    {
      int sample_index = (i - mcsequiv) / mcsinterval;
      if (sample_index < nsamples)
      {
        energy_lst[sample_index] = energy;

        double Ex_current = 0.0, Ix_current = 0.0;
        for (int k_y = 0; k_y < f.ny; ++k_y)
        {
          for (int k_x = 0; k_x < f.nx; ++k_x)
          {
            double d_theta = F(f, k_x, k_y) - F(f, k_x + 1, k_y);
            Ex_current += cos(2.0 * M_PI * d_theta);
            Ix_current += sin(2.0 * M_PI * d_theta);
          }
        }
        Ex_lst[sample_index] = Ex_current;
        Ix_lst[sample_index] = Ix_current;

        // 「瞬時ヘリシティ係数」を計算して時系列リストに保存
        double L2 = (double)(f.nx * f.ny);
        helicity_inst_lst[sample_index] = (Ex_current / L2) - beta * (Ix_current * Ix_current / L2);

        if (verbose)
          fprintf(fp, "%d %f %f\n", sample_index, energy_lst[sample_index], helicity_inst_lst[sample_index]);
      }
    }
  }

  double mean_energy = 0.0, mean_Ex = 0.0, mean_Ix_sq = 0.0, mean_helicity_inst = 0.0;
  for (int i = 0; i < nsamples; i++)
  {
    mean_energy += energy_lst[i];
    mean_Ex += Ex_lst[i];
    mean_Ix_sq += Ix_lst[i] * Ix_lst[i];
    mean_helicity_inst += helicity_inst_lst[i];
  }
  mean_energy /= nsamples;
  mean_Ex /= nsamples;
  mean_Ix_sq /= nsamples;
  mean_helicity_inst /= nsamples;

  double *autocorr_helicity = (double *)malloc(lag_max * sizeof(double));
  for (int lag = 0; lag < lag_max; lag++)
  {
    autocorr_helicity[lag] = autocorrelation(helicity_inst_lst, mean_helicity_inst, nsamples, lag);
    if (verbose)
      fprintf(fpa, "%d %f\n", lag * mcsinterval, autocorr_helicity[lag]);
  }

  *energy_ptr = mean_energy;
  double L2 = (double)(f.nx * f.ny);
  *helicity_ptr = (mean_Ex / L2) - beta * (mean_Ix_sq / L2);
  *autocorr_ptr = autocorrelation_time(autocorr_helicity, lag_max);
  if (*autocorr_ptr < 0)
    *autocorr_ptr = 0;

  if (verbose)
  {
    fprintf(fpm, "Mean Energy: %f\n", *energy_ptr);
    fprintf(fpm, "Helicity Modulus (from mean values): %f\n", *helicity_ptr);
    fprintf(fpm, "Mean of instantaneous Helicity Modulus: %f\n", mean_helicity_inst);
    fprintf(fpm, "Autocorrelation time of Helicity Modulus: %f\n", *autocorr_ptr);
  }

  free(f.v);
  free(energy_lst);
  free(helicity_inst_lst);
  free(Ex_lst);
  free(Ix_lst);
  free(autocorr_helicity);

  if (verbose)
  {
    fclose(fp);
    fclose(fpa);
    fclose(fpm);
  }

  return 0;
}
