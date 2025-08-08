#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils/stat.h"

#define F(f, x, y) (f).v[((x) + (f).nx) % (f).nx + (((y) + (f).ny) % (f).ny) * (f).nx]

typedef struct
{
  int *v;
  int nx;
  int ny;
} field2d;

int main()
{
  field2d f;
  double beta;
  double energy, *energy_lst, *mag_lst;
  int n, mcs, nsamples, mcsequiv, mcsinterval, lag_max;
  FILE *fp, *fpa, *fpm;

  n = 32;
  nsamples = 10000;
  mcsequiv = 10000;
  mcsinterval = 1;
  lag_max = 1000;

  beta = 0.44;
  mcs = n * n;

  fp = fopen("ising2d_mh.dat", "w");
  fpa = fopen("ising2d_mh_autocorr.dat", "w");
  fpm = fopen("ising2d_mh_meta.dat", "w");

  fprintf(fpm, "Lattice size: %d x %d\n", n, n);
  fprintf(fpm, "Monte Carlo steps per sweep: %d\n", mcs);
  fprintf(fpm, "Number of samples: %d\n", nsamples);
  fprintf(fpm, "Equilibration steps: %d\n", mcsequiv);
  fprintf(fpm, "MCS interval: %d\n", mcsinterval);
  fprintf(fpm, "Beta: %f\n", beta);

  srand48(0);

  f.nx = n;
  f.ny = n;
  f.v = (int *)malloc(f.nx * f.ny * sizeof(int));
  if (f.v == NULL)
  {
    fprintf(stderr, "Memory allocation failed\n");
    return 1;
  }

  for (int i = 0; i < f.nx * f.ny; i++)
    f.v[i] = (lrand48() % 2) * 2 - 1;

  energy_lst = (double *)malloc(nsamples * sizeof(double));
  mag_lst = (double *)malloc(nsamples * sizeof(double));
  if (energy_lst == NULL || mag_lst == NULL)
  {
    fprintf(stderr, "Memory allocation failed for energy or magnetization lists\n");
    free(f.v);
    return 1;
  }
  for (int i = 0; i < nsamples; i++)
  {
    energy_lst[i] = 0.0;
    mag_lst[i] = 0.0;
  }

  energy = 0.0;
  for (int i = 0; i < f.nx; i++)
  {
    for (int j = 0; j < f.ny; j++)
    {
      int spin = F(f, i, j);
      int down = F(f, i, j + 1);
      int right = F(f, i + 1, j);
      energy += -spin * (down + right);
    }
  }
  energy /= 2.0;
  energy /= (double)(f.nx * f.ny);

  for (int i = 0; i < mcsequiv + mcsinterval * nsamples; i++)
  {
    for (int j = 0; j < mcs; j++)
    {
      int x = lrand48() % f.nx;
      int y = lrand48() % f.ny;
      int spin = F(f, x, y);
      int up = F(f, x, y + 1);
      int down = F(f, x, y - 1);
      int left = F(f, x - 1, y);
      int right = F(f, x + 1, y);
      double deltaE = 2.0 * spin * (up + down + left + right);
      if (deltaE < 0 || drand48() < exp(-beta * deltaE))
      {
        f.v[((x) + (f.nx)) % (f.nx) + (((y) + (f.ny)) % (f.ny)) * (f.nx)] = -spin;
        energy += deltaE / (double)(f.nx * f.ny);
      }
    }
    if (i >= mcsequiv && (i - mcsequiv) % mcsinterval == 0)
    {
      int sample_index = (i - mcsequiv) / mcsinterval;
      energy_lst[sample_index] = energy;
      double magnetization = 0.0;
      for (int k = 0; k < f.nx * f.ny; k++)
        magnetization += f.v[k];
      mag_lst[sample_index] = fabs(magnetization / (f.nx * f.ny));
      fprintf(fp, "%d %f %f\n", sample_index,
              energy_lst[sample_index], mag_lst[sample_index]);
    }
  }

  double mean_energy = 0.0, mean_mag = 0.0;
  for (int i = 0; i < nsamples; i++)
  {
    mean_energy += energy_lst[i];
    mean_mag += mag_lst[i];
  }
  mean_energy /= nsamples;
  mean_mag /= nsamples;
  fprintf(fpm, "Mean Energy: %f\n", mean_energy);
  fprintf(fpm, "Mean Magnetization: %f\n", mean_mag);

  double *autocorr_mag;
  autocorr_mag = (double *)malloc(lag_max * sizeof(double));
  if (autocorr_mag == NULL)
  {
    fprintf(stderr, "Memory allocation failed for autocorrelation array\n");
  }
  for (int lag = 0; lag < lag_max; lag++)
  {
    autocorr_mag[lag] = autocorrelation(mag_lst, mean_mag, nsamples, lag);
    fprintf(fpa, "%d %f\n", lag * mcsinterval, autocorr_mag[lag]);
  }

  free(f.v);
  free(energy_lst);
  free(mag_lst);
  free(autocorr_mag);
  fclose(fp);
  fclose(fpa);
  fclose(fpm);

  return 0;
}
