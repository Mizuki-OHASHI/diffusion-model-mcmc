#include "xy2d_mh.h"

int main()
{
  double beta_min = 0.80, beta_max = 1.20, beta_step = 0.01;
  int n_beta = (int)((beta_max - beta_min) / beta_step) + 1;
  FILE *fp;

  fp = fopen("xy2d_mh_sweep.dat", "w");
  if (fp == NULL)
  {
    fprintf(stderr, "Failed to open output file.\n");
    return 1;
  }

  for (int i = 0; i < n_beta; i++)
  {
    double beta, energy, mag, autocorr;
    beta = beta_min + (double)i * beta_step;
    if (xy2d_mh(beta, &energy, &mag, &autocorr, 0) != 0)
    {
      fprintf(stderr, "Error in xy2d_mh for beta = %f\n", beta);
      fclose(fp);
      return 1;
    }
    fprintf(fp, "%f %f %f %f\n", beta, energy, mag, autocorr);
  }

  fclose(fp);
  return 0;
}
