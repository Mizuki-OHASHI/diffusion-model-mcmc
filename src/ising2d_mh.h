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

int ising2d_mh(double beta, double *energy_ptr, double *mag_ptr, double *autocorr_ptr, int verbose);
