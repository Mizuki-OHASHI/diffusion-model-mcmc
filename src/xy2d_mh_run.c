#include "xy2d_mh.h"

int main()
{
  double beta = 0.8;
  double energy, mag, autocorr;
  return xy2d_mh(beta, &energy, &mag, &autocorr, 1);
}
