#include "ising2d_mh.h"

int main()
{
  double beta = 0.44;
  double energy, mag, autocorr;
  return ising2d_mh(beta, &energy, &mag, &autocorr, 1);
}
