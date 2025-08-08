#include <stddef.h>

// utils/autocorrelation.c
// 各種統計関数を定義する

/** @brief
 * 配列 data の自己相関関数を計算します。
 * @param data 配列
 * @param mean 平均値
 * @param n    配列の長さ
 * @param lag  ラグ (遅れ)
 * @return double 計算された自己相関関数の値
 */
double autocorrelation(const double *data, double mean, int n, int lag)
{
  if (lag >= n || lag < 0)
    return 0.0;
  double num = 0.0, denom = 0.0;
  for (int i = 0; i < n - lag; ++i)
  {
    num += (data[i] - mean) * (data[i + lag] - mean);
  }
  for (int i = 0; i < n; ++i)
  {
    denom += (data[i] - mean) * (data[i] - mean);
  }
  if (denom == 0.0)
    return 0.0;
  return num * (double)n / ((double)(n - lag) * denom);
}

/** @brief
 * 自己相関関数 ρ(t) の配列から、積分自己相関時間 τ を計算します。
 * 計算式: τ = 1 + 2 * Σ_{t=1}^{W} ρ(t)
 * 自己整合法 (Self-Consistent Windowing) を用います。
 * これは、ラグ t がそれまでに計算された自己相関時間 τ の M 倍を
 * 超えた時点で和を打ち切る手法です。
 * @param corr 自己相関関数の配列 (t=0, 1, ..., n-1 の値を持つ)
 * @param n    配列 corr の長さ
 * @return double 計算された自己相関時間 τ
 */
double autocorrelation_time(const double *corr, int n)
{
  if (corr == NULL || n <= 1)
  {
    return -1.0;
  }

  double tau = corr[0];
  const double M = 6.0;

  for (int t = 1; t < n; ++t)
  {
    if ((double)t >= M * tau)
      break;
    tau += 2.0 * corr[t];
  }

  return tau;
}
