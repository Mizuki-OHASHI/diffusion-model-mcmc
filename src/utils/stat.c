// utils/autocorrelation.c
// 各種統計関数を定義する

// 配列 data の自己相関関数を計算する関数
// data: 配列
// mean: 平均値
// n: 配列の長さ
// lag: ラグ (遅れ)
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
