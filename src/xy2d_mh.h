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
int xy2d_mh(double beta, double *energy_ptr, double *helicity_ptr, double *autocorr_ptr, int verbose);
