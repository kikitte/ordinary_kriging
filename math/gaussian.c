#include <stdio.h>
#include <math.h>

/**
 * 使用高斯消元法求解线性方程组Ax=B
 * 将求解结果写入x数组
 * @param N 方阵阶数
 * @param mat N阶方阵
 * @param Tol small tolerance number to detect failure when the matrix is near degenerate
 * @return 真值：方程有唯一解 假值：方程无解或解不为一
 */
int gaussian_elimination_solve(int N, double **A, double *b, double *x, double Tol)
{
  int i, j, k, imax;
  double maxA, absA, f, *ptr;

  for (i = 0; i < N; ++i)
  {
    maxA = 0.0;
    imax = i;

    // 寻找当前列当下及以下所有行中最大值者, 记录绝对值以及行号(k)
    for (k = i; k < N; k++)
      if ((absA = fabs(A[k][i])) > maxA)
      {
        maxA = absA;
        imax = k;
      }

    if (maxA < Tol)
      return 0;

    if (imax != i)
    {
      // 交换行：i , imax
      ptr = A[i];
      A[i] = A[imax];
      A[imax] = ptr;
    }

    /* Do for all rows below pivot: */
    for (j = i + 1; j < N; ++j)
    {
      f = A[j][i] / A[i][i];
      A[j][i] = 0;

      for (k = i + 1; k < N; ++k)
        A[j][k] -= f * A[i][k];
    }
  }

  // back substitution
  for (j = N - 1; j >= 0; --j)
  {
    x[j] = b[j];
    for (i = j + 1; i < N; i++)
      x[j] -= A[j][i] * x[i];

    x[i] /= A[i][i];
  }

  return 1;
}