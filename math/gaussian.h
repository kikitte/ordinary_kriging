/**
 * 使用高斯消元法求解线性方程组Ax=B
 * A为增广矩阵（包含系数矩阵与等号右边的结果）
 * 将求解结果写入b数组
 * @param N 方阵阶数
 * @param A N X (N+1)阶方阵
 * @param x 结果输出
 * @param Tol small tolerance number to detect failure when the matrix is near degenerate
 * @return 真值：方程有唯一解 假值：方程无解或解不为一
 */
int gaussian_elimination_solve(int N, double **A, double *x, double Tol);