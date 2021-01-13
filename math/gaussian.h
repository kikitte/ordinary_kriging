/**
 * 使用高斯消元法求解线性方程组Ax=B
 * 将求解结果写入x数组
 * @param N 方阵阶数
 * @param mat N阶方阵
 * @param Tol small tolerance number to detect failure when the matrix is near degenerate
 * @return 真值：方程有唯一解 假值：方程无解或解不为一
 */
int gaussian_elimination_solve(int N, double **A, double *b, double *x, double Tol);