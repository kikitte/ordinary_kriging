#include "mat_ops.h"

// Need to pass by reference
void free_ptr(int r, double **p)
{
    for (int i = 0; i < r; i++)
    {
        free(p[i]);
    }
    free(p);
}

// Create Matrix of Zeros
double **mat_zeros(int row, int col)
{
    double **mat = (double **)malloc(row * sizeof(double *));
    for (int i = 0; i < row; i++)
        mat[i] = (double *)malloc(col * sizeof(double));
    return mat;
}

void print_mat(int r, int c, double **mat)
{
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            printf("%f  ", mat[i][j]);
        }
        printf("\n");
    }
}