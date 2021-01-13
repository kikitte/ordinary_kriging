#ifndef MAT_OPS
#define MAT_OPS
#include <inttypes.h>
#include <stdbool.h>
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "signal.h"

void free_ptr(int r, double **p);

double **mat_zeros(int row, int col);

void print_mat(int r, int c, double **mat);

#endif
