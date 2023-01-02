#ifndef CONSTS_H
#define CONSTS_H

#include <gsl/gsl_matrix_long_double.h>
#include <gsl/gsl_vector_long_double.h>

typedef long double ldouble_t;
typedef gsl_matrix_long_double* matrix_t;
typedef gsl_vector_long_double* vector_t;

#define matrix_get gsl_matrix_long_double_get
#define matrix_set gsl_matrix_long_double_set
#define matrix_alloc gsl_matrix_long_double_calloc
#define matrix_free gsl_matrix_long_double_free
#define matrix_set_row gsl_matrix_long_double_set_row
#define matrix_set_col gsl_matrix_long_double_set_col
#define matrix_set_all gsl_matrix_long_double_set_all
#define matrix_swap_rows gsl_matrix_long_double_swap_rows
#define matrix_swap_columns gsl_matrix_long_double_swap_columns
#define matrix_column gsl_matrix_long_double_column
#define matrix_row gsl_matrix_long_double_row
#define matrix_scale gsl_matrix_long_double_scale
#define matrix_add gsl_matrix_long_double_add

#define vector_get gsl_vector_long_double_get
#define vector_set gsl_vector_long_double_set
#define vector_alloc gsl_vector_long_double_calloc
#define vector_free gsl_vector_long_double_free
#define vector_set_all gsl_vector_long_double_set_all

// General functions
#define LEN(x) (sizeof(x) / sizeof((x)[0]))
#define MAX(a, b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define MIN(a, b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

// General contants
#define ZERO_EPS 1e-10L
#define SAME_REL_EPS 1e-5L

// Numerical contants
#define PI 3.1415926535897932385L

// Related to differential equation
#define N_EQS 3
#define X_I -1.0L
#define X_F 1.0L
#define BC_I 1
#define BC_F 2

// Related to the relaxation method
#define N_PTS 401
#define STEP (1.0L * ((X_F) - (X_I)) / (N_PTS - 1))
#define CONV_EPS 1e-8L
#define CONTROL_ERR 0.1L
#define ITER_MAX 1000

//LOG FLAGS
#define NO_LOG 0
#define LOG_SMATRIX (1<<0)
#define LOG_SMATRIX_GAUSS (1<<1)
#define LOG_DELTA_Y (1<<2)
#define PRINT_ERROR (1<<3)

#endif // CONSTS_H