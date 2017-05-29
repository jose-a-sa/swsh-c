#ifndef CONSTS_H
#define CONSTS_H

#include <gsl/gsl_matrix_long_double.h>
#include <gsl/gsl_vector_long_double.h>

typedef long double real_t;
typedef gsl_matrix_long_double *matrix_t;
typedef gsl_vector_long_double *vector_t;

#define MATRIX_GET gsl_matrix_long_double_get
#define MATRIX_SET gsl_matrix_long_double_set
#define MATRIX_ALLOC gsl_matrix_long_double_alloc
#define MATRIX_FREE gsl_matrix_long_double_free

#define VECTOR_GET gsl_vector_long_double_get
#define VECTOR_SET gsl_vector_long_double_set
#define VECTOR_ALLOC gsl_vector_long_double_calloc
#define VECTOR_FREE gsl_vector_long_double_free

// General functions
#define LEN(x) (sizeof(x) / sizeof((x)[0]))
#define MAX(a, b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define MIN(a, b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

// Numerical contants
#define PI 3.1415926535897932385L

// Related to differential equation
#define N_EQS 3
#define X_I -1.0
#define X_F 1.0
#define BC_I 1
#define BC_F 2 // (N_EQS - BC_I)

// Related to the relaxation method
#define N_PTS 401
#define STEP ((X_F - X_I) / (N_PTS - 1))

#endif // CONSTS_H