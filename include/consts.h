#ifndef CONSTS_h
#define CONSTS_h


#include <gsl/gsl_matrix_long_double.h>
#include <gsl/gsl_vector_long_double.h>

typedef long double ldouble_t;
typedef gsl_matrix_long_double* matrix_t;
typedef gsl_vector_long_double* vector_t;

#define MATRIX_GET gsl_matrix_long_double_get
#define MATRIX_SET gsl_matrix_long_double_set

#define VECTOR_GET gsl_vector_long_double_get
#define VECTOR_SET gsl_vector_long_double_set


#define PI 3.1415926535897932385L

#define N_PTS 401

#define kp fabsl(m+s)/2.0
#define km fabsl(m-s)/2.0


#endif // CONSTS_h