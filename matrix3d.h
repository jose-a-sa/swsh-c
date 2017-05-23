#include <gsl/gsl_matrix_long_double.h>
#include <stddef.h>

long double gsl_3dmatrix_get(gsl_matrix_long_double **, int, int, int);
void gsl_3dmatrix_set(gsl_matrix_long_double **, int, int, int, long double);
gsl_matrix_long_double **gsl_3dmatrix_alloc(size_t, size_t, size_t);
gsl_matrix_long_double **gsl_3dmatrix_calloc(size_t, size_t, size_t);
void gsl_3dmatrix_free(gsl_matrix_long_double **, size_t);