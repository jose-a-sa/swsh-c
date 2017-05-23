#include <stdlib.h>
#include <stddef.h>
#include <gsl/gsl_matrix.h>
#include "matrix3d.h"

long double gsl_3dmatrix_get(gsl_matrix_long_double **m_vec ,int i, int j, int k)
{
	return gsl_matrix_long_double_get(m_vec[i],j,k);
}

void gsl_3dmatrix_set(gsl_matrix_long_double **m_vec, int i, int j, int k, long double x)
{
	gsl_matrix_long_double_set(m_vec[i],j,k,x);
}

gsl_matrix_long_double **gsl_3dmatrix_alloc(size_t n1, size_t n2, size_t n3)
{
	int p;
	gsl_matrix_long_double **m_vec = malloc((size_t) n1*sizeof(gsl_matrix_long_double)-1);
	for ( p=0; p<=n1; p++)
		m_vec[p] = gsl_matrix_long_double_alloc(n2 ,n3);
	return m_vec;
}

gsl_matrix_long_double **gsl_3dmatrix_calloc(size_t n1, size_t n2, size_t n3)
{
	int p;
	gsl_matrix_long_double **m_vec = malloc((size_t) n1*sizeof(gsl_matrix_long_double)-1);
	for ( p=0; p<=n1; p++)
		m_vec[p] = gsl_matrix_long_double_calloc(n2, n3);
	return m_vec;
}

void gsl_3dmatrix_free(gsl_matrix_long_double **m, size_t n1)
{
	int p;
	for (p=0; p<=n1; p++)
		gsl_matrix_long_double_free(m[p]);
	free(m);
}