#ifndef GSL_TENSOR3_H
#define GSL_TENSOR3_H

#include <gsl/gsl_block.h>

typedef struct
{
	size_t size1;
	size_t size2;
	size_t size3;
	size_t tda;
	size_t tdb;
	double *data;
	gsl_block *block;
	int owner;
} gsl_tensor3;

gsl_tensor3 *gsl_tensor3_alloc(const size_t, const size_t, const size_t);
gsl_tensor3 *gsl_tensor3_calloc(const size_t, const size_t, const size_t);
void gsl_tensor3_set_zero(gsl_tensor3 *);
void gsl_tensor3_set_all(gsl_tensor3 *, double);
void gsl_tensor3_set(gsl_tensor3 *, const size_t, const size_t, const size_t, double);
double gsl_tensor3_get(gsl_tensor3 *, const size_t, const size_t, const size_t);


#endif // GSL_TENSOR3_H
