#include <gsl/gsl_block_long_double.h>
#include <gsl/gsl_errno.h>

typedef struct
{
	size_t size1;
	size_t size2;
	size_t size3;
	size_t tda;
	size_t tdb;
	long double *data;
	gsl_block_long_double *block;
	int owner;
} gsl_tensor3_long_double;

gsl_tensor3_long_double *gsl_tensor3_long_double_alloc(const size_t, const size_t, const size_t);
gsl_tensor3_long_double *gsl_tensor3_long_double_calloc(const size_t, const size_t, const size_t);
void gsl_tensor3_long_double_set_zero(gsl_tensor3_long_double *);
void gsl_tensor3_long_double_set_all(gsl_tensor3_long_double *, long double);
void gsl_tensor3_long_double_set(gsl_tensor3_long_double *, const size_t, const size_t, const size_t, long double);
long double gsl_tensor3_long_double_get(gsl_tensor3_long_double *, const size_t, const size_t, const size_t);
