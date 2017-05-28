// fix to long long double problem with MinGW
#define __USE_MINGW_ANSI_STDIO 1

#include <gsl/gsl_tensor3.h>
#include <gsl/gsl_errno.h>
#include <string.h>

gsl_tensor3 *gsl_tensor3_alloc(const size_t n1, const size_t n2, const size_t n3)
{
	gsl_block *block;
	gsl_tensor3 *t;

	if (n1 == 0)
		GSL_ERROR_VAL("tensor dimension n1 must be positive integer", GSL_EINVAL, 0);
	else if (n2 == 0)
		GSL_ERROR_VAL("tensor dimension n2 must be positive integer", GSL_EINVAL, 0);
	else if (n3 == 0)
		GSL_ERROR_VAL("tensor dimension n3 must be positive integer", GSL_EINVAL, 0);

	t = (gsl_tensor3 *)malloc(sizeof(gsl_tensor3));

	if (t == 0)
		GSL_ERROR_VAL("failed to allocate space for tensor struct", GSL_ENOMEM, 0);

	block = gsl_block_alloc(n1 * n2 * n3);

	if (block == 0)
		GSL_ERROR_VAL("failed to allocate space for block", GSL_ENOMEM, 0);

	t->data = block->data;
	t->size1 = n1;
	t->size2 = n2;
	t->size3 = n3;
	t->tda = n2;
	t->tdb = n3;
	t->block = block;
	t->owner = 1;

	return t;
}

gsl_tensor3 *gsl_tensor3_calloc(const size_t n1, const size_t n2, const size_t n3)
{
	size_t i;

	gsl_tensor3 * t = gsl_tensor3_alloc(n1, n2, n3);

	if (t == 0)
		return 0;

	/* initialize matrix to zero */
	memset(t->data, 0,  n1*n2*n3*sizeof(double));

	for(i = 0; i < n1*n2*n3; i++)
		t->data[i] = 0;

	return t;
}

void gsl_tensor3_set_zero(gsl_tensor3 * t)
{
    size_t i, j, k;
    double * const data = t->data;
    const size_t n1 = t->size1;
    const size_t n2 = t->size2;
    const size_t n3 = t->size3;
    const size_t tda = t->tda;
    const size_t tdb = t->tdb;

    for (i = 0; i < n1; i++)
        for (j = 0; j < n2; j++)
            for (k = 0; k < n3; k++)
                *(double *) (data + (i*tda*tdb + j*tdb + k)) = 0.0;
}

void gsl_tensor3_set_all(gsl_tensor3 * t, double x)
{
    size_t i, j, k;
    double * const data = t->data;
    const size_t n1 = t->size1;
    const size_t n2 = t->size2;
    const size_t n3 = t->size3;
    const size_t tda = t->tda;
    const size_t tdb = t->tdb;
    
    for (i = 0; i < n1; i++)
        for (j = 0; j < n2; j++)
            for (k = 0; k < n3; k++)
                *(double *) (data + (i*tda*tdb + j*tdb + k)) = x;
}

void gsl_tensor3_set(gsl_tensor3 * t, const size_t i, const size_t j, const size_t k, double x)
{
    double * const data = t->data;
    const size_t n1 = t->size1;
    const size_t n2 = t->size2;
    const size_t n3 = t->size3;
    const size_t tda = t->tda;
    const size_t tdb = t->tdb;
    
    *(double *) (data + (i*tda*tdb + j*tdb + k)) = x;
}

double gsl_tensor3_get(gsl_tensor3 * t, const size_t i, const size_t j, const size_t k)
{
    double * const data = t->data;
    const size_t n1 = t->size1;
    const size_t n2 = t->size2;
    const size_t n3 = t->size3;
    const size_t tda = t->tda;
    const size_t tdb = t->tdb;
    
    return *(double *) (data + (i*tda*tdb + j*tdb + k));
}