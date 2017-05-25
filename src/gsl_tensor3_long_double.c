#include <gsl/gsl_errno.h>
#include <gsl/gsl_tensor3_long_double.h>
#include <string.h>

gsl_tensor3_long_double *gsl_tensor3_long_double_alloc(const size_t n1, const size_t n2, const size_t n3)
{
    gsl_block_long_double *block;
    gsl_tensor3_long_double *t;

    if (n1 == 0)
        GSL_ERROR_VAL("tensor dimension n1 must be positive integer", GSL_EINVAL, 0);
    else if (n2 == 0)
        GSL_ERROR_VAL("tensor dimension n2 must be positive integer", GSL_EINVAL, 0);
    else if (n3 == 0)
        GSL_ERROR_VAL("tensor dimension n3 must be positive integer", GSL_EINVAL, 0);

    t = (gsl_tensor3_long_double *)malloc(sizeof(gsl_tensor3_long_double));

    if (t == 0)
        GSL_ERROR_VAL("failed to allocate space for tensor struct", GSL_ENOMEM, 0);

    block = gsl_block_long_double_alloc(n1 * n2 * n3);

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

gsl_tensor3_long_double *gsl_tensor3_long_double_calloc(const size_t n1, const size_t n2, const size_t n3)
{
    size_t i;

    gsl_tensor3_long_double *t = gsl_tensor3_long_double_alloc(n1, n2, n3);

    if (t == 0)
        return 0;

    /* initialize matrix to zero */
    memset(t->data, 0, n1 * n2 * n3 * sizeof(long double));

    for (i = 0; i < n1 * n2 * n3; i++)
        t->data[i] = 0;

    return t;
}

void gsl_tensor3_long_double_set_zero(gsl_tensor3_long_double *t)
{
    size_t i, j, k;
    long double *const data = t->data;
    const size_t n1 = t->size1;
    const size_t n2 = t->size2;
    const size_t n3 = t->size3;
    const size_t tda = t->tda;
    const size_t tdb = t->tdb;

    for (i = 0; i < n1; i++)
        for (j = 0; j < n2; j++)
            for (k = 0; k < n3; k++)
                *(long double *)(data + (i * tda * tdb + j * tdb + k)) = 0.0L;
}

void gsl_tensor3_long_double_set_all(gsl_tensor3_long_double *t, long double x)
{
    size_t i, j, k;
    long double *const data = t->data;
    const size_t n1 = t->size1;
    const size_t n2 = t->size2;
    const size_t n3 = t->size3;
    const size_t tda = t->tda;
    const size_t tdb = t->tdb;

    for (i = 0; i < n1; i++)
        for (j = 0; j < n2; j++)
            for (k = 0; k < n3; k++)
                *(long double *)(data + (i * tda * tdb + j * tdb + k)) = x;
}

void gsl_tensor3_long_double_set(gsl_tensor3_long_double *t, const size_t i, const size_t j, const size_t k, long double x)
{
    long double *const data = t->data;
    const size_t n1 = t->size1;
    const size_t n2 = t->size2;
    const size_t n3 = t->size3;
    const size_t tda = t->tda;
    const size_t tdb = t->tdb;

    *(long double *)(data + (i * tda * tdb + j * tdb + k)) = x;
}

long double gsl_tensor3_long_double_get(gsl_tensor3_long_double *t, const size_t i, const size_t j, const size_t k)
{
    long double *const data = t->data;
    const size_t n1 = t->size1;
    const size_t n2 = t->size2;
    const size_t n3 = t->size3;
    const size_t tda = t->tda;
    const size_t tdb = t->tdb;

    return *(long double *)(data + (i * tda * tdb + j * tdb + k));
}