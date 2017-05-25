#include <funcs.h>
#include <math.h>

#include <gsl/gsl_errno.h>


uint64_t factorial(int n)
{
    if (n == 0)
        return 1;
    else if (n < 0)
        GSL_ERROR("\nfactorial function input domain error.\n", GSL_EDOM);
    else
    {
        int i;
        uint64_t p = 1;
        for (i = 1; i <= n; i++)
            p *= i;
        return p;
    }
}

uint64_t binomial(int n, int r)
{
    if (n < 0 || r < 0 || r > n)
        return 0;
    else if (r == 0 || r == n)
        return 1;
    else
        return factorial(n) / (factorial(n - r) * factorial(r));
}

ldouble_t Yslm(int s, int l, int m, ldouble_t x)
{
    ldouble_t result, sum = 0;
    int r;

    for (r = 0; r <= l - s; r++)
        sum += binomial(l - s, r) * binomial(l + s, r + s - m) * powl(-1.0, l - r - s + m) * powl((1.0 + x) / (1.0 - x), r + (s - m) / 2.0);

    result = sqrtl((factorial(l + m) * factorial(l - m) * (2 * l + 1)) / (factorial(l + s) * factorial(l - s) * 4.0 * PI)) * powl((1.0 - x) / 2.0, 1.0 * l) * sum;

    return result;
}

ldouble_t YslmPrime(int32_t s, int32_t l, int32_t m, ldouble_t x)
{
    ldouble_t result, sum = 0;
    int r;

    for (r = 0; r <= l - s; r++)
        sum += binomial(l - s, r) * binomial(l + s, r + s - m) * powl(-1.0, l - r - s + m + 1.0) * (l + m - 2.0 * r - s + l * x) * powl((1.0 + x) / (1.0 - x), r + (s - m) / 2.0) * powl((1.0 - x) / 2.0, l) / (1.0 - x * x);

    result = sqrtl((factorial(l + m) * factorial(l - m) * (2 * l + 1)) / (factorial(l + s) * factorial(l - s) * 4.0 * PI)) * sum;

    return result;
}