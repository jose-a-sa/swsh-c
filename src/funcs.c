#include <funcs.h>
#include <math.h>

#include <gsl/gsl_errno.h>


int64_t powi(int32_t base, int32_t exp)
{
    int64_t result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

uint64_t factorial(int32_t n)
{
    if(n == 0)
        return 1;
    else if(n < 0)
        GSL_ERROR("\nfactorial function input domain error.\n", GSL_EDOM);
    else 
    {
        uint32_t i;
        uint64_t p = 1;
        for(i = 1; i <= n; i++)
            p *= i;
        return p;
    }
}


uint64_t binomial(int32_t n, int32_t r)
{
    if(n < 0 || r < 0 || r > n)
        return 0;
    else if (r == 0 || r == n)
        return 1;
    else
        return factorial(n)/(factorial(n-r)*factorial(r));
}


ldouble_t Yslm(int32_t s, int32_t l, int32_t m, ldouble_t x)
{
    ldouble_t result, sum = 0.0L;
    uint32_t r;
    for(r = 0; r <= l-s ; r++)
        sum += binomial(l-s,r) * binomial(l+s, r+s-m) * powi(-1,l-r-s) * powl( (1+x)/(1-x), r + (s-m)/2)
    result = powi(-1,m) * sqrtl((fac(l+m)*fac(l-m)*(2*l+1)/(fac(l+s)*fac(l-s)*4.0*M_PI))*powl( (1-x)/2, 1.0L*l) * sum;
    return ans;
}