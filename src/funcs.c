// fix to long long double problem with MinGW
#define __USE_MINGW_ANSI_STDIO 1

#include <math.h>
#include <stdint.h>

#include <gsl/gsl_errno.h>

#include <consts.h>
#include <funcs.h>

uint64_t factorial(int n)
{
    if (n == 0)
        return 1;
    else if (n < 0)
        GSL_ERROR("\nfactorial function input domain error.\n", GSL_EDOM);
    else
    {
        uint64_t i, p = 1;
        for (i = 1; i <= n; i++)
            p *= i;
        return p;
    }
}

int64_t binomial(int n, int r)
{
    if (r == 0 || r == n)
        return 1;

    if (n < 0)
    {
        if (r >= 0)
            return (2 * ((r + 1) % 2) - 1) * binomial(-n + r - 1, r);
        else if (r <= n)
            return (2 * ((n - r + 1) % 2) - 1) * binomial(-r - 1, n - r);
        else
            return 0;
    }
    else
    {
        if (r < 0 || r > n)
            return 0;
        else
            return factorial(n) / factorial(r) / factorial(n - r);
    }
}

ldouble_t Yslm(int s, int l, int m, ldouble_t x)
{
    if (l < MAX(abs(s), abs(m)))
        GSL_ERROR("\nYslm function s,l,m input error.\n", GSL_EDOM);

    ldouble_t result, sum = 0.0;
    int r;

    for (r = 0; r <= l - s; r++)
        sum += binomial(l - s, r) * binomial(l + s, r + s - m) * powl(-1.0, l - r - s + m) * powl((1.0 + x) / (1.0 - x), r + (s - m) / 2.0);

    result = sqrtl((factorial(l + m) * factorial(l - m) * (2 * l + 1)) / (factorial(l + s) * factorial(l - s) * 4.0 * PI)) * powl((1.0 - x) / 2.0, 1.0 * l) * sum;

    return result;
}

ldouble_t YslmPrime(int s, int l, int m, ldouble_t x)
{
    if (l < MAX(abs(s), abs(m)))
        GSL_ERROR("\nYslmPrime function s,l,m input error.\n", GSL_EDOM);

    ldouble_t result, sum = 0;
    int r;

    for (r = 0; r <= l - s; r++)
        sum += binomial(l - s, r) * binomial(l + s, r + s - m) * powl(-1.0, l - r - s + m + 1.0) * (l + m - 2.0 * r - s + l * x) * powl((1.0 + x) / (1.0 - x), r + (s - m) / 2.0) * powl((1.0 - x) / 2.0, l) / (1.0 - x * x);

    result = sqrtl((factorial(l + m) * factorial(l - m) * (2 * l + 1)) / (factorial(l + s) * factorial(l - s) * 4.0 * PI)) * sum;

    return result;
}

ldouble_t BClimYp(int s, int l, int m)
{
    if (l < MAX(abs(s), abs(m)))
        GSL_ERROR("\nBClimYm function s,l,m input error.\n", GSL_EDOM);

    ldouble_t result;

    if (m >= s && m + s >= 0)
        result = powl(-1.0, m) * powl(2.0, -m) * sqrtl(((1.0 + 2.0 * l) * factorial(m + l) * factorial(s + l)) / (4.0 * PI * factorial(l - m) * factorial(l - s))) / factorial(m + s);
    else if (m < s && m + s >= 0)
        result = powl(-1.0, m) * powl(2.0, -s) * sqrtl(((1.0 + 2.0 * l) * factorial(m + l) * factorial(s + l)) / (4.0 * PI * factorial(l - m) * factorial(l - s))) / factorial(m + s);
    else if (m >= s && m + s < 0)
        result = powl(-1.0, s) * powl(2.0, s) * sqrtl(((1.0 + 2.0 * l) * factorial(-m + l) * factorial(-s + l)) / (4.0 * PI * factorial(l + m) * factorial(l + s))) / factorial(-m - s);
    else
        result = powl(-1.0, s) * powl(2.0, m) * sqrtl(((1.0 + 2.0 * l) * factorial(-m + l) * factorial(-s + l)) / (4.0 * PI * factorial(l + m) * factorial(l + s))) / factorial(-m - s);

    return result;
}

ldouble_t BClimYm(int s, int l, int m)
{
    if (l < MAX(abs(s), abs(m)))
        GSL_ERROR("\nBClimYm function s,l,m input error.\n", GSL_EDOM);

    ldouble_t result;

    if (m >= s && m + s >= 0)
        result = powl(-1.0, l) * powl(2.0, -m) * sqrtl(((1.0 + 2.0 * l) * factorial(m + l) * factorial(-s + l)) / (4.0 * PI * factorial(l - m) * factorial(l + s))) / factorial(m - s);
    else if (m >= s && m + s < 0)
        result = powl(-1.0, l) * powl(2.0, +s) * sqrtl(((1.0 + 2.0 * l) * factorial(m + l) * factorial(-s + l)) / (4.0 * PI * factorial(l - m) * factorial(l + s))) / factorial(m - s);
    else if (m + s >= 0 && m < s)
        result = powl(-1.0, m + s + l) * powl(2.0, -s) * sqrtl(((1.0 + 2.0 * l) * factorial(-m + l) * factorial(s + l)) / (4.0 * PI * factorial(l + m) * factorial(l - s))) / factorial(-m + s);
    else
        result = powl(-1.0, m + s + l) * powl(2.0, +m) * sqrtl(((1.0 + 2.0 * l) * factorial(-m + l) * factorial(s + l)) / (4.0 * PI * factorial(l + m) * factorial(l - s))) / factorial(-m + s);

    return result;
}