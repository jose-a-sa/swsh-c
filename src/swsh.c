// fix to long long double problem with MinGW
#define __USE_MINGW_ANSI_STDIO 1

#include <math.h>
#include <stdint.h>

#include <gsl/gsl_errno.h>

#include <consts.h>
#include <swsh.h>

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

real_t Yslm(int s, int l, int m, real_t x)
{
    if (l < MAX(abs(s), abs(m)))
        GSL_ERROR("\nYslm function s,l,m input error.\n", GSL_EDOM);

    real_t sum = 0.0;
    int r;

    for (r = 0; r <= l - s; r++)
        sum += binomial(l - s, r) * binomial(l + s, r + s - m) * powl(-1.0, l - r - s + m) * powl((1.0 + x) / (1.0 - x), r + (s - m) / 2.0);

    return sqrtl((factorial(l + m) * factorial(l - m) * (2 * l + 1)) / (factorial(l + s) * factorial(l - s) * 4.0 * PI)) * powl((1.0 - x) / 2.0, 1.0 * l) * sum;
}

real_t YslmPrime(int s, int l, int m, real_t x)
{
    if (l < MAX(abs(s), abs(m)))
        GSL_ERROR("\nYslmPrime function s,l,m input error.\n", GSL_EDOM);

    real_t sum = 0.0;
    int r;

    for (r = 0; r <= l - s; r++)
        sum += binomial(l - s, r) * binomial(l + s, r + s - m) * powl(-1.0, l - r - s + m + 1.0) * (l + m - 2.0 * r - s + l * x) * powl((1.0 + x) / (1.0 - x), r + (s - m) / 2.0) * powl((1.0 - x) / 2.0, l) / (1.0 - x * x);

    return sqrtl((factorial(l + m) * factorial(l - m) * (2 * l + 1)) / (factorial(l + s) * factorial(l - s) * 4.0 * PI)) * sum;
}

// Limit of (1-x)^kp (1+x)^km Yslm(x) at x=+1
real_t BClimZp(int s, int l, int m)
{
    if (l < MAX(abs(s), abs(m)))
        GSL_ERROR("\nBClimZm function s,l,m input error.\n", GSL_EDOM);

    if (m >= s && m + s >= 0)
        return powl(-1.0, m) * powl(2.0, -m) * sqrtl(((1.0 + 2.0 * l) * factorial(m + l) * factorial(s + l)) / (4.0 * PI * factorial(l - m) * factorial(l - s))) / factorial(m + s);
    else if (m < s && m + s >= 0)
        return powl(-1.0, m) * powl(2.0, -s) * sqrtl(((1.0 + 2.0 * l) * factorial(m + l) * factorial(s + l)) / (4.0 * PI * factorial(l - m) * factorial(l - s))) / factorial(m + s);
    else if (m >= s && m + s < 0)
        return powl(-1.0, s) * powl(2.0, s) * sqrtl(((1.0 + 2.0 * l) * factorial(-m + l) * factorial(-s + l)) / (4.0 * PI * factorial(l + m) * factorial(l + s))) / factorial(-m - s);
    else
        return powl(-1.0, s) * powl(2.0, m) * sqrtl(((1.0 + 2.0 * l) * factorial(-m + l) * factorial(-s + l)) / (4.0 * PI * factorial(l + m) * factorial(l + s))) / factorial(-m - s);
}

// Limit of (1-x)^kp (1+x)^km Yslm(x) at x=-1
real_t BClimZm(int s, int l, int m)
{
    if (l < MAX(abs(s), abs(m)))
        GSL_ERROR("\nBClimZm function s,l,m input error.\n", GSL_EDOM);

    if (m >= s && m + s >= 0)
        return powl(-1.0, l) * powl(2.0, -m) * sqrtl(((1.0 + 2.0 * l) * factorial(m + l) * factorial(-s + l)) / (4.0 * PI * factorial(l - m) * factorial(l + s))) / factorial(m - s);
    else if (m >= s && m + s < 0)
        return powl(-1.0, l) * powl(2.0, +s) * sqrtl(((1.0 + 2.0 * l) * factorial(m + l) * factorial(-s + l)) / (4.0 * PI * factorial(l - m) * factorial(l + s))) / factorial(m - s);
    else if (m + s >= 0 && m < s)
        return powl(-1.0, m + s + l) * powl(2.0, -s) * sqrtl(((1.0 + 2.0 * l) * factorial(-m + l) * factorial(s + l)) / (4.0 * PI * factorial(l + m) * factorial(l - s))) / factorial(-m + s);
    else
        return powl(-1.0, m + s + l) * powl(2.0, +m) * sqrtl(((1.0 + 2.0 * l) * factorial(-m + l) * factorial(s + l)) / (4.0 * PI * factorial(l + m) * factorial(l - s))) / factorial(-m + s);
}

// Ratio of Zslm'(x)/Zslm(x) at x=-1
real_t BCratioZm(int s, int l, int m, real_t A, real_t c)
{
    if (l < MAX(abs(s), abs(m)))
        GSL_ERROR("\nBCratioZm function s,l,m input error.\n", GSL_EDOM);

    real_t kp, km;
    kp = fabsl(m + s) / 2.0;
    km = fabsl(m - s) / 2.0;

    return -(A + powl(c,2.0) + 2.0*c*s - (km + kp - s)*(1.0 + km + kp + s))/(2.0 + 4.0*km);
}

// Ratio of Zslm'(x)/Zslm(x) at x=+1
real_t BCratioZp(int s, int l, int m, real_t A, real_t c)
{
    if (l < MAX(abs(s), abs(m)))
        GSL_ERROR("\nBCratioZp function s,l,m input error.\n", GSL_EDOM);

    real_t kp, km;
    kp = fabsl(m + s) / 2.0;
    km = fabsl(m - s) / 2.0;

    return -(A + powl(c,2.0) - 2.0*c*s - (km + kp - s)*(1.0 + km + kp + s))/(2.0 + 4.0*kp);
}

// g contains all information about dif. equation: yi'=gi(x,y0,y1,..)
real_t g(int j, vector_t y, real_t x, int s, int l, int m, real_t c)
{
    if (j < 0 || j >= N_EQS || l < MAX(abs(s), abs(m)))
        GSL_ERROR("\ng function input error.\n", GSL_EDOM);

    real_t result, kp, km;
    kp = fabsl(m + s) / 2.0;
    km = fabsl(m - s) / 2.0;

    if (j == 0)
        result = VECTOR_GET(y, 1);
    else if (j == 1)
    {
        result = 2.0 * (kp - km + (1.0 + km + kp) * x) * VECTOR_GET(y, 1) + VECTOR_GET(y, 0) * ((km + kp) * (1.0 + km + kp) - powl(c * x, 2.0) + 2 * c * s * x - VECTOR_GET(y, 2) - s * (s + 1.0));
        result /= (1.0 - powl(x, 2.0));
    }
    else if (j == 2)
        result = 0.0;

    return result;
}

// jacobian of g
real_t Dg(int j, int i, vector_t y, real_t x, int s, int l, int m, real_t c)
{
    if (j < 0 || j >= N_EQS || i < 0 || i >= N_EQS || l < MAX(abs(s), abs(m)))
        GSL_ERROR("\ng function input error.\n", GSL_EDOM);

    real_t result, kp, km;
    kp = fabsl(m + s) / 2.0;
    km = fabsl(m - s) / 2.0;

    if (j == 0)
    {
        if (i == 1)
            result = 1.0;
        else
            result = 0.0;
    }
    else if (j == 1)
    {
        if (i == 0)
            result = (km + kp) * (1.0 + km + kp) - powl(c * x, 2.0) + 2 * c * s * x - VECTOR_GET(y, 2) - s * (s + 1.0);
        else if (i == 1)
            result = 2.0 * (kp - km + (1.0 + km + kp) * x);
        else if (i == 2)
            result = -VECTOR_GET(y, 0);

        result /= (1.0 - powl(x, 2.0));
    }
    else if (j == 2)
        result = 0.0;

    return result;
}