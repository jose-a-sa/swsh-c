// fix to long long double problem with MinGW
#define __USE_MINGW_ANSI_STDIO 1

#include <math.h>

#include <consts.h>
#include <test.h>
#include <swsh.h>

#define MAX_ERROR 0.00001

bool test_unit(bool test, char *type, bool print)
{
    if (print)
    {
        printf("\t - checking %s : ", type);

        if (test)
            printf("PASS\n");
        else
            printf("FAIL\n");
    }

    return test;
}

void test_title(char *title, bool verbose)
{
    printf(" > %s", title);
    if (verbose)
        printf("\n");
}

void test_subtitle(char *subtitle, bool verbose)
{
    if (verbose)
        printf("\t * %s\n", subtitle);
}

void test_print(bool test, bool verbose)
{
    if (!verbose)
    {
        if (test)
            printf(" : PASSED\n");
        else
            printf(" : FAILED\n");
    }
}

bool test_atomic(bool verbose)
{
    bool test = true;

    test_title("Size of atomic types", verbose);

    test &= test_unit(sizeof(real_t) == 16, "sizeof(long long double)==16", verbose);
    test &= test_unit(sizeof(real_t) >= 12, "sizeof(long long double)>=12", verbose);

    test_print(test, verbose);
    return test;
}

bool test_factorial(bool verbose)
{
    bool test = true;

    test_title("Factorial", verbose);

    test &= test_unit(factorial(0) == 1, "factorial(0)", verbose);
    test &= test_unit(factorial(5) == 120, "factorial(5)", verbose);
    test &= test_unit(factorial(10) == 3628800, "factorial(10)", verbose);

    test_print(test, verbose);
    return test;
}

bool test_binomial(bool verbose)
{
    bool test = true;

    test_title("Binomial", verbose);

    test &= test_unit(binomial(0, 0) == 1, "binomial(0,0)", verbose);
    test &= test_unit(binomial(10, 3) == 120, "binomial(10,3)", verbose);
    test &= test_unit(binomial(-8, 5) == -792, "binomial(-8,5)", verbose);

    test_print(test, verbose);
    return test;
}

bool test_Yslm(bool verbose)
{
    bool test = true, check_buf;
    char buf[120];
    int i;

    const int s_test[] = { 1,  1, -4, -1};
    const int l_test[] = { 3,  6,  5,  3};
    const int m_test[] = { 2,  5,  2, -3};
    real_t x_test[] = {0.1, -0.2, 0.3, -0.4};
    real_t Y_value[] = {0.343445, 0.057114, -0.155762, -0.42492};
    real_t YPrime_value[] = {0.376269, -1.71342, 1.42068, -0.101171};

    test_title("Spin-Weighted Spherical Harmonics", verbose);

    test_subtitle("Yslm(s,l,m,x) function", verbose);
    for (i = 0; i < LEN(s_test); i++)
    {
        if (verbose) printf("\t");
        sprintf(buf, "Yslm(%2d,%2d,%2d,%4.1Lf)", s_test[i], l_test[i], m_test[i], x_test[i]);
        check_buf = test_unit(fabsl(Yslm(s_test[i], l_test[i], m_test[i], x_test[i]) / Y_value[i] - 1) < MAX_ERROR, buf, verbose);
        test &= check_buf;
    }

    test_subtitle("YslmPrime(s,l,m,x) function", verbose);
    for (i = 0; i < LEN(s_test); i++)
    {
        if (verbose) printf("\t");
        sprintf(buf, "YslmPrime(%2d,%2d,%2d,%4.1Lf)", s_test[i], l_test[i], m_test[i], x_test[i]);
        check_buf = test_unit(fabsl(YslmPrime(s_test[i], l_test[i], m_test[i], x_test[i]) / YPrime_value[i] - 1) < MAX_ERROR, buf, verbose);
        test &= check_buf;
    }

    test_print(test, verbose);

    return test;
}

bool test_BC(bool verbose)
{
    bool test = true, check_buf;
    char buf[120];
    int i;

    const int s_test[] = { 1,  1, -4, -1};
    const int l_test[] = { 3,  6,  5,  3};
    const int m_test[] = { 2,  5,  2, -3};
    real_t BClimZp_value[] = {1.18009, -1.80754, 0.607692, -0.361326};
    real_t BClimZm_value[] = {-0.590044, 1.2911, -1.41795, -0.361326};

    test_title("Boundary conditions", verbose);

    test_subtitle("Limit of (1-x)^kp (1+x)^km Yslm(x) at x=-1", verbose);
    for (i = 0; i < LEN(s_test); i++)
    {
        if (verbose) printf("\t");
        sprintf(buf, "BClimZm(%2d,%2d,%2d)", s_test[i], l_test[i], m_test[i]);
        check_buf = test_unit(fabsl(BClimZm(s_test[i], l_test[i], m_test[i]) / BClimZm_value[i] - 1) < MAX_ERROR, buf, verbose);
        test &= check_buf;
    }

    test_subtitle("Limit of (1-x)^kp (1+x)^km Yslm(x) at x=+1", verbose);
    for (i = 0; i < LEN(s_test); i++)
    {
        if (verbose) printf("\t");
        sprintf(buf, "BClimZp(%2d,%2d,%2d)", s_test[i], l_test[i], m_test[i]);
        check_buf &= test_unit(fabsl(BClimZp(s_test[i], l_test[i], m_test[i]) / BClimZp_value[i] - 1) < MAX_ERROR, buf, verbose);
        test &= check_buf;
    }

    test_print(test, verbose);

    return test;
}