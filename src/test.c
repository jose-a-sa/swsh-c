#include <math.h>

#include "const.h"
#include "test.h"
#include "swsh.h"

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
        printf("     * %s\n", subtitle);
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
    printf("\t - checking sizeof(long double) : %lu\n", sizeof(ldouble_t));
    test &= test_unit(sizeof(ldouble_t) >= 16, "sizeof(long double)>=16", verbose);
    test &= test_unit(sizeof(ldouble_t) >= 12, "sizeof(long double)>=12", verbose);
    test &= test_unit(sizeof(ldouble_t) >= 8, "sizeof(long double)>=8", verbose);

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
    test &= test_unit(binomial(-5, 9) == -715, "binomial(-5,9)", verbose);

    test_print(test, verbose);
    return test;
}

bool test_BC(bool verbose)
{
    bool test = true, check_buf;
    char buf[120];
    int i;

    const int s_test[] = {1, 1, -4, -1};
    const int l_test[] = {3, 6, 5, 3};
    const int m_test[] = {2, 5, 2, -3};
    ldouble_t limYp_value[] = {1.18009, -1.80754, 0.607692, -0.361326};
    ldouble_t limYm_value[] = {-0.590044, 1.2911, -1.41795, -0.361326};
    ldouble_t limYpPrime_value[] = {0.885065, -1.54932, 1.01282, 0.};
    ldouble_t limYmPrime_value[] = {0.885065, -1.54932, 1.01282, 0.};

    test_title("Boundary conditions", verbose);

    test_subtitle("Limit of (1-x)^kp (1+x)^km Yslm(x) at x=-1", verbose);
    for (i = 0; i < LEN(s_test); i++)
    {
        sprintf(buf, "lim_Yslm(%2d,%2d,%2d, -1)", s_test[i], l_test[i], m_test[i]);
        check_buf = test_unit(possibly_same(limYm_value[i], lim_Yslm(s_test[i], l_test[i], m_test[i], -1)), buf, verbose);
        test &= check_buf;
    }

    test_subtitle("Limit of (1-x)^kp (1+x)^km Yslm(x) at x=+1", verbose);
    for (i = 0; i < LEN(s_test); i++)
    {
        sprintf(buf, "lim_Yslm(%2d,%2d,%2d, +1)", s_test[i], l_test[i], m_test[i]);
        check_buf &= test_unit(possibly_same(limYp_value[i], lim_Yslm(s_test[i], l_test[i], m_test[i], +1)), buf, verbose);
        test &= check_buf;
    }

    test_subtitle("Limit of d/dx((1-x)^kp (1+x)^km Yslm(x)) at x=-1", verbose);
    for (i = 0; i < LEN(s_test); i++)
    {
        sprintf(buf, "lim_YslmPrime(%2d,%2d,%2d, -1)", s_test[i], l_test[i], m_test[i]);
        check_buf = test_unit(possibly_same(limYmPrime_value[i], lim_YslmPrime(s_test[i], l_test[i], m_test[i], -1)), buf, verbose);
        test &= check_buf;
    }

    test_subtitle("Limit of d/dx((1-x)^kp (1+x)^km Yslm(x)) at x=+1", verbose);
    for (i = 0; i < LEN(s_test); i++)
    {
        sprintf(buf, "lim_YslmPrime(%2d,%2d,%2d, +1)", s_test[i], l_test[i], m_test[i]);
        check_buf &= test_unit(possibly_same(limYpPrime_value[i], lim_YslmPrime(s_test[i], l_test[i], m_test[i], +1)), buf, verbose);
        test &= check_buf;
    }

    test_print(test, verbose);

    return test;
}

bool test_Yslm(bool verbose)
{
    bool test = true, check_buf;
    char buf[120];
    int i;

    const int s_test[] = {1, 1, -3, -1, -2};
    const int l_test[] = {3, 5, 5, 3, 6};
    const int m_test[] = {2, 1, 2, -3, 1};
    ldouble_t x_test[] = {-0.1392, 0.8345, 0.6122, -0.3613, 0.54256};
    ldouble_t Yslm_value[] = {0.171821, -2.81394, 0.103783, -0.361326, -0.18028};
    ldouble_t YslmPrime_value[] = {0.885065, -18.798, -0.91183, 0., -1.61855};

    test_title("Function values", verbose);

    test_subtitle("Value of Yslm(x)", verbose);
    for (i = 0; i < LEN(s_test); i++)
    {
        sprintf(buf, "Yslm(%2d,%2d,%2d, %+Lf)", s_test[i], l_test[i], m_test[i], x_test[i]);
        check_buf = test_unit(possibly_same(Yslm_value[i], Yslm(s_test[i], l_test[i], m_test[i], x_test[i])), buf, verbose);
        test &= check_buf;
    }

    test_subtitle("Value of YslmPrime(x)", verbose);
    for (i = 0; i < LEN(s_test); i++)
    {
        sprintf(buf, "YslmPrime(%2d,%2d,%2d, %+Lf)", s_test[i], l_test[i], m_test[i], x_test[i]);
        check_buf = test_unit(possibly_same(YslmPrime_value[i], YslmPrime(s_test[i], l_test[i], m_test[i], x_test[i])), buf, verbose);
        test &= check_buf;
    }

    test_print(test, verbose);

    return test;
}

void show_const()
{
    printf("DEFINED CONSTANTS:\n");

    printf(" > Size of the mesh: %d\n", N_PTS);
    printf(" > Mech configuration Xi, Xf, step: %Lf, %Lf, %Lf\n", 1.0L * X_I, 1.0L * X_F, STEP);

    printf(" > Number of equations: %d\n", N_EQS);
    printf(" > Boundary conditions at Xi, Xf: %d, %d\n", BC_I, BC_F);

    printf("\n");
}