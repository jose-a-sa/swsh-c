#include <consts.h>
#include <test.h>
#include <funcs.h>

bool check_unit(bool test, char * type)
{
    printf("\t- checking %s : ", type);
    
    if(test)
        printf("PASS\n");
    else
        printf("FAIL\n");

    return test;  
}

bool test_atomic()
{
    bool test = true;

    printf("> TEST size of atomic types \n");

    test &= check_unit(sizeof(ldouble_t)==16, "sizeof(long long double)==16");
    test &= check_unit(sizeof(ldouble_t)>=12, "sizeof(long long double)>=12");

    return test;
}

bool test_factorial()
{
    bool test = true;

    printf("> TEST factorial \n");

    test &= check_unit(factorial(0)==1, "factorial(0)");
    test &= check_unit(factorial(5)==120, "factorial(5)");
    test &= check_unit(factorial(10)==3628800, "factorial(10)");

    return test;
}

bool test_binomial()
{
    bool test = true;

    printf("> TEST binomial \n");

    test &= check_unit(binomial(0,0)==1, "binomial(0,0)");
    test &= check_unit(binomial(10,3)==120, "binomial(10,3)");

    return test;
}