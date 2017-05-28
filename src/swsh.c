// fix to long float problem win MinGW
#define __USE_MINGW_ANSI_STDIO 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <consts.h>
#include <test.h>
#include <funcs.h>

int main(int argc, char **argv)
{
	ldouble_t d;

	printf("\nUNIT TESTS:\n");
	test_atomic(1);
	test_factorial(1);
	test_binomial(1);
	test_Yslm(1);
	test_BC(1);

	return EXIT_SUCCESS;
}