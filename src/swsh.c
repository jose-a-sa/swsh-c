// fix to long float problem win MinGW
#define __USE_MINGW_ANSI_STDIO 1
// still sizeof(long double)=12

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <consts.h>
#include <test.h>
#include <funcs.h>

int main(int argc, char **argv)
{
	ldouble_t d;

	test_atomic();
	test_factorial();
	test_binomial();

	return EXIT_SUCCESS;
}