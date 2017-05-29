// fix to long float problem win MinGW
#define __USE_MINGW_ANSI_STDIO 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <consts.h>
#include <test.h>
#include <diffeq.h>

int main(int argc, char **argv)
{
	real_t d;

	printf("\nUNIT TESTS:\n");
	test_atomic(1);
	test_factorial(1);
	test_binomial(1);
	test_Yslm(1);
	test_BC(1);

	vector_t y = VECTOR_ALLOC(3);
	VECTOR_SET(y,0,1.2);
	VECTOR_SET(y,1,2.0);
	VECTOR_SET(y,2,1.5);
	
	printf("\n\n%Lf\n\n", Dg(1, 2, y, 0.0, 1, 1, 1, 0.1) );

	return EXIT_SUCCESS;
}