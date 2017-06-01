// fix to long float problem win MinGW
#define __USE_MINGW_ANSI_STDIO 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <const.h>
#include <test.h>
#include <swsh.h>

int main(int argc, char **argv)
{
	real_t d;
	int i, j;

	show_const();

	printf("\nUNIT TESTS:\n");
	test_atomic(1);
	test_factorial(1);
	test_binomial(1);
	test_Yslm(1);
	test_BC(1);
	

	matrix_t Y = MATRIX_ALLOC(N_EQS, N_PTS);
	for (i = 0; i < N_EQS; i++)
	{
		for (j = 0; j < N_PTS; j++)
			MATRIX_SET(Y, i, j, 1.0L * (i + 1) * (j + 1));
	}

	for (j = 0; j < N_PTS; j++)
		printf("%Lf ", Sk(j, 1, 2, Y, -1, 1, 1, 0.1));

	return EXIT_SUCCESS;
}