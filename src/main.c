#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "const.h"
#include "test.h"
#include "swsh.h"
#include "utils.h"

int main(int argc, char **argv)
{
	ldouble_t d;
	int i, j;

	printf("\nUNIT TESTS:\n");
	test_atomic(1);
	test_factorial(1);
	test_binomial(1);
	test_BC(1);
	test_Yslm(1);
	printf("\n");

	int s = -1, l = 5, m = 1;
	ldouble_t c_step = 0.1L, c_target = 10.0L+0.001L;
	show_const();

	vector_t X = vector_alloc(N_PTS);
	for (size_t i = 0; i < N_PTS; i++)
		vector_set(X, i, X_I + i * STEP);

	vector_t C = vector_alloc(N_PTS);

	char fname[50];
	sprintf(fname, "swsh_%d_%d_%d.csv", s, l, m);
	FILE *fp;
	fp = fopen(fname, "w");

	matrix_t Y = matrix_alloc(N_EQS, N_PTS);
	ldouble_t err_scale[N_EQS] = {1.0L, 1.0L, 1.0L};
	generate_initial_Y(Y, X, s, l, m, err_scale);

	ldouble_t c = 0.0L;
	while (c < c_target)
	{
		vector_set_all(C, c);
		desolve_relaxation(X, Y, s, l, m, c, err_scale, PRINT_ERROR);

		matrix_t m = matrix_alloc(2 + N_EQS, N_PTS);
		matrix_set_row(m, 0, C);
		matrix_set_row(m, 1, X);
		for(size_t i = 0; i < N_EQS; i++)
		{
			gsl_vector_long_double yi = matrix_row(Y, i).vector;
			matrix_set_row(m, 2+i, &yi);
		}
		output_matrix(fp, m);
		matrix_free(m);

		c += c_step;
	}

	fclose(fp);
	vector_free(X);
	vector_free(C);
	matrix_free(Y);

	printf("\nDone !\n\n");
	return EXIT_SUCCESS;
}