#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include "matrix3d.h"

int main(void)
{
	int i, j, k;
	gsl_matrix_float * m = gsl_matrix_float_calloc(10, 10);
	gsl_matrix_float_set(m, 3, 3, 33.3);

	// for (i = 0; i < 10; i++)
	// 	for (j = 0; j < 10; j++)
	// 		for (k = 0; k < 10; k++)
	// 			matrix3d_set(m, i, j, k, i * j * k);

	printf("m(%d,%d) = %f\n", 3, 3, gsl_matrix_float_get(m, 3, 3) );

	//system("pause");

	return 0;
}