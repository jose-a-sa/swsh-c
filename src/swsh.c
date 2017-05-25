// fix to long float problem win MinGW
#define __USE_MINGW_ANSI_STDIO 1
// still sizeof(long double)=12

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <gsl/gsl_tensor3_long_double.h>

int main(int argc, char **argv)
{
	long double x = 3.33L;
	gsl_tensor3_long_double *t = gsl_tensor3_long_double_alloc(10, 10, 10);

	gsl_tensor3_long_double_set(t,2,3,5,2.1);

	printf("%Lf\n", gsl_tensor3_long_double_get(t,2,3,5));
	printf("%Lf\n", gsl_tensor3_long_double_get(t,0,0,0));
	printf("%Lf\n", x);
	printf("%lu\n", sizeof(long double));
	printf("%lu\n", sizeof(double));

	//system("pause");

	return EXIT_SUCCESS;
}