//fix to long float problem win MinGW
//still sizeof(long double)=12
#define __USE_MINGW_ANSI_STDIO 1 

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <gsl/gsl_tensor3.h>
#include <gsl/gsl_matrix.h>  

int main(int argc, char **argv)
{
	long double x = 3.33L;
	gsl_tensor3 *t = gsl_tensor3_alloc(10, 10, 10);

	gsl_tensor3_set(t,2,3,5,2.1);

	printf("%f\n", gsl_tensor3_get(t,2,3,5));
	printf("%f\n", gsl_tensor3_get(t,0,0,0));
	printf("%Lf\n", x);
	printf("%lu\n", sizeof(long double));
	printf("%lu\n", sizeof(double));

	return EXIT_SUCCESS;
}