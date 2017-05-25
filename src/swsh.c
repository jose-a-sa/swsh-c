// fix to long float problem win MinGW
#define __USE_MINGW_ANSI_STDIO 1
// still sizeof(long double)=12

#include <stdio.h>
#include <stdlib.h>
#include <consts.h>

#include <funcs.h>

int main(int argc, char **argv)
{
	FILE *f1;
	ldouble_t d;
	
	f1 = fopen("data/Y011.csv", "w+");	

	printf("%.12Lf\n", PI);
	printf("%lu\n", sizeof(ldouble_t));
	printf("%lu\n", sizeof(double));

	for(d = -0.999999999999; d<0.99999999999; d+=0.001)
		fprintf(f1,"%Lf, %Lf, %Lf\n", d, Yslm(1, 7, 2, d), YslmPrime(1, 7, 2, d));

	fclose(f1);

	return EXIT_SUCCESS;
}