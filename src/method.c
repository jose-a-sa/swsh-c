#include <consts.h>
#include <method.h>
#include <diffeq.h>

#define N_PTS 401
#define STEP ((XF - XI) / (N_PTS - 1))

real_t Ek(int k, int j, matrix_t Y, int s, int l, int m, real_t c)
{
    int i;

    vector_t yk = VECTOR_ALLOC( N_EQS );
    for( i = 0; i< N_EQS; i++)
        VECTOR_SET(yk, i, (MATRIX_GET(Y, 0, k - 1) + MATRIX_GET(Y, 0, k - 1)) / 2.0 )
    
    real_t xk = (k - 0.5) * STEP + X_I;


    return
}
real_t Sk(int, int, int, vector_t, int, int, int, real_t);