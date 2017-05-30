#include <consts.h>
#include <method.h>
#include <swsh.h>

#define N_PTS 401
#define STEP ((XF - XI) / (N_PTS - 1))

// E_k definition that gives "closeness" to a solution
real_t Ek(int k, int j, matrix_t Y, int s, int l, int m, real_t c)
{
    int r;
    real_t xk, result;
    vector_t yk;

    // Trapezoid method for x and y
    xk = (k - 0.5) * STEP + X_I; 
    yk = VECTOR_ALLOC(N_EQS);
    for (r = 0; r < N_EQS; r++)
        VECTOR_SET(yk, r, (MATRIX_GET(Y, r, k - 1) + MATRIX_GET(Y, r, k)) / 2.0);   

    // Internal mesh points

    return
}

// Concatenated matrix of the jacobian of E_k relative to y_{k-1} and y_{k}
real_t Sk(int k, int j, int i, matrix_t Y, int s, int l, int m, real_t c)
{
    int r;
    real_t xk, result;
    vector_t yk;

    // Trapezoid method for x and y
    xk = (k - 0.5) * STEP + X_I; 
    yk = VECTOR_ALLOC(N_EQS);
    for (r = 0; r < N_EQS; r++)
        VECTOR_SET(yk, r, (MATRIX_GET(Y, r, k - 1) + MATRIX_GET(Y, r, k)) / 2.0);    

    // Internal mesh points
    result = -0.5 * STEP * Dg(j, i, yk, xk, s, l, m, c);

    if (j == i) // derivative relative to y_{k-1} 
        result += -1.0;
    else if (j == i - N_EQS) // derivative relative to y_{k} 
        result += 1.0;
    
    return result;
}