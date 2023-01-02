#include <math.h>
#include <stdint.h>

#include <gsl/gsl_errno.h>

#include "const.h"
#include "utils.h"
#include "swsh.h"

bool possibly_zero(ldouble_t x)
{
    return fabsl(x) < ZERO_EPS;
}

bool possibly_same(ldouble_t x, ldouble_t y)
{
    if (possibly_zero(x))
        return possibly_zero(y);

    return fabsl(x / y - 1) < SAME_REL_EPS;
}

uint64_t factorial(int n)
{
    if (n == 0)
        return 1;
    else if (n < 0)
        GSL_ERROR("\nfactorial function input domain error.\n", GSL_EDOM);
    else
    {
        uint64_t i, p = 1;
        for (i = 1; i <= n; i++)
            p *= i;
        return p;
    }
}

int64_t binomial(int n, int r)
{
    if (r == 0 || r == n)
        return 1;

    if (n < 0)
    {
        if (r >= 0)
            return (2 * ((r + 1) % 2) - 1) * binomial(-n + r - 1, r);
        else if (r <= n)
            return (2 * ((n - r + 1) % 2) - 1) * binomial(-r - 1, n - r);
        else
            return 0;
    }
    else
    {
        if (r < 0 || r > n)
            return 0;
        else
            return factorial(n) / factorial(r) / factorial(n - r);
    }
}

ldouble_t powl_t(ldouble_t x, ldouble_t p)
{
    if (possibly_zero(p))
        return 1.0L;
    return powl(x, p);
}

ldouble_t sgn_t(int p)
{
    return (p % 2 == 0) ? 1.0L : -1.0L;
}

// Function (1-x)^-kp (1+x)^-km Yslm(x)
ldouble_t Yslm(int s, int l, int m, ldouble_t x)
{
    int rMax = MIN(l - s, l + m);
    int rMin = MAX(m - s, 0);

    if (rMax < rMin)
        GSL_ERROR("\nl < max(abs(s),abs(m)), in YslmPrime.\n", GSL_EDOM);

    ldouble_t common_sq = (1.0L + 2.0L * l) * factorial(l - m) * factorial(l + m);
    common_sq /= powl_t(4.0L, 1.0L + l) * PI * factorial(l - s) * factorial(l + s);

    ldouble_t sum = 0.0L;
    for (int r = rMin; r <= rMax; r++)
        sum += powl_t(-1.0L + x, -r + rMax) * powl_t(1.0L + x, r - rMin) * binomial(l - s, r) * binomial(l + s, -m + r + s);

    return sgn_t(l - s + m + rMax) * sqrtl(common_sq) * sum;
}

// Limit of (1-x)^-kp (1+x)^-km Yslm(x) at x=-+1
ldouble_t lim_Yslm(int s, int l, int m, int sgn)
{
    int rMax = MIN(l - s, l + m);
    int rMin = MAX(m - s, 0);

    if (rMax < rMin)
        GSL_ERROR("\n l < max(abs(s),abs(m)), in YslmPrime.\n", GSL_EDOM);

    ldouble_t sum, common_sq, rest;
    common_sq = (1.0L + 2.0L * l) * factorial(l - m) * factorial(l + m);
    common_sq /= powl_t(4.0L, 1.0L + l + rMin - rMax) * PI * factorial(l - s) * factorial(l + s);

    int r = sgn < 0 ? rMin : rMax;
    rest = binomial(l - s, r) * binomial(l + s, -m + r + s);
    return sgn_t(l + s + m + r) * sqrtl(common_sq) * rest;
}

// Function d/dx( (1-x)^-kp (1+x)^-km Yslm(x) )
ldouble_t YslmPrime(int s, int l, int m, ldouble_t x)
{
    int r;
    int rMax = MIN(l - s, l + m);
    int rMin = MAX(m - s, 0);
    int h = rMax - rMin;

    if (h < 0)
        GSL_ERROR("\n l < max(abs(s),abs(m)), in YslmPrime.\n", GSL_EDOM);

    if (h == 0)
        return 0.0L;

    ldouble_t common_sq = (1.0L + 2.0L * l) * factorial(l - m) * factorial(l + m);
    common_sq /= powl_t(4.0L, 1.0L + l) * PI * factorial(l - s) * factorial(l + s);

    ldouble_t sum = 0.0L;

    for (r = rMin; r <= rMax - 1; r++)
    {
        sum += (1.0L + r - rMin) * powl_t(-1.0L + x, -1 - r + rMax) * powl_t(1.0L + x, r - rMin) * binomial(l - s, 1 + r) * binomial(l + s, 1 - m + r + s);
        sum += (-r + rMax) * powl_t(-1.0L + x, -1 - r + rMax) * powl_t(1.0L + x, r - rMin) * binomial(l - s, r) * binomial(l + s, -m + r + s);
    }

    return sgn_t(l - s + m + rMax) * sqrtl(common_sq) * sum;
}

ldouble_t lim_YslmPrime(int s, int l, int m, int sgn)
{
    int rMax = MIN(l - s, l + m);
    int rMin = MAX(m - s, 0);

    if (rMax < rMin)
        GSL_ERROR("\n l < max(abs(s),abs(m)), in YslmPrime.\n", GSL_EDOM);

    if (rMax == rMin)
        return 0.0L;

    ldouble_t common_sq;
    common_sq = (1.0L + 2.0L * l) * factorial(l - m) * factorial(l + m);
    common_sq /= powl_t(4.0L, 2.0L + l - rMax + rMin) * PI * factorial(l - s) * factorial(l + s);

    int r = sgn < 0 ? rMax - 1 : rMin;
    int64_t res = (rMax - r) * binomial(l - s, r) * binomial(l + s, -m + r + s) + (1 + r - rMin) * binomial(l - s, 1 + r) * binomial(l + s, 1 - m + r + s);

    return sgn_t(l - s + m + r + 1) * sqrtl(common_sq) * res;
}

// Ratio of Zslm'(x)/Zslm(x) at x=+1
ldouble_t limratio_prime_Yslm(int s, int l, int m, ldouble_t A, ldouble_t c, int sgn)
{
    if (l < MAX(abs(s), abs(m)))
        GSL_ERROR("\nBCratioZp function s,l,m input error.\n", GSL_EDOM);

    ldouble_t kp, km;
    kp = abs(m + s) / 2.0L;
    km = abs(m - s) / 2.0L;

    int sign = sgn < 0 ? -1 : 1;

    ldouble_t res = A + c * c - 2.0L * sign * c * s - (km + kp - s) * (1.0L + km + kp + s);
    res /= (sgn < 0 ? (2.0L + 4.0L * km) : (2.0L + 4.0L * kp));

    return sign * res;
}

// g contains all information about dif. equation: yi'=gi(x,y0,y1,..)
ldouble_t g(int i, vector_t y, ldouble_t x, int s, int l, int m, ldouble_t c)
{
    if (i < 0 || i >= N_EQS || l < MAX(abs(s), abs(m)))
        GSL_ERROR("\ng function input error.\n", GSL_EDOM);

    ldouble_t result, kp, km;
    kp = abs(m + s) / 2.0L;
    km = abs(m - s) / 2.0L;

    if (i == 0)
        result = vector_get(y, 1);
    else if (i == 1)
    {
        result = 2.0L * (kp - km + (1.0L + km + kp) * x) * vector_get(y, 1) + vector_get(y, 0) * ((km + kp) * (1.0L + km + kp) - c * c * x * x + 2 * c * s * x - vector_get(y, 2) - s * (s + 1.0L));
        result /= (1.0L - x * x);
    }
    else if (i == 2)
        result = 0.0L;

    return result;
}

// jacobian of g
ldouble_t Dg(int i, int j, vector_t y, ldouble_t x, int s, int l, int m, ldouble_t c)
{
    if (j < 0 || j >= N_EQS || i < 0 || i >= N_EQS || l < MAX(abs(s), abs(m)))
        GSL_ERROR("\ng function input error.\n", GSL_EDOM);

    ldouble_t result, kp, km;
    kp = abs(m + s) / 2.0L;
    km = abs(m - s) / 2.0L;

    if (i == 0)
    {
        if (j == 1)
            result = 1.0L;
        else
            result = 0.0L;
    }
    if (i == 1)
    {
        result = 0.0L;

        if (j == 0)
            result = (km + kp) * (1.0L + km + kp) - c * c * x * x + 2 * c * s * x - vector_get(y, 2) - s * (s + 1.0L);
        if (j == 1)
            result = 2.0L * (kp - km + (1.0L + km + kp) * x);
        if (j == 2)
            result = -vector_get(y, 0);

        result /= (1.0L - x * x);
    }
    if (i == 2)
        result = 0.0L;

    return result;
}

// E_k definition that gives "closeness" to a solution
ldouble_t Ek(int k, int i, vector_t X, matrix_t Y, int s, int l, int m, ldouble_t c)
{
    size_t N = Y->size1;
    size_t M = Y->size2;
    // boundary conditions at k=M (x=+1)
    if (k == M)
    {
        if (i == 0)
            return matrix_get(Y, 1, M - 1) - matrix_get(Y, 0, M - 1) * limratio_prime_Yslm(s, l, m, matrix_get(Y, 2, M - 1), c, +1);
        else if (i == 1)
            return matrix_get(Y, 0, M - 1) - lim_Yslm(s, l, m, +1);
        else
            return 0.0L;
    }

    // boundary conditions at k=0 (x=-1)
    if (k == 0)
    {
        if (i == 2)
            return matrix_get(Y, 1, 0) - matrix_get(Y, 0, 0) * limratio_prime_Yslm(s, l, m, matrix_get(Y, 2, 0), c, -1);
        else
            return 0.0L;
    }

    // boundary conditions at k=0 (x=-1)

    ldouble_t xk = (vector_get(X, k - 1) + vector_get(X, k)) / 2.0L;
    vector_t yk = vector_alloc(Y->size1);
    for (int r = 0; r < N; r++)
        vector_set(yk, r, (matrix_get(Y, r, k - 1) + matrix_get(Y, r, k)) / 2.0L);

    ldouble_t result = matrix_get(Y, i, k) - matrix_get(Y, i, k - 1) - STEP * g(i, yk, xk, s, l, m, c);

    vector_free(yk);

    return result;
}

// Concatenated matrix of the jacobian of E_k relative to y_{k-1} and y_{k}
ldouble_t Sk(int k, int i, int j, vector_t X, matrix_t Y, int s, int l, int m, ldouble_t c)
{
    size_t N = Y->size1;
    size_t M = Y->size2;

    ldouble_t kp, km;
    kp = abs(m + s) / 2.0L;
    km = abs(m - s) / 2.0L;

    // boundary conditions at k=M-1 (x=+1)
    if (k == M)
    {
        if (i == 0 && j == 0)
            return -limratio_prime_Yslm(s, l, m, matrix_get(Y, 2, M - 1), c, +1);
        if (i == 0 && j == 1)
            return 1.0L;
        if (i == 0 && j == 2)
            return matrix_get(Y, 0, M - 1) / (2.0L + 4.0L * kp);
        if (i == 1 && j == 0)
            return 1.0L;

        return 0.0L;
    }

    // boundary conditions at k=0 (x=-1)
    if (k == 0)
    {
        if (i == 2 && j == N + 0)
            return -limratio_prime_Yslm(s, l, m, matrix_get(Y, 2, 0), c, -1);
        if (i == 2 && j == N + 1)
            return 1.0L;
        if (i == 2 && j == N + 2)
            return matrix_get(Y, 0, 0) / (2.0L + 4.0L * km);

        return 0.0L;
    }

    // internal mesh points k=1,2,..,M-2
    ldouble_t xk = (vector_get(X, k - 1) + vector_get(X, k)) / 2.0L;
    vector_t yk = vector_alloc(Y->size1);
    for (int r = 0; r < N; r++)
        vector_set(yk, r, (matrix_get(Y, r, k - 1) + matrix_get(Y, r, k)) / 2.0L);

    // Internal mesh points
    ldouble_t result = -0.5L * STEP * Dg(i, j % N, yk, xk, s, l, m, c);

    if (i == j) // derivative relative to y_{k-1}
        result += -1.0L;
    if (i == j - N) // derivative relative to y_{k}
        result += 1.0L;

    vector_free(yk);

    return result;
}

void smatrix_fill(matrix_t * smatrix, vector_t X, matrix_t Y, int s, int l, int m, ldouble_t c)
{
    size_t N = Y->size1;
    size_t M = Y->size2;

    size_t bc_i = BC_I, bc_f = BC_F;

    for(size_t k = 0; k <= M; k++)
    {
        size_t n1 = N, di = 0;
        if (k == 0)
            n1 = bc_i, di = bc_f;
        if (k == M)
            n1 = bc_f;

        smatrix[k] = matrix_alloc(n1, 2 * N + 1);

        for (int i = 0; i < n1; i++)
        {
            for (int j = 0; j < 2 * N; j++)
                matrix_set(smatrix[k], i, j, Sk(k, di + i, j, X, Y, s, l, m, c));

            matrix_set(smatrix[k], i, 2 * N, -Ek(k, di + i, X, Y, s, l, m, c));
        }
    }
}

void block_gauss_pivot(matrix_t s)
{
    size_t row = 0, col = 0;

    while (row < s->size1 && col < s->size2)
    {
        size_t idx_max = row;
        ldouble_t piv = matrix_get(s, row, col);
        ldouble_t piv_abs = fabsl(piv);

        for (size_t i = row; i < s->size1; i++)
        {
            ldouble_t tmp = matrix_get(s, i, col);
            ldouble_t tmp_abs = fabsl(tmp);
            if (tmp_abs > piv_abs)
                idx_max = i, piv = tmp, piv_abs = tmp_abs;
        }

        if (possibly_zero(piv))
        {
            col++;
            continue;
        }

        if (row != idx_max)
            matrix_swap_rows(s, row, idx_max);

        for (size_t j = 0; j < s->size2; j++)
            matrix_set(s, row, j, matrix_get(s, row, j) / piv);

        for (size_t i = row + 1; i < s->size1; i++)
        {
            ldouble_t f = matrix_get(s, i, col);

            for (size_t j = 0; j < s->size2; j++)
                matrix_set(s, i, j, matrix_get(s, i, j) - f * matrix_get(s, row, j));

            matrix_set(s, i, col, 0.0L);
        }

        row++;
        col++;
    }
}

void block_reduce(matrix_t s, matrix_t s_top, size_t col_shift)
{
    if (col_shift >= s_top->size2 - 1)
        GSL_ERROR_VOID("col_shift >= sview_top->size2 - 1, in block_reduce", GSL_EBADLEN);

    // find lowest row with 0's next to a 1.0 after col > col_shift
    size_t r = 0, c = 0;
    for (r = 0; r < s_top->size1; r++)
    {
        u_int32_t flag = 1;
        for (c = 0; c < col_shift; c++)
            flag &= possibly_zero(matrix_get(s_top, r, c));

        if (flag)
            break;
    }

    c = col_shift;
    while (r < s_top->size1 && c < s_top->size2 - 1)
    {
        ldouble_t piv_top = matrix_get(s_top, r, c);

        // find column with non-zero entry (=1 after gauss elimination)
        if (possibly_zero(piv_top))
        {
            c++;
            continue;
        }

        for (size_t i = 0; i < s->size1; i++)
        {
            ldouble_t piv = matrix_get(s, i, c - col_shift);

            // skip if element in the same column is already zero
            if (possibly_zero(piv))
                continue;

            // row subtraction on non-zero columns from both blocks
            for (size_t j = 0; j < s->size2 - 1 && j + col_shift < s_top->size2 - 1; j++)
            {
                ldouble_t x = matrix_get(s, i, j) - piv / piv_top * matrix_get(s_top, r, j + col_shift);
                matrix_set(s, i, j, x);
            }
            matrix_set(s, i, c - col_shift, 0.0L);

            // update last column, with the b vector
            ldouble_t bi = matrix_get(s, i, s->size2 - 1) - piv / piv_top * matrix_get(s_top, r, s_top->size2 - 1);
            matrix_set(s, i, s->size2 - 1, bi);
        }

        r++;
    }
}

void block_back_substituition(matrix_t *s, size_t len, size_t col_shift)
{
    size_t rows = 0;
    for (size_t k = 0; k < len; k++)
        rows += s[k]->size1;

    vector_t X = vector_alloc(rows);

    size_t r = rows - 1;

    for (int k = len - 1; k >= 0; k--)
    {
        size_t col_cut = (k == len - 1) ? col_shift : 0;

        // find col of first non-zero in the last row
        size_t c_piv;
        ldouble_t piv = 0.0L;
        for (c_piv = 0; c_piv < s[k]->size2 - 1; c_piv++)
        {
            piv = matrix_get(s[k], s[k]->size1 - 1, c_piv);
            if (!possibly_zero(piv))
                break;
        }

        for (int i = s[k]->size1 - 1; i >= 0; i--)
        {
            ldouble_t x = matrix_get(s[k], i, s[k]->size2 - 1);

            size_t j_len = s[k]->size2 > 1 + c_piv + col_cut ? s[k]->size2 - 1 - c_piv - col_cut : 0; 
             
            for (size_t j = 1; j < j_len; j++)
            {
                x += -matrix_get(s[k], i, c_piv + j) * vector_get(X, r + j);
                matrix_set(s[k], i, c_piv + j, 0.0L);
            }
            vector_set(X, r, x / piv);
            matrix_set(s[k], i, s[k]->size2 - 1, x / piv);

            r--;
            c_piv--;
        }
    }

    vector_free(X);
}

void generate_initial_Y(matrix_t Y, vector_t X, int s, int l, int m, ldouble_t *err_scale)
{
    size_t N = Y->size1;
    size_t M = Y->size2;

    for (int k = 1; k < M - 1; k++)
    {
        ldouble_t xk = vector_get(X, k);
        ldouble_t y0k = Yslm(s, l, m, xk);
        ldouble_t y1k = YslmPrime(s, l, m, xk);
        ldouble_t y2k = l * (l + 1.0L) - s * (s + 1.0L);

        matrix_set(Y, 0, k, y0k);
        matrix_set(Y, 1, k, y1k);
        matrix_set(Y, 2, k, y2k);
    }

    matrix_set(Y, 0, 0, lim_Yslm(s, l, m, -1));
    matrix_set(Y, 1, 0, lim_YslmPrime(s, l, m, -1));
    matrix_set(Y, 2, 0, l * (l + 1.0L) - s * (s + 1.0L));

    matrix_set(Y, 0, M - 1, lim_Yslm(s, l, m, +1));
    matrix_set(Y, 1, M - 1, lim_YslmPrime(s, l, m, +1));
    matrix_set(Y, 2, M - 1, l * (l + 1.0L) - s * (s + 1.0L));

    err_scale[0] = MAX(matrix_get(Y, 0, 0), matrix_get(Y, 0, M - 1));
    err_scale[1] = MAX(matrix_get(Y, 1, 0), matrix_get(Y, 1, M - 1));
    err_scale[1] = MAX(err_scale[0], err_scale[1]);
    err_scale[2] = 1.0L;

    for(size_t i = 0; i < N; i++)
        if(possibly_zero(err_scale[i]))
            err_scale[i] = 1.0L;
}

void desolve_relaxation(vector_t X, matrix_t Y, int s, int l, int m, ldouble_t c, ldouble_t *err_scale, u_int32_t log_flag)
{
    size_t N = Y->size1;
    size_t M = Y->size2;

    size_t iter = 0;
    ldouble_t err = 2 * CONV_EPS + 1.0L;

    while (err > CONV_EPS)
    {
        if (iter >= ITER_MAX)
            GSL_ERROR_VOID("Reached ITER_MAX iteration, in desolve_relaxation", GSL_EMAXITER);

        matrix_t *smatrix = (matrix_t *)malloc((M + 1) * sizeof(matrix_t));
        if (smatrix == NULL)
            GSL_ERROR_VOID("failed to malloc matrix_t vector, in desolve_relaxation", GSL_ENOMEM);

        smatrix_fill(smatrix, X, Y, s, l, m, c);

        if((log_flag & LOG_SMATRIX) == LOG_SMATRIX)
        {
            char fname0[50];
            sprintf(fname0, "smat_%d_%d_%d_%Lf_%d.csv", s, l, m, iter);
            FILE *fp0;
            fp0 = fopen(fname0, "w");
            output_block_matrix(fp0, smatrix, M + 1, N);
            fclose(fp0);
        }

        block_gauss_pivot(smatrix[0]);
        for (size_t k = 1; k <= M; k++)
        {
            block_reduce(smatrix[k], smatrix[k - 1], N);
            block_gauss_pivot(smatrix[k]);
        }

        if((log_flag & LOG_SMATRIX_GAUSS) == LOG_SMATRIX_GAUSS)
        {
            char fname1[50];
            sprintf(fname1, "smatg_%d_%d_%d_%Lf_%d.csv", s, l, m, iter);
            FILE *fp1;
            fp1 = fopen(fname1, "w");
            output_block_matrix(fp1, smatrix, M + 1, N);
            fclose(fp1);
        }

        block_back_substituition(smatrix, M + 1, N);

        matrix_t delta_Y = matrix_alloc(Y->size1, Y->size2);
        size_t row = 0;
        for (int k = 0; k <= M; k++)
        {
            for (size_t i = 0; i < smatrix[k]->size1; i++)
            {
                ldouble_t dy = matrix_get(smatrix[k], i, smatrix[k]->size2 - 1);
                matrix_set(delta_Y, row % delta_Y->size1, row / delta_Y->size1, dy);
                row++;
            }
        }

        if((log_flag & LOG_DELTA_Y) == LOG_DELTA_Y)
        {
            char fname2[50];
            sprintf(fname2, "log/dy_%d_%d_%d_%Lf_%d.csv", s, l, m, c, iter);
            FILE *fp2;
            fp2 = fopen(fname2, "w");
            output_matrix(fp2, delta_Y);
            fclose(fp2);
        }

        err = 0.0L;
        for (size_t i = 0; i < Y->size1; i++)
        {
            ldouble_t err_i = 0.0L;
            for (size_t j = 0; j < Y->size2; j++)
                err_i += fabsl(matrix_get(delta_Y, i, j));
            err += err_i / err_scale[i];
        }
        err /= 1.0L * N * M;

        if((log_flag & PRINT_ERROR) == PRINT_ERROR)
            printf("(s, l, m, c): (%d, %d, %d, %Lf), err(%lu): %Lf \n", s, l, m, c, iter, err);

        ldouble_t factor = (err > CONTROL_ERR) ? (CONTROL_ERR / err) : 1.0;
        matrix_scale(delta_Y, factor);
        matrix_add(Y, delta_Y);

        for (size_t k = 0; k <= M; k++)
            matrix_free(smatrix[k]);
        free(smatrix);

        matrix_free(delta_Y);

        iter++;
    }
}