#include "utils.h"
#include <math.h>

char convert(ldouble_t x)
{
    if (fabsl(x) < 1e-8)
        return '0';
    if (fabsl(x - 1.0L) < 1e-8)
        return '1';
    return 'x';
}

void print_smatrix_block(matrix_t s, size_t spaces)
{
    for (size_t i = 0; i < s->size1; i++)
    {
        printf("%c | ", convert(matrix_get(s, i, s->size2 - 1)));
        for (size_t j = 0; j < spaces; j++)
            printf("  ");
        for (size_t j = 0; j < s->size2 - 1; j++)
            printf("%c ", convert(matrix_get(s, i, j)));
        printf("\n");
    }
}

void output_vector_h(FILE *fp, vector_t Y)
{
    for (size_t i = 0; i < Y->size - 1; i++)
        fprintf(fp, "%Le, ", vector_get(Y, i));

    fprintf(fp, "%Le\n", vector_get(Y, Y->size - 1));
}

void output_matrix(FILE *fp, matrix_t Y)
{
    if (fp == NULL)
        GSL_ERROR_VOID("Stream is NULL, in  output_block_matrix", GSL_EFAILED);

    for (size_t j = 0; j < Y->size2; j++)
    {
        for (size_t i = 0; i < Y->size1 - 1; i++)
            fprintf(fp, "%Le, ", matrix_get(Y, i, j));

        fprintf(fp, "%Le\n", matrix_get(Y, Y->size1 - 1, j));
    }
}

void output_block_matrix(FILE *fp, matrix_t *S, size_t len, size_t col_shift)
{
    if (fp == NULL)
        GSL_ERROR_VOID("Stream is NULL, in  output_block_matrix", GSL_EFAILED);

    size_t rows = 0, cols = 0;

    for (size_t k = 0; k < len; k++)
    {
        rows += S[k]->size1;
        cols += (S[k]->size2 - 1) - col_shift;
    }
    cols -= col_shift;

    // print S[0]
    if (len > 0)
    {
        for (size_t i = 0; i < S[0]->size1; i++)
        {
            for (size_t j = col_shift; j < S[0]->size2 - 1; j++)
                fprintf(fp, "%+Le, ", matrix_get(S[0], i, j));

            for (size_t j = 0; j < cols - S[0]->size2 + 1 + col_shift; j++)
                fprintf(fp, "%+Le, ", 0.0L);

            fprintf(fp, "%+Le\n", matrix_get(S[0], i, S[0]->size2 - 1));
        }
    }

    for (size_t k = 1; k < len - 1; k++)
    {
        for (size_t i = 0; i < S[k]->size1; i++)
        {
            for (size_t j = 0; j < (k - 1) * col_shift; j++)
                fprintf(fp, "%+Le, ", 0.0L);

            for (size_t j = 0; j < S[k]->size2 - 1; j++)
                fprintf(fp, "%+Le, ", matrix_get(S[k], i, j));

            for (size_t j = 0; j < cols - S[k]->size2 + 1 - (k - 1) * col_shift; j++)
                fprintf(fp, "%+Le, ", 0.0L);

            fprintf(fp, "%+Le\n", matrix_get(S[k], i, S[k]->size2 - 1));
        }
    }

    // print S[len-1]
    if (len > 1)
    {
        for (size_t i = 0; i < S[len - 1]->size1; i++)
        {
            for (size_t j = 0; j < cols - S[len - 1]->size2 + 1 + col_shift; j++)
                fprintf(fp, "%+Le, ", 0.0L);

            for (size_t j = 0; j < S[len - 1]->size2 - 1 - col_shift; j++)
                fprintf(fp, "%+Le, ", matrix_get(S[len - 1], i, j));

            fprintf(fp, "%+Le\n", matrix_get(S[len - 1], i, S[len - 1]->size2 - 1));
        }
    }
}

void output_block_matrix_mathematica(FILE *fp, matrix_t *S, size_t len, size_t col_shift)
{
    size_t rows = 0, cols = 0;

    for (size_t k = 0; k < len; k++)
    {
        rows += S[k]->size1;
        cols += (S[k]->size2 - 1) - col_shift;
    }
    cols += -col_shift;

    size_t r = 1, c = 1;

    // print S[0]
    if (len > 0)
    {
        for (size_t i = 0; i < S[0]->size1; i++)
        {
            for (size_t j = col_shift; j < S[0]->size2 - 1; j++)
                fprintf(fp, "%lu, %lu, %Le\n", r, c++, matrix_get(S[0], i, j));

            fprintf(fp, "%lu, %lu, %Le\n", r, cols + 1, matrix_get(S[0], i, S[0]->size2 - 1));

            r++;
            c -= S[0]->size2 - 1 - col_shift;
        }
    }

    for (size_t k = 1; k < len - 1; k++)
    {
        for (size_t i = 0; i < S[k]->size1; i++)
        {
            for (size_t j = 0; j < S[k]->size2 - 1; j++)
                fprintf(fp, "%lu, %lu, %Le\n", r, c++, matrix_get(S[k], i, j));

            fprintf(fp, "%lu, %lu, %Le\n", r, cols + 1, matrix_get(S[k], i, S[k]->size2 - 1));

            r++;
            c -= S[k]->size2 - 1;
        }
    }

    // print S[len-1]
    if (len > 1)
    {
        for (size_t i = 0; i < S[len - 1]->size1; i++)
        {
            for (size_t j = 0; j < S[len - 1]->size2 - 1 - col_shift; j++)
                fprintf(fp, "%lu, %lu, %Le\n", r, c++, matrix_get(S[len - 1], i, j));

            fprintf(fp, "%lu, %lu, %Le\n", r, cols + 1, matrix_get(S[len - 1], i, S[len - 1]->size2 - 1));

            r++;
            c -= S[len - 1]->size2 - 1 - col_shift;
        }
    }
}
