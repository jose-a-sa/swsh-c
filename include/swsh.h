#ifndef SWSH_H
#define SWSH_H

#include <stdint.h>
#include <consts.h>

uint64_t factorial(int);
int64_t binomial(int, int);

real_t Yslm(int, int, int, real_t);
real_t YslmPrime(int, int, int, real_t);

real_t BClimZp(int, int, int);
real_t BClimZm(int, int, int);

real_t BCratioZp(int, int, int, real_t, real_t);
real_t BCratioZm(int, int, int, real_t, real_t);

real_t g(int, vector_t, real_t, int, int, int, real_t);
real_t Dg(int, int, vector_t, real_t, int, int, int, real_t);

#endif // SWSH_H
