#ifndef DIFFEQ_H
#define DIFFEQ_H

#include <stdint.h>
#include <consts.h>

uint64_t factorial(int);
int64_t binomial(int, int);

real_t Yslm(int, int, int, real_t);
real_t YslmPrime(int, int, int, real_t);

real_t BClimYp(int, int, int);
real_t BClimYm(int, int, int);

real_t g(int, vector_t, real_t, int, int, int, real_t);
real_t Dg(int, int, vector_t, real_t, int, int, int, real_t);

#endif // DIFFEQ_H
