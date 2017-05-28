#ifndef FUNCS_H
#define FUNCS_H

#include <stdint.h>
#include <consts.h>

uint64_t factorial(int);
int64_t binomial(int, int);

ldouble_t Yslm(int, int, int, ldouble_t);
ldouble_t YslmPrime(int, int, int, ldouble_t);

ldouble_t BClimYp(int, int, int);
ldouble_t BClimYm(int, int, int);

#endif // FUNCS_H
