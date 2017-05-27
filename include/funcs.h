#ifndef FUNCS_H
#define FUNCS_H

#include <stdint.h>
#include <consts.h>


uint64_t factorial(int);
uint64_t binomial(int, int);

ldouble_t Yslm(int, int, int, ldouble_t);
ldouble_t YslmPrime(int, int, int, ldouble_t);

ldouble_t BCp(int32_t, int32_t, int32_t);
ldouble_t BCm(int32_t, int32_t, int32_t);


#endif // FUNCS_H
