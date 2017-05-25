#ifndef FUNCS_H
#define FUNCS_H

#include <stdint.h>
#include <utils.h>

#define PI 3.1415926535897932385L


uint64_t factorial(int32_t);
uint64_t binomial(int32_t, int32_t);
ldouble_t Yslm(int32_t, int32_t, int32_t, ldouble_t);
ldouble_t YslmPrime(int32_t, int32_t, int32_t, ldouble_t);


#endif // FUNCS_H
