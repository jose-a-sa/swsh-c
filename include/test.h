#ifndef TEST_H
#define TEST_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

bool test_unit(bool, char *, bool);
void test_title(char *, bool);
void test_subtitle(char *, bool);
void test_print(bool, bool);

bool test_atomic(bool);
bool test_factorial(bool);
bool test_binomial(bool);
bool test_Yslm(bool);
bool test_BC(bool);

#endif // TEST_H