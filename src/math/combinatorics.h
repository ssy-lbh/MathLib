#ifndef COMBINATORICS_H
#define COMBINATORICS_H

#include "math_base.h"

#include "number_theory.h"

uint64_t factorial(uint64_t n, uint64_t mod);
void factorial(uint64_t vals[], uint64_t n, uint64_t mod);
uint64_t inv_factorial(uint64_t invs[], uint64_t n, uint64_t mod);
uint64_t inv_factorial(uint64_t n, uint64_t mod);
void inv_factorial(uint64_t invs[], uint64_t vals[], uint64_t n, uint64_t mod);
uint64_t binomial_coefficient(uint64_t n, uint64_t k, uint64_t invs[], uint64_t mod);
uint64_t binomial_coefficient(uint64_t n, uint64_t k, uint64_t mod);
void binomial_coefficient(uint64_t vals[], uint64_t invs[], uint64_t n, uint64_t k, uint64_t mod);
uint64_t stirling_number(uint64_t n, uint64_t k, uint64_t mod);
uint64_t catalan_number(uint64_t n, uint64_t mod);
uint64_t bell_number(uint64_t n, uint64_t mod);
uint64_t cantor_number(uint64_t vals[], uint64_t n, uint64_t mod);

#endif /* COMBINATORICS_H */