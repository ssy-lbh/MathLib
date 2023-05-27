#ifndef SIEVES_H
#define SIEVES_H

#include "math_base.h"

#include <unordered_map>

// [0, n)
uint32_t euler_sieve(uint32_t n, uint32_t fac[], uint32_t prime[]);
uint64_t euler_sieve(uint64_t n, uint64_t fac[], uint64_t prime[]);
uint32_t egypt_sieve(uint32_t n, bool tag[], uint32_t prime[]);
uint64_t egypt_sieve(uint64_t n, bool tag[], uint64_t prime[]);

int64_t dujiao_sum(int64_t x, std::unordered_map<int64_t, int64_t>& cache, int64_t(*h)(int64_t));
int64_t dujiao_sieve(int64_t n, int64_t(*h)(int64_t));

#endif // !SIEVES_H