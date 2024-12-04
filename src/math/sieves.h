#ifndef SIEVES_H
#define SIEVES_H

#include "math_base.h"

#include <cstdint>
#include <unordered_map>
#include <functional>

// [0, n)
uint32_t euler_sieve(uint32_t n, uint32_t fac[], uint32_t prime[]);
uint64_t euler_sieve(uint64_t n, uint64_t fac[], uint64_t prime[]);
uint32_t egypt_sieve(uint32_t n, bool tag[], uint32_t prime[]);
uint64_t egypt_sieve(uint64_t n, bool tag[], uint64_t prime[]);

int64_t dujiao_sum(uint64_t x, std::unordered_map<uint64_t, int64_t>& cache, std::function<int64_t(uint64_t)> hsum);
int64_t dujiao_sieve(uint64_t n, std::function<int64_t(uint64_t)> hsum);

#endif // !SIEVES_H