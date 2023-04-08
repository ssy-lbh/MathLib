#include "number_theory.h"

constexpr uint32_t N = 16;

uint64_t totient(uint64_t n, uint64_t prime[], uint32_t cnt){
    for (uint32_t i = 0; i < cnt; i++)
        n = n / prime[i] * (prime[i] - 1);
    return n;
}

uint64_t totient(uint64_t n){
    uint64_t prime[N];
    uint32_t exp[N];

    uint32_t cnt = pollard_rho(n, prime, exp, N);

    return totient(n, prime, cnt);
}
