#include "number_theory.h"

uint32_t euler_sieve(uint32_t n, uint32_t fac[], uint32_t prime[]){
    uint32_t cnt = 0;
    for (uint32_t i = 2; i <= n; i++){
        if (!fac[i]){
            fac[i] = i;
            prime[cnt++] = i;
        }
        for (uint32_t j = 0; j < cnt; j++){
            if (i * prime[j] >= n || prime[j] > fac[i]) break;
            fac[i * prime[j]] = prime[j];
        }
    }
    return cnt;
}

uint64_t euler_sieve(uint64_t n, uint64_t fac[], uint64_t prime[]){
    uint64_t cnt = 0;
    for (uint64_t i = 2; i <= n; i++){
        if (!fac[i]){
            fac[i] = i;
            prime[cnt++] = i;
        }
        for (uint64_t j = 0; j < cnt; j++){
            if (i * prime[j] >= n || prime[j] > fac[i]) break;
            fac[i * prime[j]] = prime[j];
        }
    }
    return cnt;
}

uint32_t egypt_sieve(uint32_t n, bool tag[], uint32_t prime[]){
    uint32_t cnt = 0;
    for (uint32_t i = 2; i <= n; i++){
        if (!tag[i]){
            prime[cnt++] = i;
            for (uint64_t j = i * i; j < n; j += i)
                tag[j] = true;
        }
    }
    return cnt;
}

uint64_t egypt_sieve(uint64_t n, bool tag[], uint64_t prime[]){
    uint64_t cnt = 0;
    for (uint64_t i = 2; i <= n; i++){
        if (!tag[i]){
            prime[cnt++] = i;
            if (i >= (n + i - 1) / i)
                continue;
            for (uint64_t j = i * i; j < n; j += i)
                tag[j] = true;
        }
    }
    return cnt;
}