#include "number_theory.h"

void inv(uint32_t vals[], uint32_t n, uint32_t mod){
    vals[0] = 0;
    vals[1] = 1;
    for (uint32_t i = 2; i < n; i++)
        vals[i] = mul(mod - mod / i, vals[mod % i], mod);
}

void inv(uint64_t vals[], uint64_t n, uint64_t mod){
    vals[0] = 0;
    vals[1] = 1;
    for (uint64_t i = 2; i < n; i++)
        vals[i] = mul(mod - mod / i, vals[mod % i], mod);
}
