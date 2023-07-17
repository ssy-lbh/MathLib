#include "number_theory.h"

uint64_t crt_calc(uint64_t b[], uint64_t m[], uint32_t n){
    uint64_t M = 1, ans = 0;
    for (uint32_t i = 0; i < n; i++)
        M *= m[i];
    for (uint32_t i = 0; i < n; i++){
        uint64_t w = M / m[i];
        ans = (ans + b[i] * inv(w, m[i]) % M * w % M) % M;
    }
    return ans % M;
}

uint64_t crt_calc(uint64_t b[], uint64_t m[], uint32_t n, uint64_t mod){
    uint64_t ans = 0;
    for (uint32_t i = 0; i < n; i++){
        uint64_t x = mul(b[i] % mod, mod / m[i], mod);
        x = mul(x, inv(mod / m[i], m[i]) % mod, mod);
        ans = (ans + x) % mod;
    }
    return ans;
}