#include "math_base.h"
#include "number_theory.h"

uint64_t lagrange_interpolate(uint64_t x[], uint64_t y[], uint64_t n, uint64_t x0, uint64_t mod){
    uint64_t res = 0;
    uint64_t prod = 1;
    for (size_t i = 0; i < n; i++)
        prod = prod * sub(x0, x[i], mod) % mod;
    for (size_t i = 0; i < n; i++){
        uint64_t num = prod * inv(sub(x0, x[i], mod), mod) % mod;
        uint64_t den = 1;
        for (size_t j = 0; j < n; j++){
            if (i == j)
                continue;
            den = den * sub(x[i], x[j], mod) % mod;
        }
        res = (res + y[i] * num % mod * inv(den, mod)) % mod;
    }
    return res;
}

uint64_t lagrange_interpolate2(uint64_t x[], uint64_t y[], uint64_t* z[], uint64_t n, uint64_t x0, uint64_t y0, uint64_t mod){
    uint64_t res = 0;
    for (uint64_t i = 0; i < n; i++){
        for (uint64_t j = 0; j < n; j++){
            uint64_t t = z[i][j];
            for (uint64_t k = 0; k < n; k++){
                if (k == i)
                    continue;
                t = t * sub(x0, x[k], mod) % mod * inv(sub(x[i], x[k], mod), mod) % mod;
            }
            for (uint64_t k = 0; k < n; k++){
                if (k == j)
                    continue;
                t = t * sub(y0, y[k], mod) % mod * inv(sub(y[i], y[k], mod), mod) % mod;
            }
            res = (res + t) % mod;
        }
    }
    return res;
}