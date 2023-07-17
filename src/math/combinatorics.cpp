#include "combinatorics.h"

uint64_t factorial(uint64_t n, uint64_t mod){
    uint64_t res = 1;
    for (uint64_t i = 2; i <= n; i++)
        res = mul(res, i, mod);
    return res;
}

void factorial(uint64_t vals[], uint64_t n, uint64_t mod){
    vals[0] = 1;
    for (uint64_t i = 1; i <= n; i++)
        vals[i] = mul(vals[i - 1], i, mod);
}

uint64_t inv_factorial(uint64_t invs[], uint64_t n, uint64_t mod){
    uint64_t res = 1;
    for (uint64_t i = 2; i <= n; i++)
        res = mul(res, invs[i], mod);
    return res;
}

uint64_t inv_factorial(uint64_t n, uint64_t mod){
    uint64_t* invs = new uint64_t[n + 1];
    inv(invs, n, mod);
    uint64_t res = inv_factorial(invs, n, mod);
    delete[] invs;
    return res;
}

void inv_factorial(uint64_t invs[], uint64_t vals[], uint64_t n, uint64_t mod){
    vals[0] = 1;
    for (uint64_t i = 1; i <= n; i++)
        vals[i] = mul(vals[i - 1], invs[i], mod);
}

uint64_t binomial_coefficient(uint64_t n, uint64_t k, uint64_t invs[], uint64_t mod){
    if (k > n)
        return 0;
    return mul(factorial(n, mod), mul(inv_factorial(invs, k, mod), inv_factorial(invs, n - k, mod), mod), mod);
}

uint64_t binomial_coefficient(uint64_t n, uint64_t k, uint64_t mod){
    if (k > n)
        return 0;
    uint64_t res = 1;
    for (uint64_t i = 1; i <= k; i++)
        res = mul(res, mul(n - i + 1, inv(i, mod), mod), mod);
    return res;
}

// [0..k] binomial coefficients
void binomial_coefficient(uint64_t vals[], uint64_t invs[], uint64_t n, uint64_t k, uint64_t mod){
    vals[0] = 1;
    for (uint64_t i = 1; i <= k; i++)
        vals[i] = mul(vals[i - 1], mul(n - i + 1, invs[i], mod), mod);
}

uint64_t stirling_number(uint64_t n, uint64_t k, uint64_t mod){
    uint64_t res = 0;
    for (uint64_t i = 0; i <= k; i++){
        uint64_t tmp = mul(pow(k - i, n, mod), inv(factorial(i, mod), mod), mod);
        if ((k - i) & 1)
            res = (res + mod - tmp) % mod;
        else
            res = (res + tmp) % mod;
    }
    return mul(res, inv(factorial(k, mod), mod), mod);
}

uint64_t catalan_number(uint64_t n, uint64_t mod){
    return mul(binomial_coefficient(n << 1, n, mod), inv(n + 1, mod), mod);
}

uint64_t bell_number(uint64_t n, uint64_t mod){
    uint64_t res = 0;
    for (uint64_t i = 0; i <= n; i++)
        res = (res + mul(stirling_number(n, i, mod), factorial(i, mod), mod)) % mod;
    return res;
}

uint64_t cantor_number(uint64_t vals[], uint64_t n, uint64_t mod){
    uint64_t res = 0;
    for (uint64_t i = 0; i < n; i++){
        uint64_t cnt = 0;
        for (uint64_t j = i + 1; j < n; j++)
            if (vals[j] < vals[i])
                cnt++;
        res = (res + mul(cnt, factorial(n - i - 1, mod), mod)) % mod;
    }
    return res;
}
