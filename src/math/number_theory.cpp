#include "number_theory.h"

#include <random>

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

std::mt19937_64 rng(std::random_device{}());

uint64_t randmod(uint64_t mod){
    std::uniform_int_distribution<uint64_t> dist(0, mod - 1);
    return dist(rng);
}

uint64_t rand(uint64_t l, uint64_t h){
    std::uniform_int_distribution<uint64_t> dist(l, h);
    return dist(rng);
}

uint32_t randmod(uint32_t mod){
    std::uniform_int_distribution<uint32_t> dist(0, mod - 1);
    return dist(rng);
}

uint32_t rand(uint32_t l, uint32_t h){
    std::uniform_int_distribution<uint32_t> dist(l, h);
    return dist(rng);
}

void init_random(){
    rng.seed(std::random_device{}());
}

void init_random(uint64_t seed){
    rng.seed(seed);
}

uint64_t pi_limit(uint64_t x){
    if (x <= 1)
        return 0;
    double log_x = log(x);
    double log_log_x = log(log_x);
    return (uint64_t)ceil(x / log_x * (1 + 1.2762 / log_x + 1.2762 * log_log_x / log_x / log_x));
}

int legendre(uint64_t a, uint64_t p){
    uint64_t res = pow(a, (p - 1) >> 1, p);
    return res == p - 1 ? -1 : res;
}

int jacobi(uint64_t a, uint64_t n){
    int symbol = 1;
    if (n < 1 || (n & 1) == 0)
        return 0;
    while(true){
        if (n == 1)
            return symbol;
        a %= n;
        if (a == 0)
            return 0;
        if (a == 1)
            return symbol;

        if ((a & 1) == 0){
            if (n % 8 == 3 || n % 8 == 5)
                symbol = -symbol;
            a >>= 1;
            continue;
        }

        if (a % 4 == 3 && n % 4 == 3)
            symbol = -symbol;

        uint64_t t = a;
        a = n;
        n = t;
    }
}

// 狄利克雷卷积
void dirichlet_convolution(uint64_t a[], uint64_t b[], uint64_t c[], uint64_t n){
    for (uint64_t i = 1; i < n; i++)
        c[i] = 0;
    for (uint64_t d = 1; d < n; d++)
        for (uint64_t i = d, j = 1; i < n; i += d, j++)
            c[i] += a[d] * b[j];
}
