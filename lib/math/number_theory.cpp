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

        if ((a & 1) == 0)
            if (n % 8 == 3 || n % 8 == 5)
                symbol = -symbol;
            a >>= 1;
            continue;

        if (a % 4 == 3 && n % 4 == 3)
            symbol = -symbol;

        uint64_t t = a;
        a = n;
        n = t;
    }
}
