#include "math_base.h"
#include "big_num.h"

#include <vector>

BigFrac sqrt_approx1(BigInt x, uint32_t iter = 10){
    BigFrac a = sqrt(x);
    for (uint32_t i = 0; i < iter; i++){
        BigFrac a2 = a * a;
        a = a * (a2 + 3 * x) / (3 * a2 + x);
    }
    return a;
}

BigFrac sqrt_approx2(BigInt x, uint32_t iter = 10){
    BigInt a = sqrt(x);
    // (y - a) mod (y^2 - x)
    BigInt f0 = -a, f1 = 1;
    for (uint32_t i = 0; i < iter; i++){
        BigInt nf0 = f0 * f0 + x * f1 * f1;
        BigInt nf1 = 2 * f0 * f1;
        f0 = nf0;
        f1 = nf1;
    }
    return BigFrac(f0, f1);
}

// q * q - t * p * p = 1
// @return p / q
BigFrac pell_equation1(BigInt n){
    BigInt m = sqrt(n);
    std::vector<BigInt> f;
    f.push_back(m);
    BigInt a = m, b = 1;
    while (true){
        b = (n - a * a) / b;
        if (b == 1) break;
        BigInt c = (m + a) / b;
        a = c * b - a;
        f.push_back(c);
    }
    BigInt p = 0, q = 1;
    for (int i = f.size() - 1; i >= 0; --i){
        BigInt t = p;
        p = q;
        q = f[i] * q + t;
    }
    return BigFrac(p, q);
}