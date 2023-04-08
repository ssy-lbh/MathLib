#include "number_theory.h"

// 返回已得到质数的乘积
static uint64_t totient_step(uint64_t x, uint64_t& n){
    // 进行分解
    uint64_t r = pollard_rho(x);
    // x 为质数, 处理质数并返回
    if (r == x){
        n = n / r * (r - 1);
        return r;
    }
    // x 不是质数, r1可能为合数, 递归处理
    uint64_t s = totient_step(r, n);
    uint64_t g = gcd(x, s);
    while (g > 1) {
        x /= g;
        g = gcd(x, s);
    }
    if (x == 1)
        return s;
    uint64_t t = totient_step(x, n);
    return s * t;
}

uint64_t totient(uint64_t n){
    if (n == 1)
        return 1;
    totient_step(n, n);
    return n;
}
