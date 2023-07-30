#include "number_theory.h"
#include "big_num.h"

uint64_t cpow(uint64_t x, uint64_t w2, uint64_t n, uint64_t mod) {
    uint64_t rx = 1, ry = 0;
    uint64_t y = 1, t = 0;
    for (; n; n >>= 1){
        if (n & 1){
            t = mul(rx, x, mod) + mul(ry, mul(y, w2, mod), mod);
            ry = mul(ry, x, mod) + mul(rx, y, mod);
            rx = t >= mod ? t - mod : t;
            ry = ry >= mod ? ry - mod : ry;
        }
        t = mul(x, x, mod) + mul(y, mul(y, w2, mod), mod);
        y = mul(y, x, mod); y = y + y;
        x = t >= mod ? t - mod : t;
        y = y >= mod ? y - mod : y;
    }
    assert(ry == 0 && "lagrange theorem");
    return rx;
}

uint64_t cipolla(uint64_t x, uint64_t mod) {
    if (x == 0)
        return 0;
    uint64_t modm1d2 = (mod - 1) >> 1;
    if (pow(x, modm1d2, mod) != 1)
        return -1;
    for (uint32_t i = 1;; i++){
        uint64_t w2 = i * i;
        w2 = w2 >= x ? w2 - x : w2 + mod - x;
        if (w2 == 0)
            return i;
        uint64_t y = pow(w2, modm1d2, mod);
        if (y == mod - 1)
            return cpow(i, w2, modm1d2 + 1, mod);
    }
    assert(false && "unreachable code");
    return -1;
}

BigInt cpow(BigInt x, BigInt w2, BigInt n, BigInt mod) {
    BigInt rx = 1, ry = 0;
    BigInt y = 1, t = 0;
    for (; n; n >>= 1){
        if (n & 1){
            t = mul(rx, x, mod) + mul(ry, mul(y, w2, mod), mod);
            ry = mul(ry, x, mod) + mul(rx, y, mod);
            rx = t >= mod ? t - mod : t;
            ry = ry >= mod ? ry - mod : ry;
        }
        t = mul(x, x, mod) + mul(y, mul(y, w2, mod), mod);
        y = mul(y, x, mod); y = y + y;
        x = t >= mod ? t - mod : t;
        y = y >= mod ? y - mod : y;
    }
    assert(ry == 0 && "lagrange theorem");
    return rx;
}

BigInt cipolla(BigInt x, BigInt mod) {
    if (x == 0)
        return 0;
    if (jacobi(x, mod) != 1)
        return -1;
    for (uint64_t i = 1;; i++){
        BigInt w2 = i * i;
        w2 = w2 >= x ? w2 - x : w2 + mod - x;
        if (w2 == 0)
            return i;
        int y = jacobi(w2, mod);
        if (y == -1)
            return cpow(i, w2, (mod + 1) >> 1, mod);
    }
    assert(false && "unreachable code");
    return -1;
}