#ifndef NUMBER_THEORY_H
#define NUMBER_THEORY_H

#include "math_base.h"

#include <cstdlib>
#include <cassert>
#include <cstdint>

constexpr uint64_t pow(uint64_t x, uint32_t y){
    uint64_t r = 1;
    while (y){
        if (y & 1)
            r *= x;
        x *= x;
        y >>= 1;
    }
    return r;
}

constexpr uint64_t pow(uint64_t x, uint32_t y, uint32_t mod){
    uint64_t r = 1;
    while (y){
        if (y & 1)
            r = (r * x) % mod;
        x = (x * x) % mod;
        y >>= 1;
    }
    return r;
}

constexpr uint64_t pow(uint64_t x, uint64_t y, uint64_t mod){
    uint64_t r = 1;
    while (y){
        if (y & 1)
            r = mul(r, x, mod);
        x = mul(x, x, mod);
        y >>= 1;
    }
    return r;
}

constexpr uint64_t fact(uint32_t x){
    uint64_t r = 1;
    for (uint32_t i = 2; i <= x; i++)
        r *= i;
    return r;
}

constexpr uint32_t gcd(uint32_t x, uint32_t y){
    return (x == 0) ? y : gcd(y % x, x);
}

constexpr uint64_t gcd(uint64_t x, uint64_t y){
    return (x == 0) ? y : gcd(y % x, x);
}

constexpr uint32_t lcm(uint32_t x, uint32_t y){
    return x / gcd(x, y) * y;
}

constexpr uint64_t lcm(uint64_t x, uint64_t y){
    return x / gcd(x, y) * y;
}

// 扩展欧几里得, ax + by = gcd(a, b)
constexpr uint64_t exgcd(uint64_t a, uint64_t b, int64_t& x, int64_t& y){
    if (b == 0){
        x = 1;
        y = 0;
        return a;
    }
    uint64_t r = exgcd(b, a % b, y, x);
    y -= a / b * x;
    return r;
}

constexpr uint64_t inv(uint64_t x, uint64_t mod){
    int64_t y = 0, z = 0;
    exgcd(x, mod, y, z);
    return (y % mod + mod) % mod;
}

inline uint64_t rand(uint64_t mod){
    return (((uint64_t)rand() << 60) | ((uint64_t)rand() << 45) | ((uint64_t)rand() << 30) | (rand() << 15) | rand()) % mod;
}

inline uint64_t rand(uint64_t l, uint64_t h){
    return l + rand(h - l + 1);
}

inline uint32_t rand(uint32_t mod){
    return (((uint32_t)rand() << 30) | (rand() << 15) | rand()) % mod;
}

inline uint32_t rand(uint32_t l, uint32_t h){
    return l + rand(h - l + 1);
}

constexpr uint64_t mul(uint64_t a, uint64_t b, uint64_t mod){
    uint64_t c = (long double)a / mod * b;
    uint64_t res = (uint64_t)a * b - (uint64_t)c * mod;
    return (res + mod) % mod;
}

constexpr uint64_t pow(uint64_t x, uint64_t y, uint64_t mod){
    uint64_t r = 1;
    while (y){
        if (y & 1)
            r = mul(r, x, mod);
        x = mul(x, x, mod);
        y >>= 1;
    }
    return r;
}

uint64_t cpow(uint64_t x, uint64_t w2, uint64_t n, uint64_t mod);
uint64_t cipolla(uint64_t x, uint64_t mod);

// @return is_prime?
bool miller_rabin(uint32_t n);
bool miller_rabin(uint64_t n);

uint64_t pollard_rho(uint64_t n);
uint32_t pollard_rho(uint64_t n, uint64_t factor[], uint32_t exp[], uint32_t len);

uint64_t totient(uint64_t n);

bool check_root(uint64_t x, uint64_t p, uint64_t pm1_prime[], uint32_t cnt);
uint64_t find_root(uint64_t p, uint64_t pm1_prime[], uint32_t cnt);
bool has_root(uint64_t pm1_prime[], uint32_t exp[], uint32_t cnt);

// a^x === b (mod p)
uint64_t pohlig_hellman_log(uint64_t a, uint64_t b, uint64_t p);
// g^x === b (mod p)
uint64_t pohlig_hellman(uint64_t g, uint64_t b, uint64_t p);
uint64_t pohlig_hellman(uint64_t g, uint64_t b, uint64_t p, uint64_t pm1_prime[], uint32_t exp[], uint32_t cnt);

// 有限循环群, 可以获取生成元, 目前的种类有
// 1. Zn, 模N同余群, 有限循环群
// 2. 复数原根

// 可以后续添加, 模板位置不限
//TODO R, K自动计算
template <int N> struct TModRoot;
template <> struct TModRoot<3> { static constexpr uint32_t R = 1, K = 1, G = 2; };
template <> struct TModRoot<5> { static constexpr uint32_t R = 1, K = 2, G = 2; };
template <> struct TModRoot<17> { static constexpr uint32_t R = 1, K = 4, G = 3; };
template <> struct TModRoot<97> { static constexpr uint32_t R = 3, K = 5, G = 5; };
template <> struct TModRoot<193> { static constexpr uint32_t R = 3, K = 6, G = 5; };
template <> struct TModRoot<257> { static constexpr uint32_t R = 1, K = 8, G = 3; };
template <> struct TModRoot<7681> { static constexpr uint32_t R = 15, K = 9, G = 17; };
template <> struct TModRoot<12289> { static constexpr uint32_t R = 3, K = 12, G = 11; };
template <> struct TModRoot<40961> { static constexpr uint32_t R = 5, K = 13, G = 3; };
template <> struct TModRoot<65537> { static constexpr uint32_t R = 1, K = 16, G = 3; };
template <> struct TModRoot<786433> { static constexpr uint32_t R = 3, K = 18, G = 10; };
template <> struct TModRoot<5767169> { static constexpr uint32_t R = 11, K = 19, G = 3; };
template <> struct TModRoot<7340033> { static constexpr uint32_t R = 7, K = 20, G = 3; };
template <> struct TModRoot<23068673> { static constexpr uint32_t R = 11, K = 21, G = 3; };
template <> struct TModRoot<104857601> { static constexpr uint32_t R = 25, K = 22, G = 3; };
template <> struct TModRoot<167772161> { static constexpr uint32_t R = 5, K = 25, G = 3; };
template <> struct TModRoot<469762049> { static constexpr uint32_t R = 7, K = 26, G = 3; };
template <> struct TModRoot<998244353> { static constexpr uint32_t R = 119, K = 23, G = 3; };
template <> struct TModRoot<1004535809> { static constexpr uint32_t R = 479, K = 21, G = 3; };

// N为大小, G为原根(生成元)
template <int N>
class TMod {
public:
    uint32_t n;

    constexpr TMod() : n(0) {}
    constexpr TMod(uint32_t n) : n(n % N) {}
    constexpr TMod(const TMod<N>&) = default;
    constexpr TMod<N>& operator=(const TMod<N>&) = default;
    constexpr TMod<N>& operator=(uint32_t n) { this->n = n; return *this; }
};

template <int N> constexpr bool is_commutative(TMod<N>) { return true; }
template <int N> constexpr bool is_associative(TMod<N>) { return true; }
template <int N> constexpr bool is_alternative(TMod<N>) { return true; }

template <int N> constexpr bool is_unital(TMod<N>) { return true; }
template <int N> constexpr bool is_dividable(TMod<N>) { return true; }

// congruence group, 模N同余群
template <int N> constexpr TMod<N> operator+(TMod<N> x, TMod<N> y) { uint32_t s = x.n + y.n; return TMod<N>(s >= N ? s - N : s); }
template <int N> constexpr TMod<N> operator+(TMod<N> x, uint32_t y) { uint32_t s = x.n + y; return TMod<N>(s >= N ? s - N : s); }
template <int N> constexpr TMod<N> operator-(TMod<N> x, TMod<N> y) { return TMod<N>(x.n < y.n ? x.n + N - y.n : x.n - y.n); }
template <int N> constexpr TMod<N> operator-(TMod<N> x, uint32_t y) { return TMod<N>(x.n < y ? x.n + N - y : x.n - y); }
template <int N> constexpr TMod<N> operator*(TMod<N> x, TMod<N> y) { return TMod<N>((((uint64_t)x.n) * y.n) % N); }
template <int N> constexpr TMod<N> operator*(TMod<N> x, uint32_t y) { return TMod<N>((((uint64_t)x.n) * y) % N); }
template <int N> constexpr TMod<N> operator/(TMod<N> x, TMod<N> y) { return inv(y) * x; }
template <int N> constexpr TMod<N> operator/(TMod<N> x, uint32_t y) { return inv(TMod<N>(y)) * x; }
template <int N> constexpr TMod<N> operator%(TMod<N> x, uint32_t y) { return TMod<N>(x.n % y); }
template <int N> constexpr TMod<N> operator+=(TMod<N>& x, TMod<N> y) { return x = x + y; }
template <int N> constexpr TMod<N> operator-=(TMod<N>& x, TMod<N> y) { return x = x - y; }
template <int N> constexpr TMod<N> operator*=(TMod<N>& x, TMod<N> y) { return x = y * x; }
template <int N> constexpr TMod<N> operator/=(TMod<N>& x, TMod<N> y) { return x = x / y; }
template <int N> constexpr TMod<N> operator%=(TMod<N>& x, uint32_t y) { return x = x % y; }
template <int N> constexpr TMod<N> operator+(TMod<N> x) { return x; }
template <int N> constexpr TMod<N> operator-(TMod<N> x) { return TMod<N>(N - x.n); }

template <int N> constexpr bool operator==(TMod<N> x, TMod<N> y) { return x.n == y.n; }
template <int N> constexpr bool operator!=(TMod<N> x, TMod<N> y) { return x.n != y.n; }

template <int N> constexpr bool operator==(TMod<N> x, uint32_t y) { return x.n == y; }
template <int N> constexpr bool operator!=(TMod<N> x, uint32_t y) { return x.n != y; }
template <int N> constexpr bool operator<(TMod<N> x, uint32_t y) { return x.n < y; }
template <int N> constexpr bool operator<=(TMod<N> x, uint32_t y) { return x.n <= y; }
template <int N> constexpr bool operator>(TMod<N> x, uint32_t y) { return x.n > y; }
template <int N> constexpr bool operator>=(TMod<N> x, uint32_t y) { return x.n >= y; }

template <int N> constexpr TMod<N> ident(TMod<N>) { return TMod<N>(1); }
template <int N> constexpr TMod<N> zero(TMod<N>) { return TMod<N>(0); }
template <int N2, int N> constexpr TMod<N> gen(TMod<N>) { static_assert((N - 1) % N2 == 0, "root number not compatible"); return pow(TModRoot<N>::G, (N - 1) / N2, N); }
template <int N> constexpr TMod<N> gen(TMod<N>, uint32_t n) { assert((N - 1) % n == 0); return pow(TModRoot<N>::G, (N - 1) / n, N); }
template <int N> constexpr TMod<N> conj(TMod<N> x) { return x; }
template <int N> constexpr TMod<N> inv(TMod<N> x) { return pow(x.n, (uint32_t)N - 2, (uint32_t)N); }
template <int N> constexpr TMod<N> norm(TMod<N> x) { return x; }
template <int N> constexpr TMod<N> norm2(TMod<N> x) { return x * x; }
template <int N> constexpr int line(TMod<N>) { return 1; }

template <int N> constexpr bool isnan(TMod<N> x) { return x.n == -1; }
template <int N> constexpr TMod<N> nan(TMod<N>) { TMod<N> x; x.n = -1; return x; }
template <int N> constexpr TMod<N> pow(TMod<N> x, uint32_t n) { return pow(x.n, n % (uint32_t)(N - 1), (uint32_t)N); }

template <int N> constexpr TMod<N> sqrt(const TMod<N>& x) {
    TMod<N> r;
    r.n = cipolla(x.n, N);
    return r;
}

template <int N> void print(const TMod<N>& x, int) { printf("%d", x.n); }
template <int N> void print(const TMod<N>& x) { printf("%d\n", x.n); }

#endif /* NUMBER_THEORY_H */
