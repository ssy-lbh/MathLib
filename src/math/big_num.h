#ifndef BIG_NUM_H
#define BIG_NUM_H

#include "math_base.h"

#include <cstdint>

#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>

class BigInt;
class BigFrac;
class BigFloat;

class BigInt {
public:
    mpz_t n;

    friend class BigFrac;
    friend class BigFloat;

public:
    static const BigInt zero;
    static const BigInt one;

    class Random {
    public:
        gmp_randstate_t state;

    public:
        Random();
        Random(const Random& r);
        Random(unsigned long seed);
        ~Random();

        BigInt operator()(const BigInt& n); // [0, n)
        BigInt operator()(const BigInt& a, const BigInt& b); // [a, b)
    };

    BigInt();
    BigInt(BigInt&& b);
    BigInt(const BigInt& b);
    BigInt(const char* s);
    BigInt(const char* s, int radix);
    BigInt(const int& x);
    BigInt(const unsigned int& x);
    BigInt(const long& x);
    BigInt(const unsigned long& x);
    BigInt(const long long& x);
    BigInt(const unsigned long long& x);
    explicit BigInt(const BigFrac& x);
    explicit BigInt(const BigFloat& x);
    ~BigInt();
    BigInt& operator=(BigInt&& b);
    BigInt& operator=(const BigInt& b);
    BigInt& operator=(const char* s);
    BigInt& operator=(const int& x);
    BigInt& operator=(const unsigned int& x);
    BigInt& operator=(const long& x);
    BigInt& operator=(const unsigned long& x);
    BigInt& operator=(const long long& x);
    BigInt& operator=(const unsigned long long& x);
    BigInt& operator=(const BigFrac& x);
    BigInt& operator=(const BigFloat& x);
    BigInt& operator++();
    BigInt& operator--();
    BigInt operator++(int);
    BigInt operator--(int);
    BigInt operator-() const;
    BigInt operator+() const;
    BigInt operator~() const;

    explicit operator bool() const;
    explicit operator int() const;
    explicit operator unsigned int() const;
    explicit operator long() const;
    explicit operator unsigned long() const;
    explicit operator long long() const;
    explicit operator unsigned long long() const;

    void print() const;
    void print(int radix) const;
    void print(FILE* file) const;
    void print(FILE* file, int radix) const;
    void scan();
    void scan(int radix);
    void scan(FILE* file);
    void scan(FILE* file, int radix);

    void swap(BigInt& b);
    static void swap(BigInt& a, BigInt& b);

    int sgn() const;
    BigInt abs();
    BigInt pow(unsigned long exp);
    BigInt pow(const BigInt& exp, const BigInt& mod);
    BigInt sqrt();
    BigInt sqrtrem(BigInt& rem);
    BigInt root(unsigned long exp);
    BigInt rootrem(BigInt& rem, unsigned long exp);
    BigInt gcd(const BigInt& b);
    BigInt lcm(const BigInt& b);
    BigInt factorial();
    BigInt fibonacci();
    BigInt lucas();
    BigInt binomial(const BigInt& b);
    BigInt binomial(const BigInt& b, const BigInt& m);
    BigInt multinomial(const BigInt& b);
    BigInt multinomial(const BigInt& b, const BigInt& m);
    bool is_prime() const;
    bool is_probab_prime(int reps) const;
    bool is_perfect_power() const;
    bool is_perfect_square() const;
    bool testbit(unsigned long bit) const;
    void setbit(unsigned long bit);
    void clrbit(unsigned long bit);
    void combit(unsigned long bit);
    unsigned long scanbit0(unsigned long bit);
    unsigned long scanbit1(unsigned long bit);
    void gcdext(const BigInt& b, BigInt& s, BigInt& t);
    void nextprime();
};

constexpr bool is_conjugate_identical(const BigInt&) { return true; }
constexpr bool is_commutative(const BigInt&) { return true; }
constexpr bool is_associative(const BigInt&) { return true; }
constexpr bool is_alternative(const BigInt&) { return true; }

constexpr bool is_unital(const BigInt&) { return true; }
constexpr bool is_dividable(const BigInt&) { return true; }

BigInt& operator+=(BigInt& a, const BigInt& b);
BigInt& operator-=(BigInt& a, const BigInt& b);
BigInt& operator*=(BigInt& a, const BigInt& b);
BigInt& operator/=(BigInt& a, const BigInt& b);
BigInt& operator%=(BigInt& a, const BigInt& b);
BigInt& operator<<=(BigInt& a, const unsigned long& b);
BigInt& operator>>=(BigInt& a, const unsigned long& b);
BigInt& operator&=(BigInt& a, const BigInt& b);
BigInt& operator|=(BigInt& a, const BigInt& b);
BigInt& operator^=(BigInt& a, const BigInt& b);
BigInt operator<<(const BigInt& a, const unsigned long& b);
BigInt operator>>(const BigInt& a, const unsigned long& b);
BigInt operator&(const BigInt& a, const BigInt& b);
BigInt operator|(const BigInt& a, const BigInt& b);
BigInt operator^(const BigInt& a, const BigInt& b);
BigInt operator+(const BigInt& a, const BigInt& b);
BigInt operator-(const BigInt& a, const BigInt& b);
BigInt operator*(const BigInt& a, const BigInt& b);
BigInt operator/(const BigInt& a, const BigInt& b);
BigInt operator%(const BigInt& a, const BigInt& b);
bool operator<(const BigInt& a, const BigInt& b);
bool operator>(const BigInt& a, const BigInt& b);
bool operator<=(const BigInt& a, const BigInt& b);
bool operator>=(const BigInt& a, const BigInt& b);
bool operator==(const BigInt& a, const BigInt& b);
bool operator!=(const BigInt& a, const BigInt& b);
bool operator<(const BigInt& a, long b);
bool operator>(const BigInt& a, long b);
bool operator<=(const BigInt& a, long b);
bool operator>=(const BigInt& a, long b);
bool operator==(const BigInt& a, long b);
bool operator!=(const BigInt& a, long b);

inline BigInt ident(BigInt) { return BigInt::one; }
inline BigInt zero(BigInt) { return BigInt::zero; }
template <typename U> std::enable_if_t<std::is_arithmetic_v<U>, BigInt> num(BigInt, U n) { return BigInt(n); }
inline BigInt conj(BigInt x) { return x; }
inline BigInt norm(BigInt x) { return x.abs(); }
inline BigInt norm2(BigInt x) { return x * x; }
inline int line(BigInt) { return 1; }

inline void print(const BigInt& x, int l) { x.print(); }
inline void print(const BigInt& x, int radix, int l) { x.print(radix); }

int sgn(const BigInt& x);
BigInt abs(const BigInt& x);
BigInt mul(const BigInt& x, const BigInt& y, const BigInt& mod);
BigInt pow(const BigInt& x, unsigned long exp);
BigInt pow(const BigInt& x, const BigInt& exp, const BigInt& mod);
BigInt sqrt(const BigInt& x);
BigInt sqrtrem(const BigInt& x, BigInt& rem);
BigInt root(const BigInt& x, unsigned long exp);
BigInt rootrem(const BigInt& x, BigInt& rem, unsigned long exp);
BigInt gcd(const BigInt& x, const BigInt& y);
BigInt lcm(const BigInt& x, const BigInt& y);
BigInt exgcd(const BigInt& x, const BigInt& y, BigInt& a, BigInt& b);
BigInt factorial(const BigInt& x);
BigInt fibonacci(const BigInt& x);
BigInt lucas(const BigInt& x);
BigInt binomial(const BigInt& x, const BigInt& y);
BigInt binomial(const BigInt& x, const BigInt& y, const BigInt& m);
BigInt multinomial(const BigInt& x, const BigInt& y);
BigInt multinomial(const BigInt& x, const BigInt& y, const BigInt& m);
int legendre(const BigInt& x, const BigInt& p);
int jacobi(const BigInt& x, const BigInt& p);
int kronecker(const BigInt& x, const BigInt& y);
BigInt nextprime(const BigInt& x);
BigInt randmod(const BigInt& n);
BigInt rand(const BigInt& a, const BigInt& b);
BigInt randbits(unsigned long bits);
uint64_t size(const BigInt& x);
uint64_t sizeinbase(const BigInt& x, int base);
BigInt fdivq(const BigInt& n, const BigInt& d);
BigInt fdivr(const BigInt& n, const BigInt& d);
void fdivqr(const BigInt& n, const BigInt& d, BigInt& q, BigInt& r);
BigInt cdivq(const BigInt& n, const BigInt& d);
BigInt cdivr(const BigInt& n, const BigInt& d);
void cdivqr(const BigInt& n, const BigInt& d, BigInt& q, BigInt& r);
BigInt tdivq(const BigInt& n, const BigInt& d);
BigInt tdivr(const BigInt& n, const BigInt& d);
void tdivqr(const BigInt& n, const BigInt& d, BigInt& q, BigInt& r);

extern BigInt::Random default_bigint_random;

class BigFrac {
public:
    mpq_t n;

    friend class BigInt;
    friend class BigFloat;

public:
    static const BigFrac zero;
    static const BigFrac one;

    BigFrac();
    BigFrac(const BigFrac& b);
    BigFrac(const BigInt& x);
    BigFrac(const BigInt& x, const BigInt& y);
    BigFrac(const char* s);
    BigFrac(const char* s, int radix);
    BigFrac(const int& x);
    BigFrac(const unsigned int& x);
    BigFrac(const long& x);
    BigFrac(const unsigned long& x);
    BigFrac(const long long& x);
    BigFrac(const unsigned long long& x);
    ~BigFrac();
    BigFrac& operator=(const BigFrac& b);
    BigFrac& operator=(const char* s);
    BigFrac& operator=(const int& x);
    BigFrac& operator=(const unsigned int& x);
    BigFrac& operator=(const long& x);
    BigFrac& operator=(const unsigned long& x);
    BigFrac& operator=(const long long& x);
    BigFrac& operator=(const unsigned long long& x);
    BigFrac& operator++();
    BigFrac& operator--();
    BigFrac operator++(int);
    BigFrac operator--(int);
    BigFrac operator-() const;
    BigFrac operator+() const;

    explicit operator bool() const;
    explicit operator float() const;
    explicit operator double() const;
    explicit operator long double() const;

    void print() const;
    void print(int radix) const;
    void print(FILE* file) const;
    void print(FILE* file, int radix) const;
    void scan();
    void scan(int radix);
    void scan(FILE* file);
    void scan(FILE* file, int radix);

    void get_num(BigInt& num);
    void get_den(BigInt& den);

    void swap(BigFrac& b);
    static void swap(BigFrac& a, BigFrac& b);

    int sgn() const;
    BigFrac abs();
    BigFrac inv();
    BigFrac pow(unsigned long exp);
    BigFrac pow(const BigInt& exp, const BigInt& mod);

    BigFrac& simplify();
};

constexpr bool is_conjugate_identical(const BigFrac&) { return true; }
constexpr bool is_commutative(const BigFrac&) { return true; }
constexpr bool is_associative(const BigFrac&) { return true; }
constexpr bool is_alternative(const BigFrac&) { return true; }

constexpr bool is_unital(const BigFrac&) { return true; }
constexpr bool is_dividable(const BigFrac&) { return true; }

BigFrac& operator+=(BigFrac& a, const BigFrac& b);
BigFrac& operator-=(BigFrac& a, const BigFrac& b);
BigFrac& operator*=(BigFrac& a, const BigFrac& b);
BigFrac& operator/=(BigFrac& a, const BigFrac& b);
BigFrac& operator<<=(BigFrac& a, const unsigned long& b);
BigFrac& operator>>=(BigFrac& a, const unsigned long& b);
BigFrac operator+(const BigFrac& a, const BigFrac& b);
BigFrac operator-(const BigFrac& a, const BigFrac& b);
BigFrac operator*(const BigFrac& a, const BigFrac& b);
BigFrac operator/(const BigFrac& a, const BigFrac& b);
BigFrac operator<<(const BigFrac& a, const unsigned long& b);
BigFrac operator>>(const BigFrac& a, const unsigned long& b);
bool operator<(const BigFrac& a, const BigFrac& b);
bool operator>(const BigFrac& a, const BigFrac& b);
bool operator<=(const BigFrac& a, const BigFrac& b);
bool operator>=(const BigFrac& a, const BigFrac& b);
bool operator==(const BigFrac& a, const BigFrac& b);
bool operator!=(const BigFrac& a, const BigFrac& b);

inline BigFrac ident(BigFrac) { return BigFrac::one; }
inline BigFrac zero(BigFrac) { return BigFrac::zero; }
template <typename U> std::enable_if_t<std::is_arithmetic_v<U>, BigFrac> num(BigFrac, U n) { return BigFrac(n); }
inline BigFrac conj(BigFrac x) { return x; }
inline BigFrac inv(BigFrac x) { return x.inv(); }
inline BigFrac norm(BigFrac x) { BigFrac t = x; t.abs(); return t; }
inline BigFrac norm2(BigFrac x) { return x * x; }
inline int line(BigFrac) { return 1; }

inline void print(const BigFrac& x, int l) { x.print(); }
inline void print(const BigFrac& x, int radix, int l) { x.print(radix); }

int sgn(const BigFrac& x);
BigFrac abs(const BigFrac& x);
BigFrac inv(const BigFrac& x);
BigFrac pow(const BigFrac& x, unsigned long exp);
BigFrac pow(const BigFrac& x, const BigInt& exp, const BigInt& mod);

class BigFloat {
public:
    mpfr_t n;

    friend class BigInt;
    friend class BigFrac;

public:
    static const BigFloat zero;
    static const BigFloat one;

    static const BigFloat pi;
    static const BigFloat e;

    BigFloat();
    BigFloat(const BigFloat& b);
    BigFloat(const BigInt& x);
    BigFloat(const BigFrac& x);
    BigFloat(const char* s);
    BigFloat(const char* s, int radix);
    BigFloat(const int& x);
    BigFloat(const unsigned int& x);
    BigFloat(const long& x);
    BigFloat(const unsigned long& x);
    BigFloat(const long long& x);
    BigFloat(const unsigned long long& x);
    BigFloat(const float& x);
    BigFloat(const double& x);
    BigFloat(const long double& x);
    ~BigFloat();
    BigFloat& operator=(const BigFloat& b);
    BigFloat& operator=(const char* s);
    BigFloat& operator=(const int& x);
    BigFloat& operator=(const unsigned int& x);
    BigFloat& operator=(const long& x);
    BigFloat& operator=(const unsigned long& x);
    BigFloat& operator=(const long long& x);
    BigFloat& operator=(const unsigned long long& x);
    BigFloat& operator=(const float& x);
    BigFloat& operator=(const double& x);
    BigFloat& operator=(const long double& x);
    BigFloat& operator=(const mpz_t& x);
    BigFloat& operator=(const mpq_t& x);
    BigFloat& operator=(const mpf_t& x);
    BigFloat& operator=(const BigInt& x);
    BigFloat& operator=(const BigFrac& x);
    BigFloat& operator++();
    BigFloat& operator--();
    BigFloat operator++(int);
    BigFloat operator--(int);
    BigFloat operator-() const;
    BigFloat operator+() const;

    explicit operator bool() const;
    explicit operator float() const;
    explicit operator double() const;
    explicit operator long double() const;

    void print() const;
    void print(int radix) const;
    void print(FILE* file) const;
    void print(FILE* file, int radix) const;
    void scan();
    void scan(int radix);
    void scan(FILE* file);
    void scan(FILE* file, int radix);

    void swap(BigFloat& b);
    static void swap(BigFloat& a, BigFloat& b);

    int sgn() const;
    BigFloat abs();
    BigFloat pow(unsigned long exp);
    BigFloat sqrt();
    BigInt floor();
    BigInt ceil();
    BigInt trunc();
};

constexpr bool is_conjugate_identical(const BigFloat&) { return true; }
constexpr bool is_commutative(const BigFloat&) { return true; }
constexpr bool is_associative(const BigFloat&) { return true; }
constexpr bool is_alternative(const BigFloat&) { return true; }

constexpr bool is_unital(const BigFloat&) { return true; }
constexpr bool is_dividable(const BigFloat&) { return true; }

BigFloat& operator+=(BigFloat& a, const BigFloat& b);
BigFloat& operator-=(BigFloat& a, const BigFloat& b);
BigFloat& operator*=(BigFloat& a, const BigFloat& b);
BigFloat& operator/=(BigFloat& a, const BigFloat& b);
BigFloat& operator%=(BigFloat& a, const BigFloat& b);
BigFloat& operator<<=(BigFloat& a, const unsigned long& b);
BigFloat& operator>>=(BigFloat& a, const unsigned long& b);
BigFloat operator+(const BigFloat& a, const BigFloat& b);
BigFloat operator-(const BigFloat& a, const BigFloat& b);
BigFloat operator*(const BigFloat& a, const BigFloat& b);
BigFloat operator/(const BigFloat& a, const BigFloat& b);
BigFloat operator%(const BigFloat& a, const BigFloat& b);
BigFloat operator<<(const BigFloat& a, const unsigned long& b);
BigFloat operator>>(const BigFloat& a, const unsigned long& b);
bool operator<(const BigFloat& a, const BigFloat& b);
bool operator>(const BigFloat& a, const BigFloat& b);
bool operator<=(const BigFloat& a, const BigFloat& b);
bool operator>=(const BigFloat& a, const BigFloat& b);
bool operator==(const BigFloat& a, const BigFloat& b);
bool operator!=(const BigFloat& a, const BigFloat& b);

inline BigFloat ident(BigFloat) { return BigFloat::one; }
inline BigFloat zero(BigFloat) { return BigFloat::zero; }
template <typename U> std::enable_if_t<std::is_arithmetic_v<U>, BigFloat> num(BigFloat, U n) { return BigFloat(n); }
inline BigFloat conj(const BigFloat& x) { return x; }
inline BigFloat inv(const BigFloat& x) { return BigFloat::one / x; }
inline BigFloat norm(const BigFloat& x) { BigFloat t = x; t.abs(); return t; }
inline BigFloat norm2(const BigFloat& x) { return x * x; }
inline int line(BigFloat) { return 1; }

inline void print(const BigFloat& x, int l) { x.print(); }
inline void print(const BigFloat& x, int radix, int l) { x.print(radix); }

int sgn(const BigInt& x);
BigFloat abs(const BigFloat& x);
BigFloat pow(const BigFloat& x, unsigned long exp);
BigFloat pow(const BigFloat& x, const BigFloat& exp);
BigFloat sqrt(const BigFloat& x);
BigFloat cbrt(const BigFloat& x);
BigInt floor(const BigFloat& x);
BigInt ceil(const BigFloat& x);
BigInt trunc(const BigFloat& x);

BigFloat exp(const BigFloat& x);
BigFloat exp2(const BigFloat& x);
BigFloat exp10(const BigFloat& x);
BigFloat expm1(const BigFloat& x);
BigFloat log(const BigFloat& x);
BigFloat log2(const BigFloat& x);
BigFloat log10(const BigFloat& x);
BigFloat log1p(const BigFloat& x);
BigFloat sin(const BigFloat& x);
BigFloat cos(const BigFloat& x);
BigFloat tan(const BigFloat& x);
BigFloat sinh(const BigFloat& x);
BigFloat cosh(const BigFloat& x);
BigFloat tanh(const BigFloat& x);
BigFloat asin(const BigFloat& x);
BigFloat acos(const BigFloat& x);
BigFloat atan(const BigFloat& x);
BigFloat asinh(const BigFloat& x);
BigFloat acosh(const BigFloat& x);
BigFloat atanh(const BigFloat& x);
BigFloat atan2(const BigFloat& y, const BigFloat& x);
BigFloat erf(const BigFloat& x);
BigFloat gamma(const BigFloat& x);
BigFloat zeta(const BigFloat& x);

inline BigInt inv(const BigInt& x, const BigInt& mod){
    BigInt y = 0, z = 0;
    if (exgcd(x, mod, y, z) != 1)
        return -1;
    return y < 0 ? y + mod : y;
}

template <typename T>
bool miller_prime_check(const T& n, const T& a) {
    T d = n - 1;
    while ((d & 1) == 0) d >>= 1;
    T t = pow(a, d, n);
    while (d != n - 1 && t != 1 && t != n - 1) {
        t = mul(t, t, n);
        d <<= 1;
    }
    return t == n - 1 || (d & 1) == 1;
}

// Gary Lee Miller's primality proof
// Based on Riemann Hypothesis
template <typename T>
bool miller_prime_proof(const T& n) {
    if (n < 3 || (n & 1) == 0)
        return n == 2;
    auto t = log(n);
    T reps = floor(t * log2(t)) << 1;
    for (T i = 2; i <= reps; i++)
        if (!miller_prime_check(n, i))
            return false;
    return true;
}

bool miller_rabin(BigInt n);
bool miller_rabin(BigInt n, int reps);
BigInt pollard_rho(BigInt n);
uint64_t factorize(BigInt n, BigInt prime[], uint64_t exp[], uint64_t len);
uint64_t factorize(BigInt n, BigInt prime[], uint64_t exp[], uint64_t len, uint64_t filter);

BigInt cpow(BigInt x, BigInt w2, BigInt n, BigInt mod);
BigInt cipolla(BigInt x, BigInt mod);

BigInt ecm_factorize(BigInt n);

#endif