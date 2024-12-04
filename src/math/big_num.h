#ifndef BIG_NUM_H
#define BIG_NUM_H

#include "math_base.h"

#include <cstdint>
#include <string>

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
    BigInt abs() const;
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

inline BigInt::BigInt() { mpz_init(n); }
inline BigInt::BigInt(BigInt&& b) { mpz_init(n); mpz_swap(n, b.n); }
inline BigInt::BigInt(const BigInt& b) { mpz_init_set(n, b.n); }
inline BigInt::BigInt(const char* s) { mpz_init_set_str(n, s, 10); }
inline BigInt::BigInt(const char* s, int radix) { mpz_init_set_str(n, s, radix); }
inline BigInt::BigInt(const int& x) { mpz_init_set_si(n, x); }
inline BigInt::BigInt(const unsigned int& x) { mpz_init_set_ui(n, x); }
inline BigInt::BigInt(const long& x) { mpz_init_set_si(n, x); }
inline BigInt::BigInt(const unsigned long& x) { mpz_init_set_ui(n, x); }
inline BigInt::BigInt(const long long& x) { mpz_init(n); mpz_import(n, 1, -1, sizeof(long long), 0, 0, &x); }
inline BigInt::BigInt(const unsigned long long& x) { mpz_init(n); mpz_import(n, 1, -1, sizeof(unsigned long long), 0, 0, &x); }
inline BigInt::~BigInt() { mpz_clear(n); }

inline BigInt& BigInt::operator=(BigInt&& b) { mpz_swap(n, b.n); return *this; }
inline BigInt& BigInt::operator=(const BigInt& b) { mpz_set(n, b.n); return *this; }
inline BigInt& BigInt::operator=(const char* s) { mpz_set_str(n, s, 10); return *this; }
inline BigInt& BigInt::operator=(const int& x) { mpz_set_si(n, x); return *this; }
inline BigInt& BigInt::operator=(const unsigned int& x) { mpz_set_ui(n, x); return *this; }
inline BigInt& BigInt::operator=(const long& x) { mpz_set_si(n, x); return *this; }
inline BigInt& BigInt::operator=(const unsigned long& x) { mpz_set_ui(n, x); return *this; }
inline BigInt& BigInt::operator=(const long long& x) { mpz_import(n, 1, -1, sizeof(long long), 0, 0, &x); return *this; }
inline BigInt& BigInt::operator=(const unsigned long long& x) { mpz_import(n, 1, -1, sizeof(unsigned long long), 0, 0, &x); return *this; }
inline BigInt& BigInt::operator++() { mpz_add_ui(n, n, 1); return *this; }
inline BigInt& BigInt::operator--() { mpz_sub_ui(n, n, 1); return *this; }
inline BigInt BigInt::operator++(int) { BigInt t = *this; mpz_add_ui(n, n, 1); return t; }
inline BigInt BigInt::operator--(int) { BigInt t = *this; mpz_sub_ui(n, n, 1); return t; }
inline BigInt BigInt::operator-() const { BigInt t = *this; mpz_neg(t.n, t.n); return t; }
inline BigInt BigInt::operator+() const { return *this; }
inline BigInt BigInt::operator~() const { BigInt t; mpz_com(t.n, n); return t; }

inline BigInt::operator bool() const { return mpz_cmp_ui(n, 0) != 0; }
inline BigInt::operator int() const { return mpz_get_si(n); }
inline BigInt::operator unsigned int() const { return mpz_get_ui(n); }
inline BigInt::operator long() const { return mpz_get_si(n); }
inline BigInt::operator unsigned long() const { return mpz_get_ui(n); }
inline BigInt::operator long long() const { long long x; mpz_export(&x, nullptr, -1, sizeof(long long), 0, 0, n); return x; }
inline BigInt::operator unsigned long long() const { unsigned long long x; mpz_export(&x, nullptr, -1, sizeof(unsigned long long), 0, 0, n); return x; }

inline BigInt& operator+=(BigInt& a, const BigInt& b) { mpz_add(a.n, a.n, b.n); return a; }
inline BigInt& operator-=(BigInt& a, const BigInt& b) { mpz_sub(a.n, a.n, b.n); return a; }
inline BigInt& operator*=(BigInt& a, const BigInt& b) { mpz_mul(a.n, a.n, b.n); return a; }
inline BigInt& operator/=(BigInt& a, const BigInt& b) { mpz_tdiv_q(a.n, a.n, b.n); return a; }
inline BigInt& operator%=(BigInt& a, const BigInt& b) { mpz_tdiv_r(a.n, a.n, b.n); return a; }
inline BigInt& operator<<=(BigInt& a, const unsigned long& b) { mpz_mul_2exp(a.n, a.n, b); return a; }
inline BigInt& operator>>=(BigInt& a, const unsigned long& b) { mpz_tdiv_q_2exp(a.n, a.n, b); return a; }
inline BigInt& operator&=(BigInt& a, const BigInt& b) { mpz_and(a.n, a.n, b.n); return a; }
inline BigInt& operator|=(BigInt& a, const BigInt& b) { mpz_ior(a.n, a.n, b.n); return a; }
inline BigInt& operator^=(BigInt& a, const BigInt& b) { mpz_xor(a.n, a.n, b.n); return a; }
inline BigInt operator<<(const BigInt& a, const unsigned long& b) { BigInt t; mpz_mul_2exp(t.n, a.n, b); return t; }
inline BigInt operator>>(const BigInt& a, const unsigned long& b) { BigInt t; mpz_div_2exp(t.n, a.n, b); return t; }
inline BigInt operator&(const BigInt& a, const BigInt& b) { BigInt t; mpz_and(t.n, a.n, b.n); return t; }
inline BigInt operator|(const BigInt& a, const BigInt& b) { BigInt t; mpz_ior(t.n, a.n, b.n); return t; }
inline BigInt operator^(const BigInt& a, const BigInt& b) { BigInt t; mpz_xor(t.n, a.n, b.n); return t; }
inline BigInt operator+(const BigInt& a, const BigInt& b) { BigInt t; mpz_add(t.n, a.n, b.n); return t; }
inline BigInt operator-(const BigInt& a, const BigInt& b) { BigInt t; mpz_sub(t.n, a.n, b.n); return t; }
inline BigInt operator*(const BigInt& a, const BigInt& b) { BigInt t; mpz_mul(t.n, a.n, b.n); return t; }
inline BigInt operator/(const BigInt& a, const BigInt& b) { BigInt t; mpz_fdiv_q(t.n, a.n, b.n); return t; }
inline BigInt operator%(const BigInt& a, const BigInt& b) { BigInt t; mpz_fdiv_r(t.n, a.n, b.n); return t; }
inline bool operator<(const BigInt& a, const BigInt& b) { return mpz_cmp(a.n, b.n) < 0; }
inline bool operator>(const BigInt& a, const BigInt& b) { return mpz_cmp(a.n, b.n) > 0; }
inline bool operator<=(const BigInt& a, const BigInt& b) { return mpz_cmp(a.n, b.n) <= 0; }
inline bool operator>=(const BigInt& a, const BigInt& b) { return mpz_cmp(a.n, b.n) >= 0; }
inline bool operator==(const BigInt& a, const BigInt& b) { return mpz_cmp(a.n, b.n) == 0; }
inline bool operator!=(const BigInt& a, const BigInt& b) { return mpz_cmp(a.n, b.n) != 0; }
inline bool operator<(const BigInt& a, long b) { return mpz_cmp_si(a.n, b) < 0; }
inline bool operator>(const BigInt& a, long b) { return mpz_cmp_si(a.n, b) > 0; }
inline bool operator<=(const BigInt& a, long b) { return mpz_cmp_si(a.n, b) <= 0; }
inline bool operator>=(const BigInt& a, long b) { return mpz_cmp_si(a.n, b) >= 0; }
inline bool operator==(const BigInt& a, long b) { return mpz_cmp_si(a.n, b) == 0; }
inline bool operator!=(const BigInt& a, long b) { return mpz_cmp_si(a.n, b) != 0; }

inline BigInt ident(BigInt) { return BigInt::one; }
inline BigInt zero(BigInt) { return BigInt::zero; }
template <typename U> std::enable_if_t<std::is_arithmetic_v<U>, BigInt> num(BigInt, U n) { return BigInt(n); }
inline BigInt conj(const BigInt& x) { return x; }
inline BigInt inv(const BigInt& x) { return BigInt::one / x; }
inline BigInt norm(const BigInt& x) { return x.abs(); }
inline BigInt norm2(const BigInt& x) { return x * x; }
inline int line(BigInt) { return 1; }

inline void BigInt::print() const { mpz_out_str(stdout, 10, n); }
inline void BigInt::print(int radix) const { mpz_out_str(stdout, radix, n); }
inline void BigInt::print(FILE* file) const { mpz_out_str(file, 10, n); }
inline void BigInt::print(FILE* file, int radix) const { mpz_out_str(file, radix, n); }
inline void BigInt::scan() { mpz_inp_str(n, stdin, 10); }
inline void BigInt::scan(int radix) { mpz_inp_str(n, stdin, radix); }
inline void BigInt::scan(FILE* file) { mpz_inp_str(n, file, 10); }
inline void BigInt::scan(FILE* file, int radix) { mpz_inp_str(n, file, radix); }

inline void BigInt::swap(BigInt& b) { std::swap(n, b.n); }
inline void BigInt::swap(BigInt& a, BigInt& b) { std::swap(a.n, b.n); }

inline void print(const BigInt& x, int l) { x.print(); }
inline void print(const BigInt& x, int radix, int l) { x.print(radix); }

inline std::string to_string(const BigInt& x) { char* s = mpz_get_str(nullptr, 10, x.n); std::string t(s); free(s); return t; }
inline std::string to_string(const BigInt& x, int radix) { char* s = mpz_get_str(nullptr, radix, x.n); std::string t(s); free(s); return t; }

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

inline BigInt add(const BigInt& a, const BigInt& b, const BigInt& mod) { BigInt t; mpz_add(t.n, a.n, b.n); return t >= mod ? t - mod : t; }
inline BigInt sub(const BigInt& a, const BigInt& b, const BigInt& mod) { return a >= b ? a - b : a + mod - b; }
inline uint64_t binsize(const BigInt& x) { return mpz_sizeinbase(x.n, 2); }

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
    BigFrac(const long& x, const long& y);
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
    void get_num_den(BigInt& num, BigInt& den);

    uint64_t get_den_size() const;

    void swap(BigFrac& b);
    static void swap(BigFrac& a, BigFrac& b);

    int sgn() const;
    BigFrac abs() const;
    BigFrac inv() const;
    BigFrac pow(unsigned long exp) const;
    BigFrac pow(const BigInt& exp, const BigInt& mod) const;

    BigFrac& simplify();
};

constexpr bool is_conjugate_identical(const BigFrac&) { return true; }
constexpr bool is_commutative(const BigFrac&) { return true; }
constexpr bool is_associative(const BigFrac&) { return true; }
constexpr bool is_alternative(const BigFrac&) { return true; }

constexpr bool is_unital(const BigFrac&) { return true; }
constexpr bool is_dividable(const BigFrac&) { return true; }

inline BigFrac::BigFrac() { mpq_init(n); }
inline BigFrac::BigFrac(const BigFrac& b) { mpq_init(n); mpq_set(n, b.n); }
inline BigFrac::BigFrac(const BigInt& x) { mpq_init(n); mpq_set_z(n, x.n); }
inline BigFrac::BigFrac(const BigInt& x, const BigInt& y) { mpq_init(n); mpq_set_z(n, x.n); mpq_set_den(n, y.n); }
inline BigFrac::BigFrac(const char* s) { mpq_init(n); mpq_set_str(n, s, 10); }
inline BigFrac::BigFrac(const char* s, int radix) { mpq_init(n); mpq_set_str(n, s, radix); }
inline BigFrac::BigFrac(const int& x) { mpq_init(n); mpq_set_si(n, x, 1); }
inline BigFrac::BigFrac(const unsigned int& x) { mpq_init(n); mpq_set_ui(n, x, 1); }
inline BigFrac::BigFrac(const long& x) { mpq_init(n); mpq_set_si(n, x, 1); }
inline BigFrac::BigFrac(const long& x, const long& y) { mpq_init(n); mpq_set_si(n, x, y); }
inline BigFrac::BigFrac(const unsigned long& x) { mpq_init(n); mpq_set_ui(n, x, 1); }
inline BigFrac::BigFrac(const long long& x) { mpq_init(n); mpq_set_si(n, x, 1); }
inline BigFrac::BigFrac(const unsigned long long& x) { mpq_init(n); mpq_set_ui(n, x, 1); }
inline BigFrac::~BigFrac() { mpq_clear(n); }
inline BigFrac& BigFrac::operator=(const BigFrac& b) { mpq_set(n, b.n); return *this; }
inline BigFrac& BigFrac::operator=(const char* s) { mpq_set_str(n, s, 10); return *this; }
inline BigFrac& BigFrac::operator=(const int& x) { mpq_set_si(n, x, 1); return *this; }
inline BigFrac& BigFrac::operator=(const unsigned int& x) { mpq_set_ui(n, x, 1); return *this; }
inline BigFrac& BigFrac::operator=(const long& x) { mpq_set_si(n, x, 1); return *this; }
inline BigFrac& BigFrac::operator=(const unsigned long& x) { mpq_set_ui(n, x, 1); return *this; }
inline BigFrac& BigFrac::operator=(const long long& x) { mpq_set_si(n, x, 1); return *this; }
inline BigFrac& BigFrac::operator=(const unsigned long long& x) { mpq_set_ui(n, x, 1); return *this; }
inline BigFrac& BigFrac::operator++() { mpq_add(n, n, one.n); return *this; }
inline BigFrac& BigFrac::operator--() { mpq_sub(n, n, one.n); return *this; }
inline BigFrac BigFrac::operator++(int) { BigFrac t; mpq_add(t.n, n, one.n); return t; }
inline BigFrac BigFrac::operator--(int) { BigFrac t; mpq_sub(t.n, n, one.n); return t; }
inline BigFrac BigFrac::operator-() const { BigFrac t = *this; mpq_neg(t.n, t.n); return t; }
inline BigFrac BigFrac::operator+() const { return *this; }

inline BigFrac& operator+=(BigFrac& a, const BigFrac& b) { mpq_add(a.n, a.n, b.n); return a; }
inline BigFrac& operator-=(BigFrac& a, const BigFrac& b) { mpq_sub(a.n, a.n, b.n); return a; }
inline BigFrac& operator*=(BigFrac& a, const BigFrac& b) { mpq_mul(a.n, a.n, b.n); return a; }
inline BigFrac& operator/=(BigFrac& a, const BigFrac& b) { mpq_div(a.n, a.n, b.n); return a; }
inline BigFrac& operator<<=(BigFrac& a, const unsigned long& b) { mpq_mul_2exp(a.n, a.n, b); return a; }
inline BigFrac& operator>>=(BigFrac& a, const unsigned long& b) { mpq_div_2exp(a.n, a.n, b); return a; }
inline BigFrac operator+(const BigFrac& a, const BigFrac& b) { BigFrac t; mpq_add(t.n, a.n, b.n); return t; }
inline BigFrac operator-(const BigFrac& a, const BigFrac& b) { BigFrac t; mpq_sub(t.n, a.n, b.n); return t; }
inline BigFrac operator*(const BigFrac& a, const BigFrac& b) { BigFrac t; mpq_mul(t.n, a.n, b.n); return t; }
inline BigFrac operator/(const BigFrac& a, const BigFrac& b) { BigFrac t; mpq_div(t.n, a.n, b.n); return t; }
inline BigFrac operator<<(const BigFrac& a, const unsigned long& b) { BigFrac t; mpq_mul_2exp(t.n, a.n, b); return t; }
inline BigFrac operator>>(const BigFrac& a, const unsigned long& b) { BigFrac t; mpq_div_2exp(t.n, a.n, b); return t; }
inline bool operator<(const BigFrac& a, const BigFrac& b) { return mpq_cmp(a.n, b.n) < 0; }
inline bool operator>(const BigFrac& a, const BigFrac& b) { return mpq_cmp(a.n, b.n) > 0; }
inline bool operator<=(const BigFrac& a, const BigFrac& b) { return mpq_cmp(a.n, b.n) <= 0; }
inline bool operator>=(const BigFrac& a, const BigFrac& b) { return mpq_cmp(a.n, b.n) >= 0; }
inline bool operator==(const BigFrac& a, const BigFrac& b) { return mpq_cmp(a.n, b.n) == 0; }
inline bool operator!=(const BigFrac& a, const BigFrac& b) { return mpq_cmp(a.n, b.n) != 0; }

inline BigFrac::operator bool() const { return mpq_cmp_ui(n, 0, 1) != 0; }
inline BigFrac::operator float() const { return (float)mpq_get_d(n); }
inline BigFrac::operator double() const { return mpq_get_d(n); }
inline BigFrac::operator long double() const { return mpq_get_d(n); }

inline BigFrac ident(BigFrac) { return BigFrac::one; }
inline BigFrac zero(BigFrac) { return BigFrac::zero; }
template <typename U> std::enable_if_t<std::is_arithmetic_v<U>, BigFrac> num(BigFrac, U n) { return BigFrac(n); }
inline BigFrac conj(const BigFrac& x) { return x; }
inline BigFrac inv(const BigFrac& x) { return x.inv(); }
inline BigFrac norm(const BigFrac& x) { BigFrac t = x; t.abs(); return t; }
inline BigFrac norm2(const BigFrac& x) { return x * x; }
inline int line(BigFrac) { return 1; }

inline void BigFrac::print() const { mpq_out_str(stdout, 10, n); }
inline void BigFrac::print(int radix) const { mpq_out_str(stdout, radix, n); }
inline void BigFrac::print(FILE* file) const { mpq_out_str(file, 10, n); }
inline void BigFrac::print(FILE* file, int radix) const { mpq_out_str(file, radix, n); }
inline void BigFrac::scan() { mpq_inp_str(n, stdin, 10); }
inline void BigFrac::scan(int radix) { mpq_inp_str(n, stdin, radix); }
inline void BigFrac::scan(FILE* file) { mpq_inp_str(n, file, 10); }
inline void BigFrac::scan(FILE* file, int radix) { mpq_inp_str(n, file, radix); }

inline std::string to_string(const BigFrac& x) { char* s = mpq_get_str(nullptr, 10, x.n); std::string t(s); free(s); return t; }
inline std::string to_string(const BigFrac& x, int radix) { char* s = mpq_get_str(nullptr, radix, x.n); std::string t(s); free(s); return t; }

inline void BigFrac::get_num(BigInt& num) { mpz_set(num.n, mpq_numref(n)); }
inline void BigFrac::get_den(BigInt& den) { mpz_set(den.n, mpq_denref(n)); }
inline void BigFrac::get_num_den(BigInt& num, BigInt& den) { mpz_set(num.n, mpq_numref(n)); mpz_set(den.n, mpq_denref(n)); }

inline uint64_t BigFrac::get_den_size() const { return mpz_sizeinbase(mpq_denref(n), 2); }

inline void BigFrac::swap(BigFrac& b) { mpq_swap(n, b.n); }
inline void BigFrac::swap(BigFrac& a, BigFrac& b) { mpq_swap(a.n, b.n); }

inline void print(const BigFrac& x, int l) { x.print(); }
inline void print(const BigFrac& x, int radix, int l) { x.print(radix); }

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
    BigFloat(const BigFrac& x, mpfr_prec_t prec);
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
    BigFloat(const float& x, mpfr_prec_t prec);
    BigFloat(const double& x, mpfr_prec_t prec);
    BigFloat(const long double& x, mpfr_prec_t prec);
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

    void set_prec(mpfr_prec_t prec);
    mpfr_prec_t get_prec() const;

    void print() const;
    void print(int radix) const;
    void print(FILE* file) const;
    void print(FILE* file, int radix) const;
    void scan();
    void scan(int radix);
    void scan(FILE* file);
    void scan(FILE* file, int radix);

    void print(mpfr_prec_t prec) const;
    void print(mpfr_prec_t prec, int radix) const;
    void print(mpfr_prec_t prec, FILE* file) const;
    void print(mpfr_prec_t prec, FILE* file, int radix) const;

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

inline BigFloat& operator+=(BigFloat& a, const BigFloat& b) { mpfr_add(a.n, a.n, b.n, MPFR_RNDD); return a; }
inline BigFloat& operator-=(BigFloat& a, const BigFloat& b) { mpfr_sub(a.n, a.n, b.n, MPFR_RNDD); return a; }
inline BigFloat& operator*=(BigFloat& a, const BigFloat& b) { mpfr_mul(a.n, a.n, b.n, MPFR_RNDD); return a; }
inline BigFloat& operator/=(BigFloat& a, const BigFloat& b) { mpfr_div(a.n, a.n, b.n, MPFR_RNDD); return a; }
inline BigFloat& operator%=(BigFloat& a, const BigFloat& b) { mpfr_fmod(a.n, a.n, b.n, MPFR_RNDD); return a; }
inline BigFloat& operator<<=(BigFloat& a, const unsigned long& b) { mpfr_mul_2exp(a.n, a.n, b, MPFR_RNDD); return a; }
inline BigFloat& operator>>=(BigFloat& a, const unsigned long& b) { mpfr_div_2exp(a.n, a.n, b, MPFR_RNDD); return a; }
inline BigFloat operator+(const BigFloat& a, const BigFloat& b) { BigFloat t; t.set_prec(std::max(a.get_prec(), b.get_prec())); mpfr_add(t.n, a.n, b.n, MPFR_RNDD); return t; }
inline BigFloat operator-(const BigFloat& a, const BigFloat& b) { BigFloat t; t.set_prec(std::max(a.get_prec(), b.get_prec())); mpfr_sub(t.n, a.n, b.n, MPFR_RNDD); return t; }
inline BigFloat operator*(const BigFloat& a, const BigFloat& b) { BigFloat t; t.set_prec(std::max(a.get_prec(), b.get_prec())); mpfr_mul(t.n, a.n, b.n, MPFR_RNDD); return t; }
inline BigFloat operator/(const BigFloat& a, const BigFloat& b) { BigFloat t; t.set_prec(std::max(a.get_prec(), b.get_prec())); mpfr_div(t.n, a.n, b.n, MPFR_RNDD); return t; }
inline BigFloat operator%(const BigFloat& a, const BigFloat& b) { BigFloat t; t.set_prec(std::max(a.get_prec(), b.get_prec())); mpfr_fmod(t.n, a.n, b.n, MPFR_RNDD); return t; }
inline BigFloat operator<<(const BigFloat& a, const unsigned long& b) { BigFloat t(a); mpfr_mul_2exp(t.n, t.n, b, MPFR_RNDD); return t; }
inline BigFloat operator>>(const BigFloat& a, const unsigned long& b) { BigFloat t(a); mpfr_div_2exp(t.n, t.n, b, MPFR_RNDD); return t; }
inline bool operator<(const BigFloat& a, const BigFloat& b) { return mpfr_cmp(a.n, b.n) < 0; }
inline bool operator>(const BigFloat& a, const BigFloat& b) { return mpfr_cmp(a.n, b.n) > 0; }
inline bool operator<=(const BigFloat& a, const BigFloat& b) { return mpfr_cmp(a.n, b.n) <= 0; }
inline bool operator>=(const BigFloat& a, const BigFloat& b) { return mpfr_cmp(a.n, b.n) >= 0; }
inline bool operator==(const BigFloat& a, const BigFloat& b) { return mpfr_cmp(a.n, b.n) == 0; }
inline bool operator!=(const BigFloat& a, const BigFloat& b) { return mpfr_cmp(a.n, b.n) != 0; }

inline BigFloat::operator bool() const { return mpfr_cmp_ui(n, 0) != 0; }
inline BigFloat::operator float() const { return mpfr_get_d(n, MPFR_RNDD); }
inline BigFloat::operator double() const { return mpfr_get_d(n, MPFR_RNDD); }
inline BigFloat::operator long double() const { return mpfr_get_d(n, MPFR_RNDD); }

inline BigFloat ident(const BigFloat& x) { return BigFloat(1.0, x.get_prec()); }
inline BigFloat zero(const BigFloat& x) { return BigFloat(0.0, x.get_prec()); }
template <typename U> std::enable_if_t<std::is_arithmetic_v<U>, BigFloat> num(BigFloat, U n) { return BigFloat(n); }
inline BigFloat conj(const BigFloat& x) { return x; }
inline BigFloat inv(const BigFloat& x) { return BigFloat::one / x; }
inline BigFloat norm(const BigFloat& x) { BigFloat t = x; t.abs(); return t; }
inline BigFloat norm2(const BigFloat& x) { return x * x; }
inline int line(BigFloat) { return 1; }

inline void BigFloat::swap(BigFloat& b) { mpfr_swap(n, b.n); }
inline void BigFloat::swap(BigFloat& a, BigFloat& b) { mpfr_swap(a.n, b.n); }

inline void print(const BigFloat& x, int l = 0) { x.print(); }
inline void print(const BigFloat& x, int radix, int l) { x.print(radix); }

inline void print(mpfr_prec_t prec, const BigFloat& x) { x.print(prec); }
inline void print(mpfr_prec_t prec, const BigFloat& x, int radix) { x.print(prec, radix); }

inline std::string to_string(const BigFloat& x) { char* s = mpfr_get_str(nullptr, nullptr, 10, 0, x.n, MPFR_RNDD); std::string t(s); free(s); return t; }
inline std::string to_string(const BigFloat& x, int radix) { char* s = mpfr_get_str(nullptr, nullptr, radix, 0, x.n, MPFR_RNDD); std::string t(s); free(s); return t; }

int sgn(const BigInt& x);
BigFloat abs(const BigFloat& x);
BigFloat pow(const BigFloat& x, unsigned long exp);
BigFloat pow(const BigFloat& x, const BigFloat& exp);
BigFloat sqrt(const BigFloat& x);
BigFloat cbrt(const BigFloat& x);
BigFloat root(const BigFloat& x, unsigned long exp);
BigInt floor(const BigFloat& x);
BigInt ceil(const BigFloat& x);
BigInt trunc(const BigFloat& x);

BigFloat get_pi();
BigFloat get_e();
BigFloat get_pi(mpfr_prec_t prec);
BigFloat get_e(mpfr_prec_t prec);

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

inline BigInt::BigInt(const BigFrac& x) { mpz_init(n); mpz_set_q(n, x.n); }
inline BigInt::BigInt(const BigFloat& x) { mpz_init(n); mpz_set_f(n, x.n); }

inline BigInt& BigInt::operator=(const BigFrac& x) { mpz_set_q(n, x.n); return *this; }
inline BigInt& BigInt::operator=(const BigFloat& x) { mpz_set_f(n, x.n); return *this; }

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

bool check_root(BigInt x, BigInt p, BigInt pm1_prime[], uint64_t cnt);
BigInt find_root(BigInt p, BigInt pm1_prime[], uint64_t cnt);
bool has_root(BigInt pm1_prime[], uint64_t exp[], uint64_t cnt);

BigInt cpow(BigInt x, BigInt w2, BigInt n, BigInt mod);
BigInt cipolla(BigInt x, BigInt mod);

// a^x = b (mod p)
BigInt bsgs(const BigInt& a, const BigInt& b, const BigInt& p);

BigInt index_calculus_log(BigInt a, BigInt b, BigInt g, BigInt p);
BigInt index_calculus_log(BigInt g, BigInt b, BigInt p);

#endif