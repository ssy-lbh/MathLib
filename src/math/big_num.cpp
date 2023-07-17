#include "big_num.h"

const BigInt BigInt::zero(0);
const BigInt BigInt::one(1);

BigInt::BigInt() { mpz_init(n); }
BigInt::BigInt(const BigInt& b) { mpz_init_set(n, b.n); }
BigInt::BigInt(const char* s) { mpz_init_set_str(n, s, 10); }
BigInt::BigInt(const char* s, int radix) { mpz_init_set_str(n, s, radix); }
BigInt::BigInt(const int& x) { mpz_init_set_si(n, x); }
BigInt::BigInt(const unsigned int& x) { mpz_init_set_ui(n, x); }
BigInt::BigInt(const long& x) { mpz_init_set_si(n, x); }
BigInt::BigInt(const unsigned long& x) { mpz_init_set_ui(n, x); }
BigInt::BigInt(const long long& x) { mpz_init_set_si(n, x); }
BigInt::BigInt(const unsigned long long& x) { mpz_init_set_ui(n, x); }
BigInt::BigInt(const BigFrac& x) { mpz_init(n); mpz_set_q(n, x.n); }
BigInt::BigInt(const BigFloat& x) { mpz_init(n); mpz_set_f(n, x.n); }
BigInt::~BigInt() { mpz_clear(n); }
BigInt& BigInt::operator=(const BigInt& b) { mpz_set(n, b.n); return *this; }
BigInt& BigInt::operator=(const char* s) { mpz_set_str(n, s, 10); return *this; }
BigInt& BigInt::operator=(const int& x) { mpz_set_si(n, x); return *this; }
BigInt& BigInt::operator=(const unsigned int& x) { mpz_set_ui(n, x); return *this; }
BigInt& BigInt::operator=(const long& x) { mpz_set_si(n, x); return *this; }
BigInt& BigInt::operator=(const unsigned long& x) { mpz_set_ui(n, x); return *this; }
BigInt& BigInt::operator=(const long long& x) { mpz_set_si(n, x); return *this; }
BigInt& BigInt::operator=(const unsigned long long& x) { mpz_set_ui(n, x); return *this; }
BigInt& BigInt::operator=(const BigFrac& x) { mpz_set_q(n, x.n); return *this; }
BigInt& BigInt::operator=(const BigFloat& x) { mpz_set_f(n, x.n); return *this; }
BigInt& BigInt::operator+=(const BigInt& b) { mpz_add(n, n, b.n); return *this; }
BigInt& BigInt::operator-=(const BigInt& b) { mpz_sub(n, n, b.n); return *this; }
BigInt& BigInt::operator*=(const BigInt& b) { mpz_mul(n, n, b.n); return *this; }
BigInt& BigInt::operator/=(const BigInt& b) { mpz_tdiv_q(n, n, b.n); return *this; }
BigInt& BigInt::operator%=(const BigInt& b) { mpz_tdiv_r(n, n, b.n); return *this; }
BigInt& BigInt::operator<<=(const unsigned long& b) { mpz_mul_2exp(n, n, b); return *this; }
BigInt& BigInt::operator>>=(const unsigned long& b) { mpz_tdiv_q_2exp(n, n, b); return *this; }
BigInt& BigInt::operator&=(const BigInt& b) { mpz_and(n, n, b.n); return *this; }
BigInt& BigInt::operator|=(const BigInt& b) { mpz_ior(n, n, b.n); return *this; }
BigInt& BigInt::operator^=(const BigInt& b) { mpz_xor(n, n, b.n); return *this; }
BigInt& BigInt::operator++() { mpz_add_ui(n, n, 1); return *this; }
BigInt& BigInt::operator--() { mpz_sub_ui(n, n, 1); return *this; }
BigInt BigInt::operator++(int) { BigInt t = *this; mpz_add_ui(n, n, 1); return t; }
BigInt BigInt::operator--(int) { BigInt t = *this; mpz_sub_ui(n, n, 1); return t; }
BigInt BigInt::operator-() const { BigInt t = *this; mpz_neg(t.n, t.n); return t; }
BigInt BigInt::operator+() const { return *this; }
BigInt BigInt::operator<<(const unsigned long& b) const { BigInt t; mpz_mul_2exp(t.n, n, b); return t; }
BigInt BigInt::operator>>(const unsigned long& b) const { BigInt t; mpz_div_2exp(t.n, n, b); return t; }
BigInt BigInt::operator&(const BigInt& b) const { BigInt t; mpz_and(t.n, n, b.n); return t; }
BigInt BigInt::operator|(const BigInt& b) const { BigInt t; mpz_ior(t.n, n, b.n); return t; }
BigInt BigInt::operator^(const BigInt& b) const { BigInt t; mpz_xor(t.n, n, b.n); return t; }
BigInt BigInt::operator~() const { BigInt t; mpz_com(t.n, n); return t; }
BigInt BigInt::operator+(const BigInt& b) const { BigInt t; mpz_add(t.n, n, b.n); return t; }
BigInt BigInt::operator-(const BigInt& b) const { BigInt t; mpz_sub(t.n, n, b.n); return t; }
BigInt BigInt::operator*(const BigInt& b) const { BigInt t; mpz_mul(t.n, n, b.n); return t; }
BigInt BigInt::operator/(const BigInt& b) const { BigInt t; mpz_tdiv_q(t.n, n, b.n); return t; }
BigInt BigInt::operator%(const BigInt& b) const { BigInt t; mpz_tdiv_r(t.n, n, b.n); return t; }
bool BigInt::operator<(const BigInt& b) const { return mpz_cmp(n, b.n) < 0; }
bool BigInt::operator>(const BigInt& b) const { return mpz_cmp(n, b.n) > 0; }
bool BigInt::operator<=(const BigInt& b) const { return mpz_cmp(n, b.n) <= 0; }
bool BigInt::operator>=(const BigInt& b) const { return mpz_cmp(n, b.n) >= 0; }
bool BigInt::operator==(const BigInt& b) const { return mpz_cmp(n, b.n) == 0; }
bool BigInt::operator!=(const BigInt& b) const { return mpz_cmp(n, b.n) != 0; }

BigInt::operator bool() const { return mpz_cmp_ui(n, 0) != 0; }
BigInt::operator int() const { return mpz_get_si(n); }
BigInt::operator unsigned int() const { return mpz_get_ui(n); }
BigInt::operator long() const { return mpz_get_si(n); }
BigInt::operator unsigned long() const { return mpz_get_ui(n); }
BigInt::operator long long() const { return mpz_get_si(n); }
BigInt::operator unsigned long long() const { return mpz_get_ui(n); }

void BigInt::print() const { mpz_out_str(stdout, 10, n); }
void BigInt::print(int radix) const { mpz_out_str(stdout, radix, n); }
void BigInt::print(FILE* file) const { mpz_out_str(file, 10, n); }
void BigInt::print(FILE* file, int radix) const { mpz_out_str(file, radix, n); }
void BigInt::scan() { mpz_inp_str(n, stdin, 10); }
void BigInt::scan(int radix) { mpz_inp_str(n, stdin, radix); }
void BigInt::scan(FILE* file) { mpz_inp_str(n, file, 10); }
void BigInt::scan(FILE* file, int radix) { mpz_inp_str(n, file, radix); }

void BigInt::swap(BigInt& b) { mpz_swap(n, b.n); }
void BigInt::swap(BigInt& a, BigInt& b) { mpz_swap(a.n, b.n); }

int BigInt::sgn() const { return mpz_sgn(n); }
BigInt BigInt::abs() { BigInt t; mpz_abs(t.n, n); return t; }
BigInt BigInt::pow(unsigned long exp) { BigInt t; mpz_pow_ui(t.n, n, exp); return t; }
BigInt BigInt::pow(const BigInt& exp, const BigInt& mod) { BigInt t; mpz_powm(t.n, n, exp.n, mod.n); return t; }
BigInt BigInt::sqrt() { BigInt t; mpz_sqrt(t.n, n); return t; }
BigInt BigInt::sqrtrem(BigInt& rem) { BigInt t; mpz_sqrtrem(t.n, rem.n, n); return t; }
BigInt BigInt::root(unsigned long exp) { BigInt t; mpz_root(t.n, n, exp); return t; }
BigInt BigInt::rootrem(BigInt& rem, unsigned long exp) { BigInt t; mpz_rootrem(t.n, rem.n, n, exp); return t; }
BigInt BigInt::gcd(const BigInt& b) { BigInt t; mpz_gcd(t.n, n, b.n); return t; }
BigInt BigInt::lcm(const BigInt& b) { BigInt t; mpz_lcm(t.n, n, b.n); return t; }
BigInt BigInt::factorial() { BigInt t; mpz_fac_ui(t.n, mpz_get_ui(n)); return t; }
BigInt BigInt::fibonacci() { BigInt t; mpz_fib_ui(t.n, mpz_get_ui(n)); return t; }
BigInt BigInt::lucas() { BigInt t; mpz_lucnum_ui(t.n, mpz_get_ui(n)); return t; }
BigInt BigInt::binomial(const BigInt& b) { BigInt t; mpz_bin_ui(t.n, n, mpz_get_ui(b.n)); return t; }
BigInt BigInt::binomial(const BigInt& b, const BigInt& m) { BigInt t; mpz_bin_ui(t.n, n, mpz_get_ui(b.n)); mpz_mod(t.n, t.n, m.n); return t; }
BigInt BigInt::multinomial(const BigInt& b) { BigInt t; mpz_sub(t.n, n, b.n); mpz_fac_ui(t.n, mpz_get_ui(t.n)); mpz_fac_ui(t.n, mpz_get_ui(t.n)); mpz_divexact_ui(t.n, t.n, mpz_get_ui(b.n)); return t; }
BigInt BigInt::multinomial(const BigInt& b, const BigInt& m) { BigInt t; mpz_sub(t.n, n, b.n); mpz_fac_ui(t.n, mpz_get_ui(t.n)); mpz_fac_ui(t.n, mpz_get_ui(t.n)); mpz_divexact_ui(t.n, t.n, mpz_get_ui(b.n)); mpz_mod(t.n, t.n, m.n); return t; }
bool BigInt::is_prime() const { return mpz_probab_prime_p(n, 25) > 0; }
bool BigInt::is_probab_prime(int reps) const { return mpz_probab_prime_p(n, reps) > 0; }
bool BigInt::is_perfect_power() const { return mpz_perfect_power_p(n) > 0; }
bool BigInt::is_perfect_square() const { return mpz_perfect_square_p(n) > 0; }
bool BigInt::testbit(unsigned long bit) const { return mpz_tstbit(n, bit) > 0; }
void BigInt::setbit(unsigned long bit) { mpz_setbit(n, bit); }
void BigInt::clrbit(unsigned long bit) { mpz_clrbit(n, bit); }
void BigInt::combit(unsigned long bit) { mpz_combit(n, bit); }
unsigned long BigInt::scanbit0(unsigned long bit) { return mpz_scan0(n, bit); }
unsigned long BigInt::scanbit1(unsigned long bit) { return mpz_scan1(n, bit); }
void BigInt::gcdext(const BigInt& b, BigInt& s, BigInt& t) { mpz_gcdext(n, s.n, t.n, n, b.n); }
void BigInt::nextprime() { mpz_nextprime(n, n); }

int sgn(const BigInt& x) { return mpz_sgn(x.n); }
BigInt abs(const BigInt& x) { BigInt t; mpz_abs(t.n, x.n); return t; }
BigInt mul(const BigInt& x, const BigInt& y, const BigInt& mod) { BigInt t; mpz_mul(t.n, x.n, y.n); mpz_mod(t.n, t.n, mod.n); return t; }
BigInt pow(const BigInt& x, unsigned long exp) { BigInt t; mpz_pow_ui(t.n, x.n, exp); return t; }
BigInt pow(const BigInt& x, const BigInt& exp, const BigInt& mod) { BigInt t; mpz_powm(t.n, x.n, exp.n, mod.n); return t; }
BigInt sqrt(const BigInt& x) { BigInt t; mpz_sqrt(t.n, x.n); return t; }
BigInt sqrtrem(const BigInt& x, BigInt& rem) { BigInt t; mpz_sqrtrem(t.n, rem.n, x.n); return t; }
BigInt root(const BigInt& x, unsigned long exp) { BigInt t; mpz_root(t.n, x.n, exp); return t; }
BigInt rootrem(const BigInt& x, BigInt& rem, unsigned long exp) { BigInt t; mpz_rootrem(t.n, rem.n, x.n, exp); return t; }
BigInt gcd(const BigInt& x, const BigInt& y) { BigInt t; mpz_gcd(t.n, x.n, y.n); return t; }
BigInt lcm(const BigInt& x, const BigInt& y) { BigInt t; mpz_lcm(t.n, x.n, y.n); return t; }
BigInt exgcd(const BigInt& x, const BigInt& y, BigInt& s, BigInt& t) { BigInt t; mpz_gcdext(t.n, s.n, t.n, x.n, y.n); return t; }
BigInt factorial(const BigInt& x) { BigInt t; mpz_fac_ui(t.n, mpz_get_ui(x.n)); return t; }
BigInt fibonacci(const BigInt& x) { BigInt t; mpz_fib_ui(t.n, mpz_get_ui(x.n)); return t; }
BigInt lucas(const BigInt& x) { BigInt t; mpz_lucnum_ui(t.n, mpz_get_ui(x.n)); return t; }
BigInt binomial(const BigInt& x, const BigInt& y) { BigInt t; mpz_bin_ui(t.n, x.n, mpz_get_ui(y.n)); return t; }
BigInt binomial(const BigInt& x, const BigInt& y, const BigInt& m) { BigInt t; mpz_bin_ui(t.n, x.n, mpz_get_ui(y.n)); mpz_mod(t.n, t.n, m.n); return t; }
BigInt multinomial(const BigInt& x, const BigInt& y) { BigInt t; mpz_sub(t.n, x.n, y.n); mpz_fac_ui(t.n, mpz_get_ui(t.n)); mpz_fac_ui(t.n, mpz_get_ui(t.n)); mpz_divexact_ui(t.n, t.n, mpz_get_ui(y.n)); return t; }
BigInt multinomial(const BigInt& x, const BigInt& y, const BigInt& m) { BigInt t; mpz_sub(t.n, x.n, y.n); mpz_fac_ui(t.n, mpz_get_ui(t.n)); mpz_fac_ui(t.n, mpz_get_ui(t.n)); mpz_divexact_ui(t.n, t.n, mpz_get_ui(y.n)); mpz_mod(t.n, t.n, m.n); return t; }
int legendre(const BigInt& x, const BigInt& p) { return mpz_legendre(x.n, p.n); }
int jacobi(const BigInt& x, const BigInt& p) { return mpz_jacobi(x.n, p.n); }

const BigFrac BigFrac::zero(0);
const BigFrac BigFrac::one(1);

BigFrac::BigFrac() { mpq_init(n); }
BigFrac::BigFrac(const BigFrac& b) { mpq_init(n); mpq_set(n, b.n); }
BigFrac::BigFrac(const BigInt& x) { mpq_init(n); mpq_set_z(n, x.n); }
BigFrac::BigFrac(const BigInt& x, const BigInt& y) { mpq_init(n); mpq_set_z(n, x.n); mpq_set_den(n, y.n); }
BigFrac::BigFrac(const char* s) { mpq_init(n); mpq_set_str(n, s, 10); }
BigFrac::BigFrac(const char* s, int radix) { mpq_init(n); mpq_set_str(n, s, radix); }
BigFrac::BigFrac(const int& x) { mpq_init(n); mpq_set_si(n, x, 1); }
BigFrac::BigFrac(const unsigned int& x) { mpq_init(n); mpq_set_ui(n, x, 1); }
BigFrac::BigFrac(const long& x) { mpq_init(n); mpq_set_si(n, x, 1); }
BigFrac::BigFrac(const unsigned long& x) { mpq_init(n); mpq_set_ui(n, x, 1); }
BigFrac::BigFrac(const long long& x) { mpq_init(n); mpq_set_si(n, x, 1); }
BigFrac::BigFrac(const unsigned long long& x) { mpq_init(n); mpq_set_ui(n, x, 1); }
BigFrac::~BigFrac() { mpq_clear(n); }
BigFrac& BigFrac::operator=(const BigFrac& b) { mpq_set(n, b.n); return *this; }
BigFrac& BigFrac::operator=(const char* s) { mpq_set_str(n, s, 10); return *this; }
BigFrac& BigFrac::operator=(const int& x) { mpq_set_si(n, x, 1); return *this; }
BigFrac& BigFrac::operator=(const unsigned int& x) { mpq_set_ui(n, x, 1); return *this; }
BigFrac& BigFrac::operator=(const long& x) { mpq_set_si(n, x, 1); return *this; }
BigFrac& BigFrac::operator=(const unsigned long& x) { mpq_set_ui(n, x, 1); return *this; }
BigFrac& BigFrac::operator=(const long long& x) { mpq_set_si(n, x, 1); return *this; }
BigFrac& BigFrac::operator=(const unsigned long long& x) { mpq_set_ui(n, x, 1); return *this; }
BigFrac& BigFrac::operator+=(const BigFrac& b) { mpq_add(n, n, b.n); return *this; }
BigFrac& BigFrac::operator-=(const BigFrac& b) { mpq_sub(n, n, b.n); return *this; }
BigFrac& BigFrac::operator*=(const BigFrac& b) { mpq_mul(n, n, b.n); return *this; }
BigFrac& BigFrac::operator/=(const BigFrac& b) { mpq_div(n, n, b.n); return *this; }
BigFrac& BigFrac::operator<<=(const unsigned long& b) { mpq_mul_2exp(n, n, b); return *this; }
BigFrac& BigFrac::operator>>=(const unsigned long& b) { mpq_div_2exp(n, n, b); return *this; }
BigFrac& BigFrac::operator++() { mpq_add(n, n, one.n); return *this; }
BigFrac& BigFrac::operator--() { mpq_sub(n, n, one.n); return *this; }
BigFrac BigFrac::operator++(int) { BigFrac t; mpq_add(t.n, n, one.n); return t; }
BigFrac BigFrac::operator--(int) { BigFrac t; mpq_sub(t.n, n, one.n); return t; }
BigFrac BigFrac::operator-() const { BigFrac t = *this; mpq_neg(t.n, t.n); return t; }
BigFrac BigFrac::operator+() const { return *this; }
BigFrac BigFrac::operator+(const BigFrac& b) const { BigFrac t; mpq_add(t.n, n, b.n); return t; }
BigFrac BigFrac::operator-(const BigFrac& b) const { BigFrac t; mpq_sub(t.n, n, b.n); return t; }
BigFrac BigFrac::operator*(const BigFrac& b) const { BigFrac t; mpq_mul(t.n, n, b.n); return t; }
BigFrac BigFrac::operator/(const BigFrac& b) const { BigFrac t; mpq_div(t.n, n, b.n); return t; }
BigFrac BigFrac::operator<<(const unsigned long& b) const { BigFrac t; mpq_mul_2exp(t.n, n, b); return t; }
BigFrac BigFrac::operator>>(const unsigned long& b) const { BigFrac t; mpq_div_2exp(t.n, n, b); return t; }
bool BigFrac::operator<(const BigFrac& b) const { return mpq_cmp(n, b.n) < 0; }
bool BigFrac::operator>(const BigFrac& b) const { return mpq_cmp(n, b.n) > 0; }
bool BigFrac::operator<=(const BigFrac& b) const { return mpq_cmp(n, b.n) <= 0; }
bool BigFrac::operator>=(const BigFrac& b) const { return mpq_cmp(n, b.n) >= 0; }
bool BigFrac::operator==(const BigFrac& b) const { return mpq_cmp(n, b.n) == 0; }
bool BigFrac::operator!=(const BigFrac& b) const { return mpq_cmp(n, b.n) != 0; }

BigFrac::operator bool() const { return mpq_cmp_ui(n, 0, 1) != 0; }
BigFrac::operator double() const { return mpq_get_d(n); }
BigFrac::operator long double() const { return mpq_get_d(n); }

void BigFrac::print() const { mpq_out_str(stdout, 10, n); }
void BigFrac::print(int radix) const { mpq_out_str(stdout, radix, n); }
void BigFrac::print(FILE* file) const { mpq_out_str(file, 10, n); }
void BigFrac::print(FILE* file, int radix) const { mpq_out_str(file, radix, n); }
void BigFrac::scan() { mpq_inp_str(n, stdin, 10); }
void BigFrac::scan(int radix) { mpq_inp_str(n, stdin, radix); }
void BigFrac::scan(FILE* file) { mpq_inp_str(n, file, 10); }
void BigFrac::scan(FILE* file, int radix) { mpq_inp_str(n, file, radix); }

void BigFrac::get_num(BigInt& num) { mpz_set(num.n, mpq_numref(n)); }
void BigFrac::get_den(BigInt& den) { mpz_set(den.n, mpq_denref(n)); }

void BigFrac::swap(BigFrac& b) { mpq_swap(n, b.n); }
void BigFrac::swap(BigFrac& a, BigFrac& b) { mpq_swap(a.n, b.n); }

int BigFrac::sgn() const { return mpq_sgn(n); }
BigFrac BigFrac::abs() { BigFrac t; mpq_abs(t.n, n); return t; }
BigFrac BigFrac::inv() { BigFrac t; mpq_inv(t.n, n); return t; }
BigFrac BigFrac::pow(unsigned long exp) { BigFrac t; mpz_pow_ui(mpq_numref(t.n), mpq_numref(n), exp); mpz_pow_ui(mpq_denref(t.n), mpq_denref(n), exp); return t; }
BigFrac BigFrac::pow(const BigInt& exp, const BigInt& mod) { BigFrac t; mpz_powm(mpq_numref(t.n), mpq_numref(n), exp.n, mod.n); mpz_powm(mpq_denref(t.n), mpq_denref(n), exp.n, mod.n); return t; }

constexpr mpfr_prec_t precision = 100;
constexpr size_t outprec = 10;

const BigFloat BigFloat::zero(0);
const BigFloat BigFloat::one(1);

const BigFloat BigFloat::pi(PI);
const BigFloat BigFloat::e(E);

BigFloat::BigFloat() { mpfr_init2(n, precision); }
BigFloat::BigFloat(const BigFloat& b) { mpfr_init2(n, precision); mpfr_set(n, b.n, MPFR_RNDD); }
BigFloat::BigFloat(const BigInt& x) { mpfr_init2(n, precision); mpfr_set_z(n, x.n, MPFR_RNDD); }
BigFloat::BigFloat(const BigFrac& x) { mpfr_init2(n, precision); mpfr_set_q(n, x.n, MPFR_RNDD); }
BigFloat::BigFloat(const char* s) { mpfr_init2(n, precision); mpfr_set_str(n, s, 10, MPFR_RNDD); }
BigFloat::BigFloat(const char* s, int radix) { mpfr_init2(n, precision); mpfr_set_str(n, s, radix, MPFR_RNDD); }
BigFloat::BigFloat(const int& x) { mpfr_init2(n, precision); mpfr_set_si(n, x, MPFR_RNDD); }
BigFloat::BigFloat(const unsigned int& x) { mpfr_init2(n, precision); mpfr_set_ui(n, x, MPFR_RNDD); }
BigFloat::BigFloat(const long& x) { mpfr_init2(n, precision); mpfr_set_si(n, x, MPFR_RNDD); }
BigFloat::BigFloat(const unsigned long& x) { mpfr_init2(n, precision); mpfr_set_ui(n, x, MPFR_RNDD); }
BigFloat::BigFloat(const long long& x) { mpfr_init2(n, precision); mpfr_set_si(n, x, MPFR_RNDD); }
BigFloat::BigFloat(const unsigned long long& x) { mpfr_init2(n, precision); mpfr_set_ui(n, x, MPFR_RNDD); }
BigFloat::BigFloat(const float& x) { mpfr_init2(n, precision); mpfr_set_d(n, x, MPFR_RNDD); }
BigFloat::BigFloat(const double& x) { mpfr_init2(n, precision); mpfr_set_d(n, x, MPFR_RNDD); }
BigFloat::BigFloat(const long double& x) { mpfr_init2(n, precision); mpfr_set_d(n, x, MPFR_RNDD); }
BigFloat::~BigFloat() { mpfr_clear(n); }
BigFloat& BigFloat::operator=(const BigFloat& b) { mpfr_set(n, b.n, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const char* s) { mpfr_set_str(n, s, 10, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const int& x) { mpfr_set_si(n, x, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const unsigned int& x) { mpfr_set_ui(n, x, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const long& x) { mpfr_set_si(n, x, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const unsigned long& x) { mpfr_set_ui(n, x, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const long long& x) { mpfr_set_si(n, x, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const unsigned long long& x) { mpfr_set_ui(n, x, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const float& x) { mpfr_set_d(n, x, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const double& x) { mpfr_set_d(n, x, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const long double& x) { mpfr_set_d(n, x, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const mpz_t& x) { mpfr_set_z(n, x, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const mpq_t& x) { mpfr_set_q(n, x, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const mpf_t& x) { mpfr_set(n, x, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const BigInt& x) { mpfr_set_z(n, x.n, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator=(const BigFrac& x) { mpfr_set_q(n, x.n, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator+=(const BigFloat& b) { mpfr_add(n, n, b.n, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator-=(const BigFloat& b) { mpfr_sub(n, n, b.n, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator*=(const BigFloat& b) { mpfr_mul(n, n, b.n, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator/=(const BigFloat& b) { mpfr_div(n, n, b.n, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator%=(const BigFloat& b) { mpfr_fmod(n, n, b.n, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator<<=(const unsigned long& b) { mpfr_mul_2exp(n, n, b, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator>>=(const unsigned long& b) { mpfr_div_2exp(n, n, b, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator++() { mpfr_add_ui(n, n, 1, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator--() { mpfr_sub_ui(n, n, 1, MPFR_RNDD); return *this; }
BigFloat BigFloat::operator++(int) { BigFloat t = *this; mpfr_add_ui(n, n, 1, MPFR_RNDD); return t; }
BigFloat BigFloat::operator--(int) { BigFloat t = *this; mpfr_sub_ui(n, n, 1, MPFR_RNDD); return t; }
BigFloat BigFloat::operator-() const { BigFloat t = *this; mpfr_neg(t.n, t.n, MPFR_RNDD); return t; }
BigFloat BigFloat::operator+() const { return *this; }
BigFloat BigFloat::operator+(const BigFloat& b) const { BigFloat t; mpfr_add(t.n, n, b.n, MPFR_RNDD); return t; }
BigFloat BigFloat::operator-(const BigFloat& b) const { BigFloat t; mpfr_sub(t.n, n, b.n, MPFR_RNDD); return t; }
BigFloat BigFloat::operator*(const BigFloat& b) const { BigFloat t; mpfr_mul(t.n, n, b.n, MPFR_RNDD); return t; }
BigFloat BigFloat::operator/(const BigFloat& b) const { BigFloat t; mpfr_div(t.n, n, b.n, MPFR_RNDD); return t; }
BigFloat BigFloat::operator%(const BigFloat& b) const { BigFloat t; mpfr_fmod(t.n, n, b.n, MPFR_RNDD); return t; }
BigFloat BigFloat::operator<<(const unsigned long& b) const { BigFloat t; mpfr_mul_2exp(t.n, n, b, MPFR_RNDD); return t; }
BigFloat BigFloat::operator>>(const unsigned long& b) const { BigFloat t; mpfr_div_2exp(t.n, n, b, MPFR_RNDD); return t; }
bool BigFloat::operator<(const BigFloat& b) const { return mpfr_cmp(n, b.n) < 0; }
bool BigFloat::operator>(const BigFloat& b) const { return mpfr_cmp(n, b.n) > 0; }
bool BigFloat::operator<=(const BigFloat& b) const { return mpfr_cmp(n, b.n) <= 0; }
bool BigFloat::operator>=(const BigFloat& b) const { return mpfr_cmp(n, b.n) >= 0; }
bool BigFloat::operator==(const BigFloat& b) const { return mpfr_cmp(n, b.n) == 0; }
bool BigFloat::operator!=(const BigFloat& b) const { return mpfr_cmp(n, b.n) != 0; }

BigFloat::operator bool() const { return mpfr_cmp_ui(n, 0) != 0; }
BigFloat::operator float() const { return mpfr_get_d(n, MPFR_RNDD); }
BigFloat::operator double() const { return mpfr_get_d(n, MPFR_RNDD); }
BigFloat::operator long double() const { return mpfr_get_d(n, MPFR_RNDD); }

void BigFloat::print() const { mpfr_out_str(stdout, 10, outprec, n, MPFR_RNDD); }
void BigFloat::print(int radix) const { mpfr_out_str(stdout, radix, outprec, n, MPFR_RNDD); }
void BigFloat::print(FILE* file) const { mpfr_out_str(file, 10, outprec, n, MPFR_RNDD); }
void BigFloat::print(FILE* file, int radix) const { mpfr_out_str(file, radix, outprec, n, MPFR_RNDD); }
void BigFloat::scan() { mpfr_inp_str(n, stdin, 10, MPFR_RNDD); }
void BigFloat::scan(int radix) { mpfr_inp_str(n, stdin, radix, MPFR_RNDD); }
void BigFloat::scan(FILE* file) { mpfr_inp_str(n, file, 10, MPFR_RNDD); }
void BigFloat::scan(FILE* file, int radix) { mpfr_inp_str(n, file, radix, MPFR_RNDD); }

void BigFloat::swap(BigFloat& b) { mpfr_swap(n, b.n); }
void BigFloat::swap(BigFloat& a, BigFloat& b) { mpfr_swap(a.n, b.n); }

int BigFloat::sgn() const { return mpfr_sgn(n); }
BigFloat BigFloat::abs() { BigFloat t; mpfr_abs(t.n, n, MPFR_RNDD); return t; }
BigFloat BigFloat::pow(unsigned long exp) { BigFloat t; mpfr_pow_ui(t.n, n, exp, MPFR_RNDD); return t; }
BigFloat BigFloat::sqrt() { BigFloat t; mpfr_sqrt(t.n, n, MPFR_RNDD); return t; }
BigFloat BigFloat::floor() { BigFloat t; mpfr_floor(t.n, n); return t; }
BigFloat BigFloat::ceil() { BigFloat t; mpfr_ceil(t.n, n); return t; }
BigFloat BigFloat::trunc() { BigFloat t; mpfr_trunc(t.n, n); return t; }

int sgn(const BigFloat& x) { return mpfr_sgn(x.n); }
BigFloat abs(const BigFloat& x) { BigFloat t; mpfr_abs(t.n, x.n, MPFR_RNDD); return t; }
BigFloat pow(const BigFloat& x, unsigned long exp) { BigFloat t; mpfr_pow_ui(t.n, x.n, exp, MPFR_RNDD); return t; }
BigFloat pow(const BigFloat& x, const BigFloat& exp) { BigFloat t; mpfr_pow(t.n, x.n, exp.n, MPFR_RNDD); return t; }
BigFloat sqrt(const BigFloat& x) { BigFloat t; mpfr_sqrt(t.n, x.n, MPFR_RNDD); return t; }
BigFloat cbrt(const BigFloat& x) { BigFloat t; mpfr_cbrt(t.n, x.n, MPFR_RNDD); return t; }
BigFloat floor(const BigFloat& x) { BigFloat t; mpfr_floor(t.n, x.n); return t; }
BigFloat ceil(const BigFloat& x) { BigFloat t; mpfr_ceil(t.n, x.n); return t; }
BigFloat trunc(const BigFloat& x) { BigFloat t; mpfr_trunc(t.n, x.n); return t; }

BigFloat exp(const BigFloat& x) { BigFloat t; mpfr_exp(t.n, x.n, MPFR_RNDD); return t; }
BigFloat exp2(const BigFloat& x) { BigFloat t; mpfr_exp2(t.n, x.n, MPFR_RNDD); return t; }
BigFloat exp10(const BigFloat& x) { BigFloat t; mpfr_exp10(t.n, x.n, MPFR_RNDD); return t; }
BigFloat expm1(const BigFloat& x) { BigFloat t; mpfr_expm1(t.n, x.n, MPFR_RNDD); return t; }
BigFloat log(const BigFloat& x) { BigFloat t; mpfr_log(t.n, x.n, MPFR_RNDD); return t; }
BigFloat log2(const BigFloat& x) { BigFloat t; mpfr_log2(t.n, x.n, MPFR_RNDD); return t; }
BigFloat log10(const BigFloat& x) { BigFloat t; mpfr_log10(t.n, x.n, MPFR_RNDD); return t; }
BigFloat log1p(const BigFloat& x) { BigFloat t; mpfr_log1p(t.n, x.n, MPFR_RNDD); return t; }
BigFloat sin(const BigFloat& x) { BigFloat t; mpfr_sin(t.n, x.n, MPFR_RNDD); return t; }
BigFloat cos(const BigFloat& x) { BigFloat t; mpfr_cos(t.n, x.n, MPFR_RNDD); return t; }
BigFloat tan(const BigFloat& x) { BigFloat t; mpfr_tan(t.n, x.n, MPFR_RNDD); return t; }
BigFloat sinh(const BigFloat& x) { BigFloat t; mpfr_sinh(t.n, x.n, MPFR_RNDD); return t; }
BigFloat cosh(const BigFloat& x) { BigFloat t; mpfr_cosh(t.n, x.n, MPFR_RNDD); return t; }
BigFloat tanh(const BigFloat& x) { BigFloat t; mpfr_tanh(t.n, x.n, MPFR_RNDD); return t; }
BigFloat asin(const BigFloat& x) { BigFloat t; mpfr_asin(t.n, x.n, MPFR_RNDD); return t; }
BigFloat acos(const BigFloat& x) { BigFloat t; mpfr_acos(t.n, x.n, MPFR_RNDD); return t; }
BigFloat atan(const BigFloat& x) { BigFloat t; mpfr_atan(t.n, x.n, MPFR_RNDD); return t; }
BigFloat asinh(const BigFloat& x) { BigFloat t; mpfr_asinh(t.n, x.n, MPFR_RNDD); return t; }
BigFloat acosh(const BigFloat& x) { BigFloat t; mpfr_acosh(t.n, x.n, MPFR_RNDD); return t; }
BigFloat atanh(const BigFloat& x) { BigFloat t; mpfr_atanh(t.n, x.n, MPFR_RNDD); return t; }
BigFloat atan2(const BigFloat& y, const BigFloat& x) { BigFloat t; mpfr_atan2(t.n, y.n, x.n, MPFR_RNDD); return t; }
BigFloat gamma(const BigFloat& x) { BigFloat t; mpfr_gamma(t.n, x.n, MPFR_RNDD); return t; }
