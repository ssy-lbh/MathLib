#include "big_num.h"

#include <random>
#include <ctime>

const BigInt BigInt::zero(0);
const BigInt BigInt::one(1);

int BigInt::sgn() const { return mpz_sgn(n); }
BigInt BigInt::abs() const { BigInt t; mpz_abs(t.n, n); return t; }
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
BigInt exgcd(const BigInt& x, const BigInt& y, BigInt& a, BigInt& b) { BigInt t; mpz_gcdext(t.n, a.n, b.n, x.n, y.n); return t; }
BigInt factorial(const BigInt& x) { BigInt t; mpz_fac_ui(t.n, mpz_get_ui(x.n)); return t; }
BigInt fibonacci(const BigInt& x) { BigInt t; mpz_fib_ui(t.n, mpz_get_ui(x.n)); return t; }
BigInt lucas(const BigInt& x) { BigInt t; mpz_lucnum_ui(t.n, mpz_get_ui(x.n)); return t; }
BigInt binomial(const BigInt& x, const BigInt& y) { BigInt t; mpz_bin_ui(t.n, x.n, mpz_get_ui(y.n)); return t; }
BigInt binomial(const BigInt& x, const BigInt& y, const BigInt& m) { BigInt t; mpz_bin_ui(t.n, x.n, mpz_get_ui(y.n)); mpz_mod(t.n, t.n, m.n); return t; }
BigInt multinomial(const BigInt& x, const BigInt& y) { BigInt t; mpz_sub(t.n, x.n, y.n); mpz_fac_ui(t.n, mpz_get_ui(t.n)); mpz_fac_ui(t.n, mpz_get_ui(t.n)); mpz_divexact_ui(t.n, t.n, mpz_get_ui(y.n)); return t; }
BigInt multinomial(const BigInt& x, const BigInt& y, const BigInt& m) { BigInt t; mpz_sub(t.n, x.n, y.n); mpz_fac_ui(t.n, mpz_get_ui(t.n)); mpz_fac_ui(t.n, mpz_get_ui(t.n)); mpz_divexact_ui(t.n, t.n, mpz_get_ui(y.n)); mpz_mod(t.n, t.n, m.n); return t; }
int legendre(const BigInt& x, const BigInt& p) { return mpz_legendre(x.n, p.n); }
int jacobi(const BigInt& x, const BigInt& p) { return mpz_jacobi(x.n, p.n); }
int kronecker(const BigInt& x, const BigInt& p) { return mpz_kronecker(x.n, p.n); }
BigInt nextprime(const BigInt& x) { BigInt t; mpz_nextprime(t.n, x.n); return t; }
uint64_t size(const BigInt& x) { return mpz_size(x.n); }
uint64_t sizeinbase(const BigInt& x, int base) { return mpz_sizeinbase(x.n, base); }
BigInt fdivq(const BigInt& n, const BigInt& d) { BigInt t; mpz_fdiv_q(t.n, n.n, d.n); return t; }
BigInt fdivr(const BigInt& n, const BigInt& d) { BigInt t; mpz_fdiv_r(t.n, n.n, d.n); return t; }
void fdivqr(const BigInt& n, const BigInt& d, BigInt& q, BigInt& r) { mpz_fdiv_qr(q.n, r.n, n.n, d.n); }
BigInt cdivq(const BigInt& n, const BigInt& d) { BigInt t; mpz_cdiv_q(t.n, n.n, d.n); return t; }
BigInt cdivr(const BigInt& n, const BigInt& d) { BigInt t; mpz_cdiv_r(t.n, n.n, d.n); return t; }
void cdivqr(const BigInt& n, const BigInt& d, BigInt& q, BigInt& r) { mpz_cdiv_qr(q.n, r.n, n.n, d.n); }
BigInt tdivq(const BigInt& n, const BigInt& d) { BigInt t; mpz_tdiv_q(t.n, n.n, d.n); return t; }
BigInt tdivr(const BigInt& n, const BigInt& d) { BigInt t; mpz_tdiv_r(t.n, n.n, d.n); return t; }
void tdivqr(const BigInt& n, const BigInt& d, BigInt& q, BigInt& r) { mpz_tdiv_qr(q.n, r.n, n.n, d.n); }

BigInt::Random::Random() { gmp_randinit_default(state); gmp_randseed_ui(state, std::random_device{}()); }
BigInt::Random::Random(const Random& r) { gmp_randinit_set(state, r.state); }
BigInt::Random::Random(unsigned long seed) { gmp_randinit_default(state); gmp_randseed_ui(state, seed); }
BigInt::Random::~Random() { gmp_randclear(state); }

BigInt BigInt::Random::operator()(const BigInt& n) { BigInt t; mpz_urandomm(t.n, state, n.n); return t; }
BigInt BigInt::Random::operator()(const BigInt& a, const BigInt& b) { BigInt t; mpz_sub(t.n, b.n, a.n); mpz_urandomm(t.n, state, t.n); mpz_add(t.n, t.n, a.n); return t; }

BigInt::Random default_bigint_random(std::random_device{}());
static time_t prev_reset_time = std::time(nullptr);

static void reset_check(){
    if (std::time(nullptr) - prev_reset_time > 10){
        default_bigint_random = BigInt::Random(std::random_device{}());
        prev_reset_time = std::time(nullptr);
    }
}

BigInt randmod(const BigInt& x) { reset_check(); return default_bigint_random(x); }
BigInt rand(const BigInt& a, const BigInt& b) { reset_check(); return default_bigint_random(a, b); }
BigInt randbits(unsigned long bits) { reset_check(); BigInt t; mpz_urandomb(t.n, default_bigint_random.state, bits); return t; }

const BigFrac BigFrac::zero(0);
const BigFrac BigFrac::one(1);

int BigFrac::sgn() const { return mpq_sgn(n); }
BigFrac BigFrac::abs() const { BigFrac t; mpq_abs(t.n, n); return t; }
BigFrac BigFrac::inv() const { BigFrac t; mpq_inv(t.n, n); return t; }
BigFrac BigFrac::pow(unsigned long exp) const { BigFrac t; mpz_pow_ui(mpq_numref(t.n), mpq_numref(n), exp); mpz_pow_ui(mpq_denref(t.n), mpq_denref(n), exp); return t; }
BigFrac BigFrac::pow(const BigInt& exp, const BigInt& mod) const { BigFrac t; mpz_powm(mpq_numref(t.n), mpq_numref(n), exp.n, mod.n); mpz_powm(mpq_denref(t.n), mpq_denref(n), exp.n, mod.n); return t; }

BigFrac& BigFrac::simplify() { mpq_canonicalize(n); return *this; }

constexpr mpfr_prec_t precision = 100;
constexpr size_t outprec = 20;

const BigFloat BigFloat::zero(0.0);
const BigFloat BigFloat::one(1.0);

const BigFloat BigFloat::pi(PI);
const BigFloat BigFloat::e(E);

BigFloat::BigFloat() { mpfr_init2(n, precision); }
BigFloat::BigFloat(const BigFloat& b) { mpfr_init2(n, mpfr_get_prec(b.n)); mpfr_set(n, b.n, MPFR_RNDD); }
BigFloat::BigFloat(const BigInt& x) { mpfr_init2(n, precision); mpfr_set_z(n, x.n, MPFR_RNDD); }
BigFloat::BigFloat(const BigFrac& x) { mpfr_init2(n, precision); mpfr_set_q(n, x.n, MPFR_RNDD); }
BigFloat::BigFloat(const BigFrac& x, mpfr_prec_t prec) { mpfr_init2(n, prec); mpfr_set_q(n, x.n, MPFR_RNDD); }
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
BigFloat::BigFloat(const float& x, mpfr_prec_t prec) { mpfr_init2(n, prec); mpfr_set_d(n, x, MPFR_RNDD); }
BigFloat::BigFloat(const double& x, mpfr_prec_t prec) { mpfr_init2(n, prec); mpfr_set_d(n, x, MPFR_RNDD); }
BigFloat::BigFloat(const long double& x, mpfr_prec_t prec) { mpfr_init2(n, prec); mpfr_set_d(n, x, MPFR_RNDD); }
BigFloat::~BigFloat() { mpfr_clear(n); }
BigFloat& BigFloat::operator=(const BigFloat& b) { mpfr_set_prec(n, mpfr_get_prec(b.n)); mpfr_set(n, b.n, MPFR_RNDD); return *this; }
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
BigFloat& BigFloat::operator++() { mpfr_add_ui(n, n, 1, MPFR_RNDD); return *this; }
BigFloat& BigFloat::operator--() { mpfr_sub_ui(n, n, 1, MPFR_RNDD); return *this; }
BigFloat BigFloat::operator++(int) { BigFloat t = *this; mpfr_add_ui(n, n, 1, MPFR_RNDD); return t; }
BigFloat BigFloat::operator--(int) { BigFloat t = *this; mpfr_sub_ui(n, n, 1, MPFR_RNDD); return t; }
BigFloat BigFloat::operator-() const { BigFloat t = *this; mpfr_neg(t.n, t.n, MPFR_RNDD); return t; }
BigFloat BigFloat::operator+() const { return *this; }

void BigFloat::set_prec(mpfr_prec_t prec) { mpfr_set_prec(n, prec); }
mpfr_prec_t BigFloat::get_prec() const { return mpfr_get_prec(n); }

void BigFloat::print() const { mpfr_out_str(stdout, 10, outprec, n, MPFR_RNDD); }
void BigFloat::print(int radix) const { mpfr_out_str(stdout, radix, outprec, n, MPFR_RNDD); }
void BigFloat::print(FILE* file) const { mpfr_out_str(file, 10, outprec, n, MPFR_RNDD); }
void BigFloat::print(FILE* file, int radix) const { mpfr_out_str(file, radix, outprec, n, MPFR_RNDD); }
void BigFloat::scan() { mpfr_inp_str(n, stdin, 10, MPFR_RNDD); }
void BigFloat::scan(int radix) { mpfr_inp_str(n, stdin, radix, MPFR_RNDD); }
void BigFloat::scan(FILE* file) { mpfr_inp_str(n, file, 10, MPFR_RNDD); }
void BigFloat::scan(FILE* file, int radix) { mpfr_inp_str(n, file, radix, MPFR_RNDD); }

void BigFloat::print(mpfr_prec_t prec) const { mpfr_out_str(stdout, 10, prec, n, MPFR_RNDD); }
void BigFloat::print(mpfr_prec_t prec, int radix) const { mpfr_out_str(stdout, radix, prec, n, MPFR_RNDD); }
void BigFloat::print(mpfr_prec_t prec, FILE* file) const { mpfr_out_str(file, 10, prec, n, MPFR_RNDD); }
void BigFloat::print(mpfr_prec_t prec, FILE* file, int radix) const { mpfr_out_str(file, radix, prec, n, MPFR_RNDD); }

int BigFloat::sgn() const { return mpfr_sgn(n); }
BigFloat BigFloat::abs() { BigFloat t(*this); mpfr_abs(t.n, t.n, MPFR_RNDD); return t; }
BigFloat BigFloat::pow(unsigned long exp) { BigFloat t; mpfr_pow_ui(t.n, n, exp, MPFR_RNDD); return t; }
BigFloat BigFloat::sqrt() { BigFloat t(*this); mpfr_sqrt(t.n, t.n, MPFR_RNDD); return t; }
BigInt BigFloat::floor() { BigFloat t(*this); mpfr_floor(t.n, t.n); return BigInt(t); }
BigInt BigFloat::ceil() { BigFloat t(*this); mpfr_ceil(t.n, t.n); return BigInt(t); }
BigInt BigFloat::trunc() { BigFloat t(*this); mpfr_trunc(t.n, t.n); return BigInt(t); }

int sgn(const BigFloat& x) { return mpfr_sgn(x.n); }
BigFloat abs(const BigFloat& x) { BigFloat t(x); mpfr_abs(t.n, t.n, MPFR_RNDD); return t; }
BigFloat pow(const BigFloat& x, unsigned long exp) { BigFloat t(x); mpfr_pow_ui(t.n, t.n, exp, MPFR_RNDD); return t; }
BigFloat pow(const BigFloat& x, const BigFloat& exp) { BigFloat t(x); mpfr_pow(t.n, t.n, exp.n, MPFR_RNDD); return t; }
BigFloat sqrt(const BigFloat& x) { BigFloat t(x); mpfr_sqrt(t.n, t.n, MPFR_RNDD); return t; }
BigFloat cbrt(const BigFloat& x) { BigFloat t(x); mpfr_cbrt(t.n, t.n, MPFR_RNDD); return t; }
BigFloat root(const BigFloat& x, unsigned long exp) { BigFloat t(x); mpfr_root(t.n, t.n, exp, MPFR_RNDD); return t; }
BigInt floor(const BigFloat& x) { BigFloat t(x); mpfr_floor(t.n, t.n); return BigInt(t); }
BigInt ceil(const BigFloat& x) { BigFloat t(x); mpfr_ceil(t.n, t.n); return BigInt(t); }
BigInt trunc(const BigFloat& x) { BigFloat t(x); mpfr_trunc(t.n, t.n); return BigInt(t); }

BigFloat get_pi() { BigFloat t; mpfr_const_pi(t.n, MPFR_RNDD); return t; }
BigFloat get_e() { BigFloat t = 1.; mpfr_exp(t.n, t.n, MPFR_RNDD); return t; }
BigFloat get_pi(mpfr_prec_t prec) { BigFloat t; mpfr_set_prec(t.n, prec); mpfr_const_pi(t.n, MPFR_RNDD); return t; }
BigFloat get_e(mpfr_prec_t prec) { BigFloat t(1., prec); mpfr_exp(t.n, t.n, MPFR_RNDD); return t; }

BigFloat exp(const BigFloat& x) { BigFloat t(x); mpfr_exp(t.n, t.n, MPFR_RNDD); return t; }
BigFloat exp2(const BigFloat& x) { BigFloat t(x); mpfr_exp2(t.n, t.n, MPFR_RNDD); return t; }
BigFloat exp10(const BigFloat& x) { BigFloat t(x); mpfr_exp10(t.n, t.n, MPFR_RNDD); return t; }
BigFloat expm1(const BigFloat& x) { BigFloat t(x); mpfr_expm1(t.n, t.n, MPFR_RNDD); return t; }
BigFloat log(const BigFloat& x) { BigFloat t(x); mpfr_log(t.n, t.n, MPFR_RNDD); return t; }
BigFloat log2(const BigFloat& x) { BigFloat t(x); mpfr_log2(t.n, t.n, MPFR_RNDD); return t; }
BigFloat log10(const BigFloat& x) { BigFloat t(x); mpfr_log10(t.n, t.n, MPFR_RNDD); return t; }
BigFloat log1p(const BigFloat& x) { BigFloat t(x); mpfr_log1p(t.n, t.n, MPFR_RNDD); return t; }
BigFloat sin(const BigFloat& x) { BigFloat t(x); mpfr_sin(t.n, t.n, MPFR_RNDD); return t; }
BigFloat cos(const BigFloat& x) { BigFloat t(x); mpfr_cos(t.n, t.n, MPFR_RNDD); return t; }
BigFloat tan(const BigFloat& x) { BigFloat t(x); mpfr_tan(t.n, t.n, MPFR_RNDD); return t; }
BigFloat sinh(const BigFloat& x) { BigFloat t(x); mpfr_sinh(t.n, t.n, MPFR_RNDD); return t; }
BigFloat cosh(const BigFloat& x) { BigFloat t(x); mpfr_cosh(t.n, t.n, MPFR_RNDD); return t; }
BigFloat tanh(const BigFloat& x) { BigFloat t(x); mpfr_tanh(t.n, t.n, MPFR_RNDD); return t; }
BigFloat asin(const BigFloat& x) { BigFloat t(x); mpfr_asin(t.n, t.n, MPFR_RNDD); return t; }
BigFloat acos(const BigFloat& x) { BigFloat t(x); mpfr_acos(t.n, t.n, MPFR_RNDD); return t; }
BigFloat atan(const BigFloat& x) { BigFloat t(x); mpfr_atan(t.n, t.n, MPFR_RNDD); return t; }
BigFloat asinh(const BigFloat& x) { BigFloat t(x); mpfr_asinh(t.n, t.n, MPFR_RNDD); return t; }
BigFloat acosh(const BigFloat& x) { BigFloat t(x); mpfr_acosh(t.n, t.n, MPFR_RNDD); return t; }
BigFloat atanh(const BigFloat& x) { BigFloat t(x); mpfr_atanh(t.n, t.n, MPFR_RNDD); return t; }
BigFloat atan2(const BigFloat& y, const BigFloat& x) { BigFloat t; mpfr_atan2(t.n, y.n, x.n, MPFR_RNDD); return t; }
BigFloat erf(const BigFloat& x) { BigFloat t(x); mpfr_erf(t.n, t.n, MPFR_RNDD); return t; }
BigFloat gamma(const BigFloat& x) { BigFloat t(x); mpfr_gamma(t.n, t.n, MPFR_RNDD); return t; }
BigFloat zeta(const BigFloat& x) { BigFloat t(x); mpfr_zeta(t.n, t.n, MPFR_RNDD); return t; }
