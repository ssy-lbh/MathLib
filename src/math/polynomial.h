#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "math_base.h"
#include "tensor.h"
#include "number_theory.h"

constexpr bool is_power_of_2(uint32_t n){
    return n && ((n & (n - 1)) == 0);
}

constexpr uint32_t ceil_power_of_2(uint32_t n){
    uint32_t fst = 31 - __builtin_clz(n);
    return (1 << fst) < n ? (1 << (fst + 1)) : (1 << fst);
}

void poly_rev_bit(uint32_t poly[], uint32_t n);
void poly_ntt(uint32_t poly[], uint32_t n, uint32_t g, uint32_t mod, bool intt = false);

void poly_sub(uint32_t poly1[], uint32_t poly2[], uint32_t n);
void poly_add(uint32_t poly1[], uint32_t poly2[], uint32_t n);

uint32_t poly_deg(uint32_t poly[], uint32_t n);
void poly_rev(uint32_t poly[], uint32_t n);

// 快速沃尔什变换
void poly_fwt_or(uint32_t poly[], uint32_t n, bool inv);
void poly_fwt_and(uint32_t poly[], uint32_t n, bool inv);
void poly_fwt_xor(uint32_t poly[], uint32_t n, bool inv);

// 多项式乘法
// len(poly1) = 2n
// len(poly2) = 2n
// poly1 <= result
// poly2 <= ntt(poly2)
void poly_mul(uint32_t poly1[], uint32_t poly2[], uint32_t n, uint32_t g, uint32_t mod);

// 多项式求逆
// len(poly) = 2n
// len(tmp) = 2n
// poly <= inv(poly)
// tmp <= ntt(poly)
void poly_inv(uint32_t poly[], uint32_t tmp[], uint32_t n, uint32_t g, uint32_t mod);

// 多项式除法
// len(poly1) = 3n
// len(poly2) = 2n
// poly1 <= result
// poly2 <= ntt(poly2)
void poly_div(uint32_t poly1[], uint32_t poly2[], uint32_t n, uint32_t g, uint32_t mod);

// 多项式开跟
// len(poly) = 2n
// poly <= sqrt(poly)
void poly_sqrt(uint32_t poly[], uint32_t n, uint32_t g, uint32_t mod, uint32_t inv2 = 0);

// 多项式微分
// len(poly) = n
// poly <= derivative(poly)
void poly_deriv(uint32_t poly[], uint32_t n, uint32_t mod);

// 多项式积分
// len(poly) = n
// invs = inv(0, 1, 2, ..., n - 1)
// poly <= int(poly)
void poly_int(uint32_t poly[], uint32_t invs[], uint32_t n, uint32_t mod);
void poly_int(uint32_t poly[], uint32_t n, uint32_t mod);

// 多项式对数
// len(poly) = 2n
// poly <= log(poly)
void poly_log(uint32_t poly[], uint32_t tmp_diff[], uint32_t tmp_inv[], uint32_t n, uint32_t g, uint32_t mod);
void poly_log(uint32_t poly[], uint32_t tmp_diff[], uint32_t n, uint32_t g, uint32_t mod);
void poly_log(uint32_t poly[], uint32_t n, uint32_t g, uint32_t mod);

// 多项式幂
// len(poly) = 2n
// len(tmp) = 2n
// poly <= pow(poly, k)
void poly_pow(uint32_t poly[], uint32_t tmp[], uint32_t n, uint32_t k, uint32_t g, uint32_t mod);

template <typename T, int... L>
class TPolynomial : public TTensor<T, L...> {};

template <int... L>
using Polynomial = TPolynomial<default_type, L...>;

template <typename T>
struct poly_dim_equal;
template <typename T, int CL, int NL, int... L>
struct poly_dim_equal<TPolynomial<T, CL, NL, L...>> { static constexpr bool value = CL == NL && poly_dim_equal<TPolynomial<T, NL, L...>>::value; };
template <typename T, int CL>
struct poly_dim_equal<TPolynomial<T, CL>> { static constexpr bool value = true; };
template <typename T>
struct poly_dim_equal<TPolynomial<T>> { static constexpr bool value = true; };
template <typename T>
constexpr bool poly_dim_equal_v = poly_dim_equal<T>::value;

template <int L, typename T>
struct poly_dim_add;
template <typename T, int CL, int... L>
struct poly_dim_add<CL, TPolynomial<T, L...>> { using type = TPolynomial<T, CL, L...>; };
template <int L, typename T>
using poly_dim_add_t = typename poly_dim_add<L, T>::type;

template <typename T, int I>
struct poly_dim_get;
template <typename T, int I, int CL, int... L>
struct poly_dim_get<TPolynomial<T, CL, L...>, I> { static constexpr int value = poly_dim_get<TPolynomial<T, L...>, I - 1>::value; };
template <typename T, int CL, int... L>
struct poly_dim_get<TPolynomial<T, CL, L...>, 0> { static constexpr int value = CL; };
template <typename T, int I>
constexpr int poly_dim_get_v = poly_dim_get<T, I>::value;

template <typename T, int... L> constexpr bool is_conjugate_identical(const TPolynomial<T, L...>&) { return is_conjugate_identical(T()); }
template <typename T, int... L> constexpr bool is_commutative(const TPolynomial<T, L...>&) { return is_commutative(T()); }
template <typename T, int... L> constexpr bool is_associative(const TPolynomial<T, L...>&) { return is_associative(T()); }
template <typename T, int... L> constexpr bool is_alternative(const TPolynomial<T, L...>&) { return is_alternative(T()); }

template <typename T, int... L> constexpr bool is_unital(const TPolynomial<T, L...>&) { return is_unital(T()); }
template <typename T, int... L> constexpr bool is_dividable(const TPolynomial<T, L...>&) { return sizeof...(L) == 1 && is_dividable(T()); }

template <typename T, int... L> constexpr TPolynomial<T, L...> ident(TPolynomial<T, L...>) { return TPolynomial<T, L...>({ident(T())}); }
template <typename T, int... L> constexpr TPolynomial<T, L...> zero(TPolynomial<T, L...>) { return TPolynomial<T, L...>({zero(T())}); }
template <typename T, int... L> constexpr TPolynomial<T, L...> conj(const TPolynomial<T, L...>& a) { return conj(a); }

template <typename T, int CL, int... L> constexpr TPolynomial<T, CL, L...> inv(const TPolynomial<T, CL, L...>& a) {
    static_assert(is_power_of_2(CL), "CL must be power of 2");
    TPolynomial<T, CL, L...> b;
    if constexpr (CL == 1){
        b[0] = inv(a[0]);
        return b;
    } else {
        TPolynomial<T, CL, L...> c{};
        b = a;
        TPolynomial<T, CL / 2, L...>& bh = reinterpret_cast<TPolynomial<T, CL / 2, L...>&>(b);
        TPolynomial<T, CL / 2, L...>& ch = reinterpret_cast<TPolynomial<T, CL / 2, L...>&>(c);
        ch = inv(bh);
        return ch * (num(a, 2) - b * ch);
    }
}

template <typename T, int... L> constexpr auto norm(const TPolynomial<T, L...>& a) { return norm(a[0]); }
template <typename T, int... L> constexpr auto norm2(const TPolynomial<T, L...>& a) { return norm2(a[0]); }

template <typename T, int CL> constexpr int line(const TPolynomial<T, CL>&) { return line(T()); }
template <typename T, int CL, int NL, int... L> constexpr int line(const TPolynomial<T, CL, NL, L...>&) { return CL * line(TPolynomial<T, NL, L...>()); }

template <typename T>
void print(const TPolynomial<T>& x, int l){
    print(x.n, l);
}

template <typename T, int CL>
void print(const TPolynomial<T, CL>& x, int l){
    putchar('[');
    if (CL > 0){
        print(x[0], l);
        for (int i = 1; i < CL; i++){
            putchar(' ');
            print(x[i], l);
        }
    }
    putchar(']');
}

template <typename T, int CL, int NL, int... L>
void print(const TPolynomial<T, CL, NL, L...>& x, int l){
    int s = line(TPolynomial<T, NL, L...>());
    int s2 = line(TPolynomial<T, CL, NL, L...>());
    int q = l / s, r = l % s;
    putchar(l == 0 ? '[' : ' ');
    print(x[q], r);
    if (l == s2 - 1)
        putchar(']');
}

template <typename T, int CL, int... L> constexpr TPolynomial<T, CL, L...> operator+(const TPolynomial<T, CL, L...>& a, const TPolynomial<T, CL, L...>& b) {
    TPolynomial<T, CL, L...> c;
    for (int i = 0; i < CL; i++)
        c[i] = a[i] + b[i];
    return c;
}

template <typename T, int CL, int... L> constexpr TPolynomial<T, CL, L...>& operator+=(TPolynomial<T, CL, L...>& a, const TPolynomial<T, CL, L...>& b) {
    for (int i = 0; i < CL; i++)
        a[i] += b[i];
    return a;
}

template <typename T, int CL, int... L> constexpr TPolynomial<T, CL, L...> operator-(const TPolynomial<T, CL, L...>& a, const TPolynomial<T, CL, L...>& b) {
    TPolynomial<T, CL, L...> c;
    for (int i = 0; i < CL; i++)
        c[i] = a[i] - b[i];
    return c;
}

template <typename T, int CL, int... L> constexpr TPolynomial<T, CL, L...>& operator-=(TPolynomial<T, CL, L...>& a, const TPolynomial<T, CL, L...>& b) {
    for (int i = 0; i < CL; i++)
        a[i] -= b[i];
    return a;
}

template <typename T, int CL, int... L> constexpr TPolynomial<T, CL, L...> operator*(const TPolynomial<T, CL, L...>& a, const T& b) {
    TPolynomial<T, CL, L...> c;
    for (int i = 0; i < CL; i++)
        c[i] = a[i] * b;
    return c;
}

template <typename T, int CL, int... L> constexpr TPolynomial<T, CL, L...> operator*(const T& b, const TPolynomial<T, CL, L...>& a) {
    TPolynomial<T, CL, L...> c;
    for (int i = 0; i < CL; i++)
        c[i] = a[i] * b;
    return c;
}

template <typename T, int CL, int... L> constexpr TPolynomial<T, CL, L...>& operator*=(TPolynomial<T, CL, L...>& a, const T& b) {
    for (int i = 0; i < CL; i++)
        a[i] *= b;
    return a;
}

template <typename T, int... L> constexpr TPolynomial<T, L...> operator*(const TPolynomial<T, L...>& a, const TPolynomial<T, L...>& b) {
    TPolynomial<T, L...> at, bt, c;
    at = a; bt = b;
    fft(at); fft(bt);
    c = hadamard(at, bt);
    fft(c, true);
    return c;
}

template <typename T, int... L> constexpr TPolynomial<T, L...>& operator*=(TPolynomial<T, L...>& a, const TPolynomial<T, L...>& b) {
    TPolynomial<T, L...> bt;
    bt = b;
    fft(a); fft(bt);
    a = hadamard(a, bt);
    fft(a, true);
    return a;
}

template <typename T, int... L> constexpr TPolynomial<T, L...> operator/(const TPolynomial<T, L...>& a, const TPolynomial<T, L...>& b) {
    return inv(b) * a;
}

template <typename T, int... L> constexpr TPolynomial<T, L...>& operator/=(TPolynomial<T, L...>& a, const TPolynomial<T, L...>& b) {
    return a = inv(b) * a;
}

template <int N, typename T, int CL, int... L> constexpr TPolynomial<T, CL, L...> deriv(const TPolynomial<T, CL, L...>& a) {
    TPolynomial<T, CL, L...> b;
    if constexpr (N == 0){
        for (int i = 1; i < CL; i++)
            b[i - 1] = a[i] * T(i);
        b[CL - 1] = zero(T());
        return b;
    } else {
        for (int i = 0; i < CL; i++)
            b[i] = deriv<N - 1>(a[i]);
        return b;
    }
}

template <int N, typename T, int CL, int... L> constexpr TPolynomial<T, CL, L...> integ(const TPolynomial<T, CL, L...>& a) {
    TPolynomial<T, CL, L...> b;
    if constexpr (N == 0){
        for (int i = 1; i < CL; i++)
            b[i] = a[i - 1] * inv(T(i));
        b[0] = zero(T());
        return b;
    } else {
        for (int i = 0; i < CL; i++)
            b[i] = integ<N - 1>(a[i]);
        return b;
    }
}

template <int N, typename T, int CL, int... L> constexpr auto subs(const TPolynomial<T, CL, L...>& a, const T& x) {
    if constexpr (N == 0){
        T b = a[CL - 1];
        for (int i = CL - 2; i >= 0; i--)
            b = a[i] + b * x;
        return b;
    } else {
        using I = decltype(subs<N - 1>(a[0], x));
        poly_dim_add_t<CL, I> b;
        for (int i = 0; i < CL; i++)
            b[i] = subs<N - 1>(a[i], x);
        return b;
    }
}

#endif /* POLYNOMIAL_H */
