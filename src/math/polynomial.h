#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "math_base.h"
#include "tensor.h"
#include "number_theory.h"

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
using TPolynomial = TTensor<T, L...>;

template <int... L>
using Polynomial = TPolynomial<default_type, L...>;

template <typename T, int... L> constexpr bool is_conjugate_identical(const TPolynomial<T, L...>&) { return is_conjugate_identical(T()); }
template <typename T, int... L> constexpr bool is_commutative(const TPolynomial<T, L...>&) { return is_commutative(T()); }
template <typename T, int... L> constexpr bool is_associative(const TPolynomial<T, L...>&) { return is_associative(T()); }
template <typename T, int... L> constexpr bool is_alternative(const TPolynomial<T, L...>&) { return is_alternative(T()); }

template <typename T, int... L> constexpr bool is_unital(const TPolynomial<T, L...>&) { return is_unital(T()); }
template <typename T, int... L> constexpr bool is_dividable(const TPolynomial<T, L...>&) { return sizeof...(L) == 1 && is_dividable(T()); }

template <typename T, int... L> constexpr TPolynomial<T, L...> ident(TPolynomial<T, L...>) { return TPolynomial<T, L...>({ident(T())}); }
template <typename T, int... L> constexpr TPolynomial<T, L...> zero(TPolynomial<T, L...>) { return TPolynomial<T, L...>({zero(T())}); }
template <typename T, int... L> constexpr TPolynomial<T, L...> conj(const TPolynomial<T, L...>& a) { return conj(a); }

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

#endif /* POLYNOMIAL_H */
