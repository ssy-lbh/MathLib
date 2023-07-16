#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "math_base.h"
#include "tensor.h"
#include "number_theory.h"

// ����ʽģ��
//TODO ����ʽ�洢�ṹӦ�������͵Ľ����ɡ�����ɡ��������Ƿ���ڽ����Ż�
template <typename T, int... L>
using TPolynomial = TTensor<T, L...>;

template <int... L>
using Polynomial = TPolynomial<default_type, L...>;

template <typename T, int... L> constexpr TPolynomial<T, L...> ident(TPolynomial<T, L...>) { return TPolynomial<T, L...>({1}); }

void poly_rev(uint32_t poly[], uint32_t n);
void poly_ntt(uint32_t poly[], uint32_t n, uint32_t g, uint32_t mod, bool intt = false);

void poly_sub(uint32_t poly1[], uint32_t poly2[], uint32_t n);
void poly_add(uint32_t poly1[], uint32_t poly2[], uint32_t n);

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

// 多项式开根
// len(poly) = 2n
// poly <= sqrt(poly)
void poly_sqrt(uint32_t poly[], uint32_t n, uint32_t g, uint32_t mod, uint32_t inv2 = 0);

// 多项式微分
// len(poly) = n
// poly <= diff(poly)
void poly_diff(uint32_t poly[], uint32_t n, uint32_t mod);

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

#endif /* POLYNOMIAL_H */
