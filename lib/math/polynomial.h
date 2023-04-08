#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "math_base.h"
#include "tensor.h"

// 多项式模板
//TODO 多项式存储结构应根据类型的交换律、结合律、交错律是否存在进行优化
template <typename T, int... L>
using TPolynomial = TTensor<T, L...>;

template <int... L>
using Polynomial = TPolynomial<default_type, L...>;

template <typename T, int... L> constexpr TPolynomial<T, L...> ident(TPolynomial<T, L...>) { return TPolynomial<T, L...>({1}); }

#endif /* POLYNOMIAL_H */
