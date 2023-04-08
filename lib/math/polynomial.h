#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "math_base.h"
#include "tensor.h"

// ����ʽģ��
//TODO ����ʽ�洢�ṹӦ�������͵Ľ����ɡ�����ɡ��������Ƿ���ڽ����Ż�
template <typename T, int... L>
using TPolynomial = TTensor<T, L...>;

template <int... L>
using Polynomial = TPolynomial<default_type, L...>;

template <typename T, int... L> constexpr TPolynomial<T, L...> ident(TPolynomial<T, L...>) { return TPolynomial<T, L...>({1}); }

#endif /* POLYNOMIAL_H */
