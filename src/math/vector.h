#ifndef VECTOR_H
#define VECTOR_H

#include "math_base.h"
#include "tensor.h"
#include "complex.h"

// 向量和矩阵模板, 向量默认为列向量
template <typename T, int L>
using TVector = TTensor<T, L, 1>;
template <typename T, int L>
using TRowVector = TTensor<T, 1, L>;

template <int L>
using Vector = TVector<default_type, L>;
template <int L>
using RowVector = TRowVector<default_type, L>;

// vector, 向量

// 叉积推广, 2^n-1维向量的叉积
// 超复数性质: 每个虚数单位平方为-1, 不同复数单位相乘得到另一个复数单位, 结果不计入实部, 八元数及以下乘法表反对称
// 1,3,7维满足性质自身与自身计算结果为0
template <typename T, int L> constexpr TVector<T, L> cross(const TVector<T, L>& x, const TVector<T, L>& y){
    static_assert((L & (L + 1)) == 0, "cross is not defined unless dimension is 2^n-1");
    TVector<T, L> r;
    typename lg2_nth_complex<T, (L + 1)>::type x1, y1;
    x1[0] = 0; y1[0] = 0;
    for (int i = 0; i < L; i++) {
        x1[i+1] = x[i];
        y1[i+1] = y[i];
    }
    auto r1 = x1 * y1;
    for (int i = 0; i < L; i++)
        r[i] = r1[i+1];
    return r;
}

template <typename T> constexpr T cross(const TVector<T, 0>& x, const TVector<T, 0>& y){ return zero(x[0]); }
template <typename T> constexpr T cross(const TVector<T, 1>& x, const TVector<T, 1>& y){ return zero(x[0]); }

template <typename T> constexpr T cross(const TVector<T, 2>& x, const TVector<T, 2>& y){
    return x[0] * y[1] - x[1] * y[0];
}

template <typename T> constexpr TVector<T, 3> cross(const TVector<T, 3>& x, const TVector<T, 3>& y){
    return TVector<T, 3>(x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]);
}

template <typename T, int L> constexpr TVector<T, L> unit(const TVector<T, L>& x){
    return x / norm(x);
}

template <typename T, int L> constexpr TVector<T, L> proj(const TVector<T, L>& x, const TVector<T, L>& t){
    return t * dot(x, t) / norm2(t);
}

template <typename T, int L> constexpr TVector<T, L> proj_unit(const TVector<T, L>& x, const TVector<T, L>& u){
    return u * dot(x, u);
}

#endif /* VECTOR_H */
