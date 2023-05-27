#ifndef MATRIX_H
#define MATRIX_H

#include "math_base.h"
#include "tensor.h"

#include <functional>

template <typename T, int H, int W>
using TMatrix = TTensor<T, H, W>;

template <int H, int W>
using Matrix = TMatrix<default_type, H, W>;

// matrix, 矩阵
template <typename T, int H, int L, int W> constexpr TMatrix<T, W, H> operator/(const TMatrix<T, L, H>& x, const TMatrix<T, L, W>& y){
    return inv(y) * x;
}

template <typename T, int H, int L, int W> constexpr TMatrix<T, W, H>& operator/=(TMatrix<T, L, H>& x, const TMatrix<T, L, W>& y){
    return x = x / y;
}

// 转置操作, tilde
template <typename T, int H, int W> constexpr TMatrix<T, W, H> operator~(const TMatrix<T, H, W>& x){
    TMatrix<T, W, H> m;
    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++)
            m[j][i] = x[i][j];
    return m;
}

template <typename T, int L> constexpr TMatrix<T, L, L> operator^(const TMatrix<T, L, L>& x, int y){
    return pow(x, y);
}

template <typename T, int L> constexpr TMatrix<T, L, L> ident(TMatrix<T, L, L>) {
    TMatrix<T, L, L> m;
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++)
            m[i][j] = zero(T());
    for (int i = 0; i < L; i++)
        m[i][i] = ident(T());
    return m;
}

template <typename T, int L> constexpr TMatrix<T, L, L> inv(const TMatrix<T, L, L>& x) {
    TMatrix<T, L, L> l, u;
    lu_decom(x, l, u);
    TMatrix<T, L, L> m = zero(TMatrix<T, L, L>());
    for (int i = 0; i < L; i++){
        m[i][i] = ident(T());
        for (int j = i + 1; j < L; j++)
            for (int k = i; k < j; k++)
                m[j][i] -= l[j][k] * m[k][i];
        for (int j = L - 1; j >= 0; j--){
            for (int k = L - 1; k > j; k--)
                m[j][i] -= u[j][k] * m[k][i];
            m[j][i] /= u[j][j];
        }
    }
    return m;
}

// pseudo-inverse, 伪逆矩阵, 逆矩阵不存在时使用
// 以后考虑使用SVD分解求解
template <typename T, int H, int W> constexpr TMatrix<T, W, H> inv(const TMatrix<T, H, W>& x){
    TMatrix<T, W, H> xt = ~x;
    return inv(xt * x) * xt;
}

// input: l, u is zeros
template <typename T, int L> constexpr void lu_decom(const TMatrix<T, L, L>& x, TMatrix<T, L, L>& l, TMatrix<T, L, L>& u) {
    static_assert(is_dividable(T()), "items in matrix must be dividable");
    for (int i = 0; i < L; i++){
        for (int j = 0; j < i; j++){
            l[i][j] = x[i][j];
            u[i][j] = zero(T());
            for (int k = 0; k < j; k++)
                l[i][j] -= l[i][k] * u[k][j];
            l[i][j] /= u[j][j];
        }
        for (int j = i; j < L; j++){
            l[i][j] = zero(T());
            u[i][j] = x[i][j];
            for (int k = 0; k < i; k++)
                u[i][j] -= l[i][k] * u[k][j];
        }
        l[i][i] = ident(T());
    }
}

template <typename T, int H, int W> constexpr void gauss_elim(TMatrix<T, H, W>& x) {
    static_assert(is_dividable(T()), "items in matrix must be dividable");
    bool tag[H] = {false};
    TVector<T, W> row;
    for (int i = 0; i < W, i++){
        for (int j = 0; j < H; j++){
            if (!tag[j] && x[j][i] != zero(T())){
                tag[j] = true;
                T finv = inv(x[j][i]);
                for (int k = i + 1; k < W; k++)
                    row[k] = finv * x[j][k];
                for (int k = 0; k < H; k++){
                    if (!tag[k] && x[k][i] != zero(T())){
                        for (int l = i + 1; l < W; l++)
                            x[k][l] -= x[k][i] * row[l];
                        x[k][i] = zero(T());
                    }
                }
                break;
            }
        }
    }
}

template <typename T, int H, int W> constexpr int rank(const TMatrix<T, H, W>& x) {
    TMatrix<T, H, W> t = x;
    gauss_elim(t);
    int res = 0;
    for (int i = 0; i < H; i++){
        if (t[i][W - 1] != zero(T()))
            res++;
    }
    return res;
}

template <typename T, int L> constexpr T det(const TMatrix<T, L, L>& x) {
    TMatrix<T, L, L> l, u;
    lu_decom(x, l, u);
    T d = ident(T());
    for (int i = 0; i < L; i++)
        d *= u[i][i];
    return d;
}

template <typename T> constexpr T det(const TMatrix<T, 1, 1>& x){
    return x[0][0];
}

template <typename T> constexpr T det(const TMatrix<T, 2, 2>& x){
    return x[0][0] * x[1][1] - x[0][1] * x[1][0];
}

template <typename T> constexpr T det(const TMatrix<T, 3, 3>& x){
    return x[0][0] * (x[1][1] * x[2][2] - x[1][2] * x[2][1]) +
            - x[0][1] * (x[1][0] * x[2][2] - x[1][2] * x[2][0]) +
            + x[0][2] * (x[1][0] * x[2][1] - x[1][1] * x[2][0]);
}

template <typename T, int L>
void calc_eigen(const TMatrix<T, L, L> A, TTensor<T, L>& E, TMatrix<T, L, L>& V) {
    // 将 A 转换为上三角矩阵 U
    TMatrix<T, L, L> LM, U;
    lu_decom(A, LM, U);
    // 从 U 中提取特征值
    for (int i = 0; i < L; i++)
        E[i] = U[i][i];
    // 递归计算每个特征向量
    TTensor<T, L> v;
    for (int i = L - 1; i >= 0; i--) {
        v[i] = 1;
        // 对角线元素为1时直接返回
        if (U[i][i] == 1) {
            V[i] = v;
        } else {
            for (int j = i + 1; j < L; j++) 
                v[j] = U[i][j] / (U[j][j] - U[i][i]);
            // 归一化特征向量
            V[i] = v / norm(v);
            // 更新 U 矩阵
            TMatrix<T, L, L> update = direct(v, v) * U;
            U -= update;
        }
    }
}

template <typename T, int H, int W>
void sort_svd(TMatrix<T, H, H>& U, TMatrix<T, H, W>& S, TMatrix<T, W, W>& V) {
    constexpr int N = (H > W ? W : H);
    for (int i = 0; i < N; i++) {
        int idx = i;
        T mval = S[i][i];
        // 查找最大奇异值所在的列
        for (int j = i + 1; j < N; j++) {
            if (S(j, j) > mval) {
                idx = j;
                mval = S[j][j];
            }
        }
        // 如果不在对角线上，则交换 U, S, V 中相应的列和行
        if (idx != i) {
            T t;
            for (int j = 0; j < N; j++)
                t = S[j][i], S[j][i] = S[j][idx], S[j][idx] = t;
            for (int j = 0; j < N; j++)
                t = U[j][i], U[j][i] = U[j][idx], U[j][idx] = t;
            for (int j = 0; j < N; j++)
                t = V[i][j], V[i][j] = V[idx][j], V[idx][j] = t;
        }
    }
}

template <typename T, int H, int W>
void calc_svd(TMatrix<T, H, W>& A, TMatrix<T, H, H>& U, TMatrix<T, W, W>& V) {
    constexpr int N = (H > W ? W : H);
    TMatrix<T, H, W> S;
    // 计算 A^T * A 的特征值和特征向量
    // 提取 U 和 V 矩阵
    TMatrix<T, W, W> ATA = ~A * A;
    TMatrix<T, H, H> AAT = A * ~A;
    TTensor<T, W> E1;
    TTensor<T, H> E2;
    calc_eigen(ATA, E1, V);
    calc_eigen(AAT, E2, U);
    // 提取 Sigma 矩阵
    for (int i = 0; i < N; i++)
        S[i][i]= sqrt(E1[i]);
    // 对 U 和 V 进行排序
    sort_svd(U, S, V);
}

// 将函数(多项式)应用于矩阵运算上
template <typename T, int L>
TMatrix<T, L, L> apply_mat(const TMatrix<T, L, L>& X, const std::function<T(T)>& f) {
    TTensor<T, L> E;
    TMatrix<T, L, L> V, IV;
    calc_eigen(X, E, V);
    IV = inv(V);
    for (int i = 0; i < L; i++)
        IV[i] *= f(E[i]);
    return V * IV;
}

template <typename T, int L> TMatrix<T, L, L> exp(const TMatrix<T, L, L>& X) { return apply_mat(X, exp); }
template <typename T, int L> TMatrix<T, L, L> log(const TMatrix<T, L, L>& X) { return apply_mat(X, log); }
template <typename T, int L> TMatrix<T, L, L> pow(const TMatrix<T, L, L>& X, T p) { return apply_mat(X, [p](T x) { return pow(x, p); }); }
template <typename T, int L> TMatrix<T, L, L> sqrt(const TMatrix<T, L, L>& X) { return apply_mat(X, sqrt); }
template <typename T, int L> TMatrix<T, L, L> sin(const TMatrix<T, L, L>& X) { return apply_mat(X, sin); }
template <typename T, int L> TMatrix<T, L, L> cos(const TMatrix<T, L, L>& X) { return apply_mat(X, cos); }
template <typename T, int L> TMatrix<T, L, L> tan(const TMatrix<T, L, L>& X) { return apply_mat(X, tan); }
template <typename T, int L> TMatrix<T, L, L> asin(const TMatrix<T, L, L>& X) { return apply_mat(X, asin); }
template <typename T, int L> TMatrix<T, L, L> acos(const TMatrix<T, L, L>& X) { return apply_mat(X, acos); }
template <typename T, int L> TMatrix<T, L, L> atan(const TMatrix<T, L, L>& X) { return apply_mat(X, atan); }
template <typename T, int L> TMatrix<T, L, L> sinh(const TMatrix<T, L, L>& X) { return apply_mat(X, sinh); }
template <typename T, int L> TMatrix<T, L, L> cosh(const TMatrix<T, L, L>& X) { return apply_mat(X, cosh); }
template <typename T, int L> TMatrix<T, L, L> tanh(const TMatrix<T, L, L>& X) { return apply_mat(X, tanh); }
template <typename T, int L> TMatrix<T, L, L> asinh(const TMatrix<T, L, L>& X) { return apply_mat(X, asinh); }
template <typename T, int L> TMatrix<T, L, L> acosh(const TMatrix<T, L, L>& X) { return apply_mat(X, acosh); }
template <typename T, int L> TMatrix<T, L, L> atanh(const TMatrix<T, L, L>& X) { return apply_mat(X, atanh); }

#endif /* MATRIX_H */
