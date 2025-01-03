#ifndef MATRIX_H
#define MATRIX_H

#include "math_base.h"
#include "tensor.h"
#include "polynomial.h"

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

// hermitian transpose, 共轭转置操作
template <typename T, int H, int W> constexpr TMatrix<T, W, H> operator~(const TMatrix<T, H, W>& x){
    TMatrix<T, W, H> m;
    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++)
            m[j][i] = conj(x[i][j]);
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

template <typename T, int L, typename U> constexpr TMatrix<T, L, L> num(const TMatrix<T, L, L>& x, U n) {
    if constexpr (std::is_same_v<TMatrix<T, L, L>, U>){
        return n;
    } else {
        TMatrix<T, L, L> m;
        for (int i = 0; i < L; i++)
            for (int j = 0; j < L; j++)
                m[i][j] = zero(T());
        for (int i = 0; i < L; i++)
            m[i][i] = num(x[i][i], n);
        return m;
    }
}

template <typename T, int L> constexpr TMatrix<T, L, L> lu_field_inv(const TMatrix<T, L, L>& x) {
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
            m[j][i] *= inv(u[j][j]);
        }
    }
    return m;
}

template <typename T, int L> constexpr TMatrix<T, L, L> gauss_field_inv(TMatrix<T, L, L> x) {
    TMatrix<T, L, L> ix = ident(TMatrix<T, L, L>());
    for (int i = 0; i < L; i++){
        if (x[i][i] == zero(T())){
            bool found = false;
            for (int j = i + 1; j < L; j++){
                if (x[j][i] != zero(T())){
                    // 交换行
                    TTensor<T, L> t = x[i]; x[i] = x[j]; x[j] = t;
                    t = ix[i]; ix[i] = ix[j]; ix[j] = t;
                    found = true;
                    break;
                }
            }
            if (!found)
                return zero(TMatrix<T, L, L>());
        }
        T finv = inv(x[i][i]);
        x[i][i] = ident(T());
        ix[i] *= finv;
        for (int j = 0; j < L; j++)
            x[i][j] *= finv;
        for (int j = 0; j < L; j++){
            if (j == i)
                continue;
            T f = x[j][i];
            if (f == zero(T()))
                continue;
            x[j] -= x[i] * f;
            ix[j] -= ix[i] * f;
        }
    }
    return ix;
}

template <typename T, int L> constexpr TMatrix<T, L, L> gauss_commutative_ring_inv(TMatrix<T, L, L> x) {
    TMatrix<T, L, L> ix = ident(TMatrix<T, L, L>());
    for (int i = 0; i < L; i++){
        if (isnan(inv(x[i][i]))){
            bool found = false;
            for (int j = i + 1; j < L; j++){
                if (!isnan(inv(x[j][i]))){
                    // 交换行
                    TTensor<T, L> t;
                    t = x[i]; x[i] = x[j]; x[j] = t;
                    t = ix[i]; ix[i] = ix[j]; ix[j] = t;
                    found = true;
                    break;
                }
            }
            if (!found)
                return zero(TMatrix<T, L, L>());
        }
        T finv = inv(x[i][i]);
        x[i][i] = ident(T());
        ix[i] *= finv;
        for (int j = 0; j < L; j++)
            x[i][j] *= finv;
        for (int j = 0; j < L; j++){
            if (j == i)
                continue;
            T f = x[j][i];
            if (f == zero(T()))
                continue;
            x[j] -= x[i] * f;
            ix[j] -= ix[i] * f;
        }
    }
    return ix;
}

// 作用于域上的逆矩阵
template <typename T, int L> constexpr TMatrix<T, L, L> inv(const TMatrix<T, L, L>& x) {
    return gauss_field_inv(x);
}

// pseudo-inverse, 伪逆矩阵, 逆矩阵不存在时使用
// 以后考虑使用SVD分解求解
template <typename T, int H, int W> constexpr TMatrix<T, W, H> inv(const TMatrix<T, H, W>& x){
    TMatrix<T, W, H> xt = ~x;
    return inv(xt * x) * xt;
}

template <typename T, int H, int W> constexpr TMatrix<T, W, H> transpose(const TMatrix<T, H, W>& x){
    TMatrix<T, W, H> m;
    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++)
            m[j][i] = x[i][j];
    return m;
}

// input: l, u is zeros
template <typename T, int L> constexpr void lu_decom(const TMatrix<T, L, L>& x, TMatrix<T, L, L>& l, TMatrix<T, L, L>& u) {
    // static_assert(is_dividable(T()), "items in matrix must be dividable");
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

template <typename T, int L> constexpr void qr_decom(const TMatrix<T, L, L>& x, TMatrix<T, L, L>& q, TMatrix<T, L, L>& r) {
    // static_assert(is_dividable(T()), "items in matrix must be dividable");
    for (int i = 0; i < L; i++){
        q[i] = x[i];
        for (int j = 0; j < i; j++){
            r[j][i] = dot(q[j], x[i]);
            r[i][j] = zero(T());
            q[i] -= r[j][i] * q[j];
        }
        q[i] /= T(norm(q[i]));
        r[i][i] = dot(q[i], x[i]);
    }
}

template <typename T, int H, int W> constexpr int gauss_elim(TMatrix<T, H, W>& x) {
    static_assert(is_dividable(T()), "items in matrix must be dividable");
    int cnt = 0;
    for (int i = 0; i < W; i++){
        if (x[cnt][i] == zero(T())){
            bool found = false;
            for (int j = cnt + 1; j < H; j++){
                if (x[j][i] != zero(T())){
                    x[i] += x[j];
                    found = true;
                }
            }
            if (!found)
                continue;
        }
        T finv = inv(x[cnt][i]);
        x[cnt][i] = ident(T());
        for (int j = i + 1; j < W; j++)
            x[cnt][j] *= finv;
        for (int j = cnt + 1; j < H; j++){
            if (x[j][i] != zero(T())){
                for (int k = i + 1; k < W; k++)
                    x[j][k] -= x[j][i] * x[cnt][k];
                x[j][i] = zero(T());
            }
        }
        cnt++;
    }
    return cnt;
}

template <typename T, int H, int W> constexpr int rank(const TMatrix<T, H, W>& x) {
    TMatrix<T, H, W> t = x;
    return gauss_elim(t);
}

template <typename T, int L> constexpr T trace(const TMatrix<T, L, L>& x) {
    T res = zero(T());
    for (int i = 0; i < L; i++)
        res += x[i][i];
    return res;
}

// tr(AB) = tr(BA)
template <typename T, int N, int L> constexpr T trace(const TMatrix<T, N, L>& x, const TMatrix<T, L, N>& y) {
    return sum(hadamard(x, transpose(y)));
}

template <typename T, int L> constexpr TTensor<T, L> diag(const TMatrix<T, L, L>& x) {
    TTensor<T, L> res;
    for (int i = 0; i < L; i++)
        res[i] = x[i][i];
    return res;
}

template <typename T, int L> constexpr TMatrix<T, L, L> diag(const TTensor<T, L>& x) {
    TMatrix<T, L, L> res = zero(TMatrix<T, L, L>());
    for (int i = 0; i < L; i++)
        res[i][i] = x[i];
    return res;
}

template <typename T, int L> constexpr T det(const TMatrix<T, L, L>& x) {
    TMatrix<T, L, L> l, u;
    lu_decom(x, l, u);
    T res = ident(T());
    for (int i = 0; i < L; i++)
        res *= u[i][i];
    return res;
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

/**
* @brief            Jacobi eigenvalue algorithm
* @param X		    L*L array
* @param V			L*L array
* @param E			L*1 array
* @param precision  precision requirements (float 1e-6, double 1e-10)
* @param max_iter	max_iter number of iterations
* @return
*/
// 求解实对称矩阵的特征值和特征向量
// 无法解决复数矩阵的特征值问题
template <typename T, int L>
constexpr void jacobi_eigen(TMatrix<T, L, L> X, TTensor<T, L>& E, TMatrix<T, L, L>& V, const decltype(norm2(T()))& precision, int max_iter){
    V = ident(X);
    auto precision2 = precision * precision;

	int cnt = 0; // current iteration
	while (cnt < max_iter) {
		// find the largest element on the off-diagonal line of the X
		auto maxv = norm2(X[0][1]);
		int ridx = 0;
		int cidx = 1;
		for (int i = 0; i < L; i++) {			// row
			for (int j = i + 1; j < L; j++) {		// column
				auto d = norm2(X[i][j]);
				if (d > maxv) {
					maxv = d;
					ridx = i;
					cidx = j;
				}
                d = norm2(X[j][i]);
                if (d > maxv) {
                    maxv = d;
                    ridx = j;
                    cidx = i;
                }
			}
		}

		if (maxv < precision2) // precision check 
			break;

		T dbApp = X[ridx][ridx];
		T dbApq = X[ridx][cidx];
		T dbAqq = X[cidx][cidx];
		// compute rotate angle
		T dbAngle = num(dbApp, 0.5) * atan2(num(dbApq, -2) * dbApq, dbAqq - dbApp);
		T dbSinTheta = sin(dbAngle);
		T dbCosTheta = cos(dbAngle);
		T dbSin2Theta = sin(num(dbAngle, 2) * dbAngle);
		T dbCos2Theta = cos(num(dbAngle, 2) * dbAngle);

		X[ridx][ridx] = dbApp * dbCosTheta * dbCosTheta +
			dbAqq * dbSinTheta * dbSinTheta + num(dbApq, 2) * dbApq * dbCosTheta * dbSinTheta;
		X[cidx][cidx] = dbApp * dbSinTheta * dbSinTheta +
			dbAqq * dbCosTheta * dbCosTheta - num(dbApq, 2) * dbApq * dbCosTheta * dbSinTheta;
		X[ridx][cidx] = num(dbAqq, 0.5) * (dbAqq - dbApp) * dbSin2Theta + dbApq * dbCos2Theta;
		X[cidx][ridx] = X[ridx][cidx];

        T tmp;

		for (int i = 0; i < L; i++) {
			if ((i != cidx) && (i != ridx)) {
                TTensor<T, L>& row = X[i];
				tmp = row[ridx];
				row[ridx] = row[cidx] * dbSinTheta + tmp * dbCosTheta;
				row[cidx] = row[cidx] * dbCosTheta - tmp * dbSinTheta;
			}
		}

		for (int j = 0; j < L; j++) {
			if ((j != cidx) && (j != ridx)) {
				tmp = X[ridx][j];
				X[ridx][j] = X[cidx][j] * dbSinTheta + tmp * dbCosTheta;
				X[cidx][j] = X[cidx][j] * dbCosTheta - tmp * dbSinTheta;
			}
		}

		// compute eigenvector
		for (int i = 0; i < L; i++) {
            TTensor<T, L>& row = V[i];
			tmp = row[ridx];
			row[ridx] = row[cidx] * dbSinTheta + tmp * dbCosTheta;
			row[cidx] = row[cidx] * dbCosTheta - tmp * dbSinTheta;
		}

        cnt++;
	}

    for (int i = 0; i < L; i++)
		E[i] = X[i][i];
    V = ~V;
}

// 备选精度和最大迭代次数
// 1. 1e-6, 9 * L * (L - 1)
// 2. 1e-8, 12 * L * (L - 1)
// 3. 1e-10, 15 * L * (L - 1)
template <typename T, int L>
constexpr void jacobi_eigen(const TMatrix<T, L, L>& X, TTensor<T, L>& E, TMatrix<T, L, L>& V){
    jacobi_eigen(X, E, V, 1e-6, 9 * L * (L - 1));
}

// QR分解法求解特征值
// 无法解决复数矩阵的特征值问题
template <typename T, int L>
constexpr void qr_eigen(TMatrix<T, L, L> X, TTensor<T, L>& E, int max_iter){
    for (int i = 0; i < max_iter; i++){
        TMatrix<T, L, L> q, r;
        qr_decom(X, q, r);
        X = r * q;
    }
    E = diag(X);
}

template <typename T, int L>
constexpr void qr_eigen(TMatrix<T, L, L> X, TTensor<T, L>& E, TMatrix<T, L, L>& V, int max_iter){
    V = ident(X);
    for (int i = 0; i < max_iter; i++){
        TMatrix<T, L, L> q, r;
        qr_decom(X, q, r);
        X = r * q;
        V = V * q;
    }
    E = diag(X);
}

template <typename T, int L>
constexpr void calc_eigen(const TMatrix<T, L, L>& A, TTensor<T, L>& E, TMatrix<T, L, L>& V) {
    jacobi_eigen(A, E, V);
}

template <typename T, int H, int W>
constexpr void sort_svd(TMatrix<T, H, H>& U, TMatrix<T, H, W>& S, TMatrix<T, W, W>& V) {
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
constexpr void calc_svd(TMatrix<T, H, W>& A, TMatrix<T, H, H>& U, TMatrix<T, W, W>& V) {
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
constexpr TMatrix<T, L, L> apply_mat(const TMatrix<T, L, L>& X, const std::function<T(T)>& f) {
    TTensor<T, L> E;
    TMatrix<T, L, L> V, IV;
    calc_eigen(X, E, V);
    IV = inv(V);
    for (int i = 0; i < L; i++)
        IV[i] *= f(E[i]);
    return V * IV;
}

template <typename T, int L> constexpr TMatrix<T, L, L> exp(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)exp(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> log(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)log(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> pow(const TMatrix<T, L, L>& X, T p) { return apply_mat(X, std::function<T(T)>([p](T x) { return (T)pow(x, p); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> sqrt(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)sqrt(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> sin(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)sin(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> cos(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)cos(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> tan(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)tan(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> asin(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)asin(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> acos(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)acos(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> atan(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)atan(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> sinh(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)sinh(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> cosh(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)cosh(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> tanh(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)tanh(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> asinh(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)asinh(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> acosh(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)acosh(x); })); }
template <typename T, int L> constexpr TMatrix<T, L, L> atanh(const TMatrix<T, L, L>& X) { return apply_mat(X, std::function<T(T)>([](T x){ return (T)atanh(x); })); }

// @ref https://doi.org/10.1088/0253-6102/49/4/01
template <typename T, int L> constexpr TPolynomial<T, L + 1> eigen_poly(const TMatrix<T, L, L>& X) {
    // A_0 = trace(X^0) = trace(I) = L
    // A_1 = trace(X) = \Sigma \lambda_i
    // A_2 = trace(X^2) = \Sigma \lambda_i^2
    // ...
    // F_0[A(L,0)] = 1
    // F_1[A(L,1)] = A(0,0) * A_1
    // F_2[A(L,2)] = A(1,0) * F_1 * A_1 - A(1,1) * A_2
    // F_3[A(L,3)] = A(2,0) * F_2 * A_1 - A(2,1) * F_1 * A_2 + A(2,2) * A_3, L(L-1)(L-2) = L(L-1)L -2 LL +2 L
    // F_4[A(L,4)] = A(3,0) * F_3 * A_1 - A(3,1) * F_2 * A_2 + A(3,2) * F_1 * A_3 - A(3,3) * A_4, L(L-1)(L-2)(L-3) = L(L-1)(L-2)L -3 L(L-1)L +6 LL -6 L
    // ...
    // 注意 F_i / i!
    TTensor<T, L> A;
    if constexpr (L >= 4){
        // 空间换时间
        // 空间复杂度 O(L^2) -> O(L^2 \sqrt L)
        // 时间复杂度 O(L^4) -> O(L^3 \sqrt L)
        int L2 = (int)ceil(sqrt(L - 0.5));
        TMatrix<T, L, L>* Xs = new TMatrix<T, L, L>[L2];
        TMatrix<T, L, L>* Fs = new TMatrix<T, L, L>[L2];
        Xs[0] = Fs[0] = ident(TMatrix<T, L, L>());
        Xs[1] = X;
        for (int i = 2; i < L2; i++)
            Xs[i] = Xs[i - 1] * X;
        Fs[1] = Xs[L2 - 1] * X;
        for (int i = 2; i < L2; i++)
            Fs[i] = Fs[i - 1] * Fs[1];
        for (int i = 0, rem = 1, idx = 0; i < L; i++, rem++){
            if (rem == L2){
                rem = 0;
                idx++;
            }
            A[i] = trace(Xs[rem], Fs[idx]);
        }
        delete[] Xs;
        delete[] Fs;
    } else {
        TMatrix<T, L, L> x = X;
        for (int i = 0; i < L; i++){
            A[i] = trace(x);
            x = x * X;
        }
    }
    TTensor<T, L> F;
    TPolynomial<T, L + 1> R;
    R[L] = ident(T());
    F[0] = A[0];
    R[L - 1] = -A[0];
    T FT = -ident(T());
    for (int i = 1; i < L; i++){
        F[i] = zero(T());
        T P = ident(T());
        for (int j = 0; j < i; j++){
            F[i] += P * F[i - j - 1] * A[j];
            P *= num(P, j - i);
        }
        F[i] += P * A[i];
        FT *= num(FT, -i - 1);
        R[L - i - 1] = F[i] / FT;
    }
    return R;
}

// 矩阵自嵌套
template <typename T, int L> constexpr TMatrix<T, L * L, L * L> nest(const TMatrix<T, L, L>& x) {
    TMatrix<T, L * L, L * L> m;
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++)
            for (int k = 0; k < L; k++)
                for (int l = 0; l < L; l++)
                    m[i * L + k][j * L + l] = x[i][j] * x[k][l];
    return m;
}

#endif /* MATRIX_H */
