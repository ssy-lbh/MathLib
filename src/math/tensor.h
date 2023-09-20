#ifndef TENSOR_H
#define TENSOR_H

#include "math_base.h"

#include <initializer_list>
#include <functional>
#include <cstdint>

// 张量模板, 维度: ..., 行, 列
template <typename T, int... L>
class TTensor;

template <typename T, int CL, int... L>
class TTensor<T, CL, L...> {
public:
    static constexpr int order = sizeof...(L) + 1;
    using base_type = std::conditional_t<order == 1, T, TTensor<T, L...>>;

    TTensor<T, L...> n[CL];

    constexpr TTensor<T, CL, L...>() = default;
    constexpr TTensor<T, CL, L...>(const TTensor<T, CL, L...>&) = default;
    constexpr TTensor<T, CL, L...> &operator=(const TTensor<T, CL, L...>&) = default;
    constexpr base_type& operator[](int i) { return *(base_type*)&n[i]; }
    constexpr const base_type& operator[](int i) const { return *(base_type*)&n[i]; }

    constexpr TTensor<T, CL, L...>(const std::initializer_list<base_type>& l) {
        for (size_t i = 0; i < l.size(); i++)
            n[i] = l.begin()[i];
    }

    template <typename U>
    constexpr TTensor<T, CL, L...> &operator=(const TTensor<U, CL, L...>& t){
        for (int i = 0; i < CL; i++)
            n[i] = t[i];
        return *this;
    }
};

template <typename T>
class TTensor<T> {
public:
    T n;

    constexpr TTensor<T>() = default;
    constexpr TTensor<T>(const TTensor<T>&) = default;
    constexpr TTensor<T> &operator=(const TTensor<T>&) = default;
    constexpr TTensor<T>(const T& t) : n(t) {}
    constexpr TTensor<T> &operator=(const T& t) { n = t; return *this; }
    constexpr operator T() const { return n; }
};

template <int... L>
using Tensor = TTensor<default_type, L...>; 

template <typename T>
struct tensor_dim_equal;
template <typename T, int CL, int NL, int... L>
struct tensor_dim_equal<TTensor<T, CL, NL, L...>> { static constexpr bool value = CL == NL && tensor_dim_equal<TTensor<T, NL, L...>>::value; };
template <typename T, int CL>
struct tensor_dim_equal<TTensor<T, CL>> { static constexpr bool value = true; };
template <typename T>
struct tensor_dim_equal<TTensor<T>> { static constexpr bool value = true; };
template <typename T>
constexpr bool tensor_dim_equal_v = tensor_dim_equal<T>::value;

template <int L, typename T>
struct tensor_dim_add;
template <typename T, int CL, int... L>
struct tensor_dim_add<CL, TTensor<T, L...>> { using type = TTensor<T, CL, L...>; };
template <int L, typename T>
using tensor_dim_add_t = typename tensor_dim_add<L, T>::type;

template <typename T, int I>
struct tensor_dim_get;
template <typename T, int I, int CL, int... L>
struct tensor_dim_get<TTensor<T, CL, L...>, I> { static constexpr int value = tensor_dim_get<TTensor<T, L...>, I - 1>::value; };
template <typename T, int CL, int... L>
struct tensor_dim_get<TTensor<T, CL, L...>, 0> { static constexpr int value = CL; };
template <typename T, int I>
constexpr int tensor_dim_get_v = tensor_dim_get<T, I>::value;

template <typename T, int... L> constexpr bool is_conjugate_identical(TTensor<T, L...>) { return is_conjugate_identical(T()); }
template <typename T, int... L> constexpr bool is_commutative(TTensor<T, L...>) { return false; }
template <typename T, int... L> constexpr bool is_associative(TTensor<T, L...>) { return is_associative(T()); }
template <typename T, int... L> constexpr bool is_alternative(TTensor<T, L...>) { return is_alternative(T()); }

template <typename T, int... L> constexpr bool is_unital(TTensor<T, L...>) { return sizeof...(L) == 2 && tensor_dim_equal_v<TTensor<T, L...>> && is_unital(T()); }
template <typename T, int... L> constexpr bool is_dividable(TTensor<T, L...>) { return sizeof...(L) <= 2 && is_dividable(T()); }

// tensor, 张量
template <typename T> constexpr TTensor<T> operator+(const TTensor<T>& x, const TTensor<T>& y) {
    return TTensor<T>(x.n + y.n);
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> operator+(const TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y){
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = x[i] + y[i];
    return t;
}

template <typename T> constexpr TTensor<T> operator-(const TTensor<T>& x, const TTensor<T>& y) {
    return TTensor<T>(x.n - y.n);
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> operator-(const TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y){
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = x[i] - y[i];
    return t;
}

template <typename T> constexpr TTensor<T> operator*(const T& x, const TTensor<T>& y) {
    return TTensor<T>(x * y.n);
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> operator*(const T& x, const TTensor<T, CL, L...>& y){
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = x * y[i];
    return t;
}

template <typename T> constexpr TTensor<T> operator*(const TTensor<T>& x, const T& y) {
    return TTensor<T>(x.n * y);
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> operator*(const TTensor<T, CL, L...>& x, const T& y){
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = x[i] * y;
    return t;
}

// 张量乘法
//TODO 设计指标轮转的转置功能，用模板参数列表实现
template <typename T, int... L> TTensor<T, L...> operator*(const TTensor<T>& x, const TTensor<T, L...>& y){
    return T(x) * y;
}

template <typename T, int CL, int... L> TTensor<T, L...> operator*(const TTensor<T, CL>& x, const TTensor<T, CL, L...>& y){
    TTensor<T, L...> t = zero(TTensor<T, L...>());
    for (int i = 0; i < CL; i++)
        t += x[i] * y[i];
    return t;
}

template <typename T, int CL, int NL, int... L1, int... L2> auto operator*(const TTensor<T, CL, NL, L1...>& x, const TTensor<T, L2...>& y){
    tensor_dim_add_t<CL, decltype(x[0] * y)> t;
    for (int i = 0; i < CL; i++)
        t[i] = x[i] * y;
    return t;
}

template <typename T, int... L> constexpr TTensor<T, L...> operator/(const TTensor<T, L...>& x, const T& y){
    return inv(y) * x;
}

template <typename T> constexpr TTensor<T>& operator+=(TTensor<T>& x, const T& y){
    x.n += y;
    return x;
}

template <typename T> constexpr TTensor<T>& operator+=(TTensor<T>& x, const TTensor<T>& y){
    x.n += y.n;
    return x;
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...>& operator+=(TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y){
    for (int i = 0; i < CL; i++)
        x[i] += y[i];
    return x;
}

template <typename T> constexpr TTensor<T>& operator-=(TTensor<T>& x, const T& y){
    x.n -= y;
    return x;
}

template <typename T> constexpr TTensor<T>& operator-=(TTensor<T>& x, const TTensor<T>& y){
    x.n -= y.n;
    return x;
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...>& operator-=(TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y){
    for (int i = 0; i < CL; i++)
        x[i] -= y[i];
    return x;
}

template <typename T> constexpr TTensor<T>& operator*=(TTensor<T>& x, const T& y){
    x.n *= y;
    return x;
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...>& operator*=(TTensor<T, CL, L...>& x, const T& y){
    for (int i = 0; i < CL; i++)
        x[i] *= y;
    return x;
}

template <typename T, int... L1, int... L2> constexpr auto& operator*=(TTensor<T, L1...>& x, const TTensor<T, L2...>& y){
    return x = y * x;
}

template <typename T> constexpr TTensor<T>& operator/=(TTensor<T>& x, const T& y){
    x.n /= y;
    return x;
}

template <typename T> constexpr TTensor<T>& operator/=(TTensor<T>& x, const TTensor<T>& y){
    x.n /= y.n;
    return x;
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...>& operator/=(TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y){
    for (int i = 0; i < CL; i++)
        x[i] /= y[i];
    return x;
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...>& operator/=(TTensor<T, CL, L...>& x, const T& y){
    return x *= inv(y);
}

template <typename T, int... L> constexpr TTensor<T, L...> operator+(const TTensor<T, L...>& x){
    return x;
}

template <typename T> constexpr TTensor<T> operator-(const TTensor<T>& x){
    return TTensor<T>(-x.n);
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> operator-(const TTensor<T, CL, L...>& x){
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = -x[i];
    return t;
}

template <typename T> constexpr bool operator!=(const TTensor<T>& x, const TTensor<T>& y){
    return x.n != y.n;
}

template <typename T, int CL, int... L> constexpr bool operator!=(const TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y){
    for (int i = 0; i < CL; i++){
        if (x[i] != y[i])
            return true;
    }
    return false;
}

template <typename T, int CL, int... L> constexpr bool operator==(const TTensor<T>& x, const TTensor<T>& y){
    return !(x != y);
}

template <typename T, int CL, int... L> constexpr bool operator==(const TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y){
    return !(x != y);
}

template <typename T> constexpr TTensor<T> ident(const TTensor<T>& x) {
    return TTensor<T>(ident(x.n));
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> ident(const TTensor<T, CL, L...>& x) {
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = ident(x[i]);
    return t;
}

template <typename T> constexpr TTensor<T> zero(const TTensor<T>& x) {
    return TTensor<T>(zero(x.n));
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> zero(const TTensor<T, CL, L...>& x) {
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = zero(x[i]);
    return t;
}

template <typename T, typename U> constexpr std::enable_if_t<std::is_arithmetic_v<U>, TTensor<T>> num(const TTensor<T>& x, U n) {
    if constexpr (std::is_same_v<TTensor<T>, U>){
        return n;
    } else {
        return TTensor<T>(num<T>(x.n, n));
    }
}

template <typename T, int CL, int... L, typename U> constexpr TTensor<T, CL, L...> num(const TTensor<T, CL, L...>& x, U n) {
    if constexpr (std::is_same_v<TTensor<T, CL, L...>, U>){
        return n;
    } else {
        TTensor<T, CL, L...> t;
        for (int i = 0; i < CL; i++)
            t[i] = num<T>(x[i], n);
        return t;
    }
}

template <int N2, typename T> constexpr TTensor<T> gen(const TTensor<T>& x) {
    return TTensor<T>(gen<N2>(x.n));
}

template <int N2, typename T, int CL, int... L> constexpr TTensor<T, CL, L...> gen(const TTensor<T, CL, L...>& x) {
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = gen<N2>(x[i]);
    return t;
}

template <typename T> constexpr TTensor<T> gen(const TTensor<T>& x, uint32_t n) {
    return TTensor<T>(gen(x.n, n));
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> gen(const TTensor<T, CL, L...>& x, uint32_t n) {
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = gen(x[i], n);
    return t;
}

template <typename T> constexpr TTensor<T> conj(const TTensor<T>& x) {
    return TTensor<T>(conj(x.n));
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> conj(const TTensor<T, CL, L...>& x) {
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = conj(x[i]);
    return t;
}

template <typename T> constexpr TTensor<T> inv(const TTensor<T>& x) {
    return TTensor<T>(inv(x.n));
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> inv(const TTensor<T, CL, L...>& x) {
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = inv(x[i]);
    return t;
}

template <typename T, int CL, int... L> constexpr auto norm(const TTensor<T, CL, L...>& x) { return sqrt(norm2(x)); }

template <typename T> constexpr auto norm2(const TTensor<T>& x) {
    return norm2(x.n);
}

template <typename T, int CL, int... L> constexpr auto norm2(const TTensor<T, CL, L...>& x) {
    decltype(norm2(T())) n;
    n = zero(n);
    for (int i = 0; i < CL; i++)
        n += norm2(x[i]);
    return n;
}

template <typename T> constexpr T dot(const TTensor<T>& x, const TTensor<T>& y) {
    return conj(x.n) * y.n;
}

template <typename T, int CL, int... L> constexpr T dot(const TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y) {
    T d = zero(T());
    for (int i = 0; i < CL; i++)
        d += conj(x[i]) * y[i];
    return d;
}

template <typename T> constexpr T sum(const TTensor<T>& x) {
    return x.n;
}

template <typename T, int L> constexpr T sum(const TTensor<T, L>& x) {
    T s = zero(T());
    for (int i = 0; i < L; i++)
        s += x[i];
    return s;
}

template <typename T, int CL1, int CL2, int... L> constexpr T sum(const TTensor<T, CL1, CL2, L...>& x) {
    T s = zero(T());
    for (int i = 0; i < CL1; i++)
        s += sum(x[i]);
    return s;
}

template <typename T> constexpr T prod(const TTensor<T>& x) {
    return x.n;
}

template <typename T, int L> constexpr T prod(const TTensor<T, L>& x) {
    T s = ident(T());
    for (int i = 0; i < L; i++)
        s *= x[i];
    return s;
}

template <typename T, int CL1, int CL2, int... L> constexpr T prod(const TTensor<T, CL1, CL2, L...>& x) {
    T s = ident(T());
    for (int i = 0; i < CL1; i++)
        s *= prod(x[i]);
    return s;
}

template <typename T, int CL, int... L> constexpr T area(const TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y) {
    return sqrt(dot(x, x) * dot(y, y) - sqr(dot(x, y)));
}

template <typename T> constexpr T hadamard(const T& x, const T& y) {
    return x * y;
}

template <typename T> constexpr T hadamard(const TTensor<T>& x, const TTensor<T>& y) {
    return TTensor<T>(x.n * y.n);
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> hadamard(const TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y) {
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = hadamard(x[i], y[i]);
    return t;
}

// 直积
template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> direct(const T& x, const TTensor<T, CL, L...>& y) {
    return x * y;
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> direct(const TTensor<T>& x, const TTensor<T, CL, L...>& y) {
    return x.n * y;
}

template <typename T, int CL, int... L, int... L2> constexpr TTensor<T, CL, L..., L2...> direct(const TTensor<T, CL, L...>& x, const TTensor<T, L2...>& y) {
    TTensor<T, CL, L..., L2...> t;
    for (int i = 0; i < CL; i++)
        t[i] = direct(x[i], y);
    return t;
}

// 楔积/外积/叉积
template <typename T, int... L1, int... L2> constexpr auto wedge(const TTensor<T, L1...>& x, const TTensor<T, L2...>& y) {
    return direct(x, y) - direct(y, x);
}

template <typename T, int... L1, int... L2> constexpr auto exterior(const TTensor<T, L1...>& x, const TTensor<T, L2...>& y) {
    return direct(x, y) - direct(y, x);
}

template <typename T> constexpr TTensor<T> apply(const TTensor<T>& x, const std::function<T(T)>& f) {
    return TTensor<T>(f(x.n));
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> apply(const TTensor<T, CL, L...>& x, const std::function<T(T)>& f) {
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = apply(x[i], f);
    return t;
}

template <typename T> constexpr TTensor<T> apply(const TTensor<T>& x, const TTensor<T> y, const std::function<T(T, T)>& f) {
    return TTensor<T>(f(x.n, y.n));
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> apply(const TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y, const std::function<T(T, T)>& f) {
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = apply(x[i], y[i], f);
    return t;
}

template <typename T> constexpr int line(TTensor<T>) { return line(T()); }
template <typename T, int CL> constexpr int line(TTensor<T, CL>) { return line(T()); }
template <typename T, int CL, int NL, int... L> constexpr int line(TTensor<T, CL, NL, L...>) { return CL * line(TTensor<T, NL, L...>()); }

template <typename T>
void print(const TTensor<T>& x, int l){
    print(x.n, l);
}

template <typename T, int CL>
void print(const TTensor<T, CL>& x, int l){
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
void print(const TTensor<T, CL, NL, L...>& x, int l){
    int s = line(TTensor<T, NL, L...>());
    int s2 = line(TTensor<T, CL, NL, L...>());
    int q = l / s, r = l % s;
    putchar(l == 0 ? '[' : ' ');
    print(x[q], r);
    if (l == s2 - 1)
        putchar(']');
}

#endif /* TENSOR_H */
