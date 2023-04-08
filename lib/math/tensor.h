#ifndef TENSOR_H
#define TENSOR_H

#include "math_base.h"

#include <initializer_list>
#include <functional>

// ����ģ��, ά��: ..., ��, ��
template <typename T, int... L>
class TTensor;

template <typename T, int CL, int... L>
class TTensor<T, CL, L...> {
public:
    static constexpr int order = sizeof...(L) + 1;
    using elem_type = std::conditional_t<order == 1, T, TTensor<T, L...>>;

    TTensor<T, L...> n[CL];

    constexpr TTensor<T, CL, L...>() = default;
    constexpr TTensor<T, CL, L...>(const TTensor<T, CL, L...>&) = default;
    constexpr TTensor<T, CL, L...> &operator=(const TTensor<T, CL, L...>&) = default;
    constexpr elem_type& operator[](int i) { return *(elem_type*)&n[i]; }
    constexpr const elem_type& operator[](int i) const { return *(elem_type*)&n[i]; }

    constexpr TTensor<T, CL, L...>(const std::initializer_list<TTensor<T, L...>>& l) {
        memcpy(n, l.begin(),  (l.size() > CL ? CL : l.size()) * sizeof(TTensor<T, L...>));
    }

    // ֱ��һά�������
    constexpr TTensor<T, CL, L...>(const std::initializer_list<T>& l) {
        memcpy(n, l.begin(), (l.size() * sizeof(T)) > sizeof(TTensor<T, CL, L...>) ? sizeof(TTensor<T, CL, L...>) : (l.size() * sizeof(T)));
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

// tensor, ����
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

// �����˷�
//TODO ���ָ����ת��ת�ù��ܣ���ģ������б�ʵ��
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

template <typename T> constexpr TTensor<T> zero(TTensor<T>) {
    return TTensor<T>(zero(T()));
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> zero(TTensor<T, CL, L...>) {
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = zero(TTensor<T, L...>());
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

template <typename T, int CL, int... L> constexpr auto norm(const TTensor<T, CL, L...>& x) { return (T)sqrt(norm2(x)); }

template <typename T> constexpr auto norm2(const TTensor<T>& x) {
    return norm2(x.n);
}

template <typename T, int CL, int... L> constexpr auto norm2(const TTensor<T, CL, L...>& x) {
    decltype(norm2(T())) n = 0;
    for (int i = 0; i < CL; i++)
        n += norm2(x[i]);
    return n;
}

template <typename T> constexpr T dot(const TTensor<T>& x, const TTensor<T>& y) {
    return x.n * y.n;
}

template <typename T, int CL, int... L> constexpr T dot(const TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y) {
    T d = zero(T());
    for (int i = 0; i < CL; i++)
        d += x[i] * y[i];
    return d;
}

template <typename T, int CL, int... L> constexpr T area(const TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y) {
    return sqrt(dot(x, x) * dot(y, y) - sqr(dot(x, y)));
}

template <typename T> constexpr TTensor<T> hadamard(const TTensor<T>& x, const TTensor<T>& y) {
    return TTensor<T>(x.n * y.n);
}

template <typename T, int CL, int... L> constexpr TTensor<T, CL, L...> hadamard(const TTensor<T, CL, L...>& x, const TTensor<T, CL, L...>& y) {
    TTensor<T, CL, L...> t;
    for (int i = 0; i < CL; i++)
        t[i] = hadamard(x[i], y[i]);
    return t;
}

// ֱ��
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

// Ш��/���/���
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
    int q = l / s, r = l % s;
    putchar('[');
    print(x[q], r);
    putchar(']');
}

#endif /* TENSOR_H */
