#ifndef COMPLEX_H
#define COMPLEX_H

#include "math_base.h"
#include "tensor.h"

template <typename T> class TComplex;

template <typename T>
struct complex_base { using type = T; };
template <typename T>
struct complex_base<TComplex<T>> { using type = typename complex_base<T>::type; };

template <int L, typename T>
struct complex_unfold { using type = T; };
template <int L, typename T>
struct complex_unfold<L, TComplex<T>> { using type = typename std::conditional<(L > 0), typename complex_unfold<L - 1, T>::type, TComplex<T>>::type; };

// 复数以及超复数生成模板
template <typename T>
class TComplex {
public:
    T w;
    T i;

    using base_type = typename complex_base<T>::type;

    constexpr TComplex<T>() : w(zero(w)), i(zero(w)) {}
    constexpr TComplex<T>(const T& x) : w(x), i(zero(w)) {}
    constexpr TComplex<T>(const T& x, const T& y) : w(x), i(y) {}
    constexpr TComplex<T>(const TComplex<T>&) = default;
    constexpr TComplex<T> &operator=(const TComplex<T>&) = default;
    constexpr TComplex<T> &operator=(const T& x) { w = x; i = zero(w); return *this; }
    constexpr base_type& operator[](int i) { return ((base_type*)this)[i]; }
};

template <typename T, int N>
struct nth_complex { using type = TComplex<typename nth_complex<T, N - 1>::type>; };
template <typename T>
struct nth_complex<T, 0> { using type = T; };
template <typename T, int N>
using nth_complex_t = typename nth_complex<T, N>::type;

// 使用前需要static_assert判断N是否为2的幂
template <typename T, int N>
struct lg2_nth_complex { using type = TComplex<typename lg2_nth_complex<T, (N >> 1)>::type>; };
template <typename T>
struct lg2_nth_complex<T, 1> { using type = T; };
template <typename T, int N>
using lg2_nth_complex_t = typename lg2_nth_complex<T, N>::type;

template <typename T>
struct is_complex { static constexpr bool value = false; };
template <typename T>
struct is_complex<TComplex<T>> { static constexpr bool value = true; };
template <typename T>
constexpr bool is_complex_v = is_complex<T>::value;

template <typename T>
struct is_nth_complex { static constexpr int value = 0; };
template <typename T>
struct is_nth_complex<TComplex<T>> { static constexpr int value = is_nth_complex<T>::value + 1; };
template <typename T>
constexpr int is_nth_complex_v = is_nth_complex<T>::value;

template <typename T> constexpr bool is_conjugate_identical(TComplex<T>) { return false; }
template <typename T> constexpr bool is_commutative(TComplex<T>) { return is_conjugate_identical(T()) && is_commutative(T()); }
template <typename T> constexpr bool is_associative(TComplex<T>) { return is_commutative(T()) && is_associative(T()); }
template <typename T> constexpr bool is_alternative(TComplex<T>) { return is_associative(T()) && is_alternative(T()); }

template <typename T> constexpr bool is_unital(TComplex<T>) { return is_unital(T()); }
template <typename T> constexpr bool is_dividable(TComplex<T>) { return true; }

template <typename T1, typename T2>
struct complex_cmp_level { static constexpr int value = 0; }; // 0: T1 == T2, 1: T1 < T2, 2: T2 < T1
template <typename T1, typename T2>
struct complex_cmp_level<TComplex<T1>, TComplex<T2>> : complex_cmp_level<T1, T2> {};
template <typename T1, typename T2>
struct complex_cmp_level<T1, TComplex<T2>> : complex_cmp_level<T1, T2> { static constexpr int value = 1; };
template <typename T1, typename T2>
struct complex_cmp_level<TComplex<T1>, T2> : complex_cmp_level<T1, T2> { static constexpr int value = 2; };

template <typename T> constexpr TComplex<T> operator+(const TComplex<T>& x, const TComplex<T>& y) { return TComplex<T>(x.w + y.w, x.i + y.i); }
template <typename T1, typename T2> constexpr std::enable_if_t<complex_cmp_level<T1, T2>::value != 2, TComplex<decltype(T1() + T2())>>
    operator+(const T1& x, const TComplex<T2>& y) { return TComplex<decltype(T1() + T2())>(x + y.w, y.i); }
template <typename T1, typename T2> constexpr std::enable_if_t<complex_cmp_level<T1, T2>::value != 1, TComplex<decltype(T1() + T2())>>
    operator+(const TComplex<T1>& x, const T2& y) { return TComplex<decltype(T1() + T2())>(x.w + y, x.i); }
template <typename T> constexpr TComplex<T> operator-(const TComplex<T>& x, const TComplex<T>& y) { return TComplex<T>(x.w - y.w, x.i - y.i); }
template <typename T1, typename T2> constexpr std::enable_if_t<complex_cmp_level<T1, T2>::value != 2, TComplex<decltype(T1() - T2())>>
    operator-(const T1& x, const TComplex<T2>& y) { return TComplex<decltype(T1() - T2())>(x - y.w, -y.i); }
template <typename T1, typename T2> constexpr std::enable_if_t<complex_cmp_level<T1, T2>::value != 1, TComplex<decltype(T1() - T2())>>
    operator-(const TComplex<T1>& x, const T2& y) { return TComplex<decltype(T1() - T2())>(x.w - y, x.i); }
// 凯莱-迪克松构造
// 其中体现的特性有
// 1. 若自己的共轭与自己不相等，则高一级无交换律
// 2. 若无交换律，则高一级无结合律
// 3. 若无结合律，则高一级无交错律
// 4. 若无交错律，继续丢掉什么性质建议自己研究
template <typename T> constexpr TComplex<T> operator*(const TComplex<T>& x, const TComplex<T>& y) { return TComplex<T>(x.w * y.w - conj(y.i) * x.i, y.i * x.w + x.i * conj(y.w)); }
template <typename T1, typename T2> constexpr std::enable_if_t<complex_cmp_level<T1, T2>::value != 2, TComplex<decltype(T1() * T2())>>
    operator*(const T1& x, const TComplex<T2>& y) { return TComplex<decltype(T1() * T2())>(x * y.w, x * y.i); }
template <typename T1, typename T2> constexpr std::enable_if_t<complex_cmp_level<T1, T2>::value != 1, TComplex<decltype(T1() * T2())>>
    operator*(const TComplex<T1>& x, const T2& y) { return TComplex<decltype(T1() * T2())>(x.w * y, x.i * y); }
template <typename T> constexpr TComplex<T> operator/(const TComplex<T>& x, const TComplex<T>& y) { return inv(y) * x; }
template <typename T1, typename T2> constexpr TComplex<T1> operator/(const TComplex<T1>& x, const T2& y) { return inv(y) * x; }
template <typename T> constexpr TComplex<T>& operator+=(TComplex<T>& x, const TComplex<T>& y) { x.w += y.w; x.i += y.i; return x; }
template <typename T> constexpr TComplex<T>& operator-=(TComplex<T>& x, const TComplex<T>& y) { x.w -= y.w; x.i -= y.i; return x; }
template <typename T> constexpr TComplex<T>& operator*=(TComplex<T>& x, const TComplex<T>& y) { x = y * x; return x; }
template <typename T1, typename T2> constexpr TComplex<T1>& operator*=(TComplex<T1>& x, const T2& y) { x = y * x; return x; }
template <typename T> constexpr TComplex<T>& operator/=(TComplex<T>& x, const TComplex<T>& y) { return x = x / y; }
template <typename T1, typename T2> constexpr TComplex<T1>& operator/=(TComplex<T1>& x, const T2& y) { return x = x / y; }
template <typename T> constexpr TComplex<T> operator+(const TComplex<T>& x) { return x; }
template <typename T> constexpr TComplex<T> operator-(const TComplex<T>& x) { return TComplex<T>(-x.w, -x.i); }

template <typename T> constexpr bool operator==(const TComplex<T>& x, const TComplex<T>& y) { return x.w == y.w && x.i == y.i; }
template <typename T> constexpr bool operator!=(const TComplex<T>& x, const TComplex<T>& y) { return !(x == y); }

// 请注意
// 1. 以下运算符的优先级较低, 请使用括号
// 2. 八元数及以上不满足结合律, 十六元数及以上不满足交错率
template <typename T> constexpr TComplex<T> operator^(const TComplex<T>& x, int y){ return pow(x, y); }

template <typename T> constexpr TComplex<T> ident(const TComplex<T>& x) { return TComplex<T>(ident(x.w), zero(x.i)); }
template <typename T> constexpr TComplex<T> zero(const TComplex<T>& x) { return TComplex<T>(zero(x.w), zero(x.i)); }
template <typename T, typename U> constexpr TComplex<T> num(const TComplex<T>& x, U n) { if constexpr (std::is_same_v<TComplex<T>, U>) return n; else return TComplex<T>(num(x.w, n), zero(x.w)); }
template <int N, typename T> constexpr TComplex<T> gen(const TComplex<T>& x){ return TComplex<T>(num(x.w, cos(2.0 * PI / N)), num(x.i, sin(2.0 * PI / N))); }
template <typename T> constexpr TComplex<T> gen(const TComplex<T>& x, int n){ return TComplex<T>(num(x.w, cos(2.0 * PI / n)), num(x.i, sin(2.0 * PI / n))); }
template <typename T> constexpr TComplex<T> conj(const TComplex<T>& x) { return TComplex<T>(conj(x.w), -x.i); }
template <typename T> constexpr TComplex<T> inv(const TComplex<T>& x) { return conj(x) / norm2(x); }
template <typename T> constexpr auto norm(const TComplex<T>& x) -> decltype(norm2(x)) { return sqrt(norm2(x)); }
template <typename T> constexpr auto norm2(const TComplex<T>& x) { return norm2(x.w) + norm2(x.i); }
template <typename T> constexpr int line(TComplex<T>) { return line(T()); }

template <typename T> void print(const TComplex<T>& x, int l) {
    putchar('(');
    print(x.w, l);
    putchar(',');
    print(x.i, l);
    putchar(')');
}

using Complex = TComplex<default_type>;
using Quaternion = TComplex<Complex>;
using Octonion = TComplex<Quaternion>;
using Sedenion = TComplex<Octonion>;
using Deduciliion = TComplex<Sedenion>;

template <typename T>
using TQuaternion = TComplex<TComplex<T>>;
template <typename T>
using TOctonion = TComplex<TQuaternion<T>>;
template <typename T>
using TSedenion = TComplex<TOctonion<T>>;
template <typename T>
using TDeducilion = TComplex<TSedenion<T>>;
// 想要更多超复数再加

constexpr Complex operator ""_i(long double x) { return Complex(0.0, (default_type)x); }
constexpr Quaternion operator ""_j(long double x) { return Quaternion({0.0, 0.0}, {(default_type)x, 0.0}); }
constexpr Quaternion operator ""_k(long double x) { return Quaternion({0.0, 0.0}, {0.0, (default_type)x}); }

template <typename T1, typename T2> constexpr T1& get(TComplex<T2>& x, int i) { return ((T1*)&x)[i]; }
template <typename T1, typename T2> constexpr const T1& get(const TComplex<T2>& x, int i) { return ((T1*)&x)[i]; }

template <typename T> constexpr TComplex<T> pow(const TComplex<T>& x, int y){
    TComplex<T> m = ident(TComplex<T>());
    TComplex<T> t = x;
    while (y){
        if (y & 1)
            m *= t;
        t *= t;
        y >>= 1;
    }
    return m;
}

// 以下为常用函数, 目前仅对复数有用, 超复数的话需要以后实现
// 十六元数及以下符合幂结合性(power associativity)，可以使用多项式
template <typename T> constexpr TComplex<T> exp(const TComplex<T>& x) {
    constexpr int nth = is_nth_complex_v<T>;
    if constexpr (nth > 0){
        TComplex<T> t = x;
        auto w = t[0]; t[0] = zero(w);
        auto n = norm(t);
        if (n > zero(n))
            t = ((decltype(n))sin(n) / n) * t;
        t[0] = cos(n) * ident(w);
        return (T)exp(w) * t;
    }
    return (T)exp(x.w) * TComplex<T>(cos(x.i), sin(x.i));
}

template <typename T> constexpr TComplex<T> log(const TComplex<T>& x) {
    constexpr int nth = is_nth_complex_v<T>;
    if constexpr (nth > 0){
        TComplex<T> t = x;
        auto w = t[0]; t[0] = zero(w);
        auto n2 = norm2(t);
        decltype(n2) n = sqrt(n2);
        if (n > zero(n))
            t *= (decltype(n2))atan2(n, w) / n;
        t[0] = (log(norm2(w) + n2) * 0.5) * ident(w);
        return t;
    }
    return TComplex<T>(log(norm2(x)) * 0.5, atan2(x.i, x.w));
}

template <typename T1, typename T2> constexpr TComplex<T1> pow(const TComplex<T1>& x, const T2& y) { return exp(y * log(x)); }
template <typename T> constexpr TComplex<T> sqrt(const TComplex<T>& x) { return exp(0.5 * log(x)); }
template <typename T> constexpr TComplex<T> sin(const TComplex<T>& x) { return TComplex<T>(sin(x.w) * cosh(x.i), cos(x.w) * sinh(x.i)); }
template <typename T> constexpr TComplex<T> cos(const TComplex<T>& x) { return TComplex<T>(cos(x.w) * cosh(x.i), -sin(x.w) * sinh(x.i)); }
template <typename T> constexpr TComplex<T> tan(const TComplex<T>& x) { return sin(x) / cos(x); }
template <typename T> constexpr TComplex<T> sinh(const TComplex<T>& x) { TComplex<T> e = exp(x); return (e - inv(e)) * 0.5; }
template <typename T> constexpr TComplex<T> cosh(const TComplex<T>& x) { TComplex<T> e = exp(x); return (e + inv(e)) * 0.5; }
template <typename T> constexpr TComplex<T> tanh(const TComplex<T>& x) { TComplex<T> e = exp(x), ei = inv(e); return (e - ei) / (e + ei); }

template <typename T> constexpr TComplex<T> asin(const TComplex<T>& x) {
    return TComplex<T>(zero(T()), -ident(T())) * log(TComplex<T>(-x.i, x.w) + sqrt(TComplex<T>(ident(T())) - x * x));
}

//TODO 以下Copilot补全的代码
template <typename T> constexpr TComplex<T> acos(const TComplex<T>& x) {
    return TComplex<T>(zero(T()), -ident(T())) * log(x + sqrt(TComplex<T>(ident(T())) - x * x) * TComplex<T>(zero(T()), ident(T())));
}

template <typename T> constexpr TComplex<T> atan(const TComplex<T>& x) {
    return TComplex<T>(zero(T()), -ident(T()) * 0.5) * log((TComplex<T>(ident(T())) + TComplex<T>(zero(T()), ident(T())) * x) / (TComplex<T>(ident(T())) - TComplex<T>(zero(T()), ident(T())) * x));
}

template <typename T> constexpr TComplex<T> asinh(const TComplex<T>& x) {
    return log(x + sqrt(x * x + TComplex<T>(ident(T()))));
}

template <typename T> constexpr TComplex<T> acosh(const TComplex<T>& x) {
    return log(x + sqrt(x * x - TComplex<T>(ident(T()))));
}

template <typename T> constexpr TComplex<T> atanh(const TComplex<T>& x) {
    return TComplex<T>(zero(T()), ident(T()) * 0.5) * log((TComplex<T>(ident(T())) + x) / (TComplex<T>(ident(T())) - x));
}

template <typename T> constexpr TComplex<T> atan2(const TComplex<T>& y, const TComplex<T>& x) {
    return TComplex<T>(zero(T()), -ident(T()) * 0.5) * log((x + TComplex<T>(zero(T()), ident(T())) * y) / (x - TComplex<T>(zero(T()), ident(T())) * y));
}

#endif /* COMPLEX_H */
