#ifndef MATH_BASE_H
#define MATH_BASE_H

#include <cstdio>
#include <cmath>
#include <cstring>
#include <cassert>
#include <type_traits>

//! 注意事项
//! 1. 本库记得开优化, O1也行, 主要是优化constexpr
//! 2. 无交换律的计算顺序与数学一致
//! 3. 除法定义: x / y => inv(y) * x, 乘等于定义: x *= y => x = y * x
//! 4. 要求性能时把assert替换掉，debug时保留

//! 数域操作定义
// 标准数域定义要求，其他可以扩展

// operator + - * / += -= *= /= =
// T operator+(const T&, const T&);
// T operator*(const T&, const T&);
// T operator/(const T&, const T&);
// T operator-(const T&, const T&);
// T& operator+=(T&, const T&);
// T& operator-=(T&, const T&);
// T& operator*=(T&, const T&);
// T& operator/=(T&, const T&);
// T operator+(const T&);
// T operator-(const T&);

// is_commutative, is_associative, is_alternative

// 幺元和零元(1, 0)(identity, zero)
// 参数仅用于重载, 为了高可扩展性
// T ident(T);
// T zero(T);

// (可选)T gen(T); // 群生成元
// (可选)T gen(T, int n); // n阶群生成元

// 求共轭
// T conj(const T&);

// 求逆元inv(性能考虑,专门的功能)
// T inv(const T&);

// 求范数norm, norm2, 矩阵定义范数为行列式
// T norm(const T&);
// auto norm(const C<T>&) -> decltype(norm(T()));
// T norm2(const T&);
// auto norm2(const C<T>&) -> decltype(norm2(T()));

// int line(T); // 返回打印的行数, 用来控制打印格式
// void print(const T&, int); // 打印, 用来控制打印格式

//! 基本常量与功能

using default_type = float;

constexpr double PI = 3.14159265358979323846;
constexpr double E = 2.71828182845904523536;

template <typename T>
constexpr T sqr(const T& x){
    return x * x;
}

template <typename T1, typename... T>
struct last_type { using type = typename last_type<T...>::type; };
template <typename T1>
struct last_type<T1> { using type = T1; };

template <int CL, int... L>
struct last_int { static constexpr int value = last_int<L...>::value; };
template <int L>
struct last_int<L> { static constexpr int value = L; };

//! 数域定义, 必须留默认构造器
// 其他文件此处定义新类型

//TODO 参照环论写性质判断

// 乘法性质: 共轭相等，交换律，结合律，交错律，幂等律
template <typename T> constexpr std::enable_if_t<std::is_arithmetic_v<T>, bool> is_conjugate_identical(T) { return true; }
template <typename T> constexpr std::enable_if_t<std::is_arithmetic_v<T>, bool> is_commutative(T) { return true; }
template <typename T> constexpr std::enable_if_t<std::is_arithmetic_v<T>, bool> is_associative(T) { return true; }
template <typename T> constexpr std::enable_if_t<std::is_arithmetic_v<T>, bool> is_alternative(T) { return true; }

// 含幺元，可除代数
template <typename T> constexpr std::enable_if_t<std::is_arithmetic_v<T>, bool> is_unital(T) { return true; }
template <typename T> constexpr std::enable_if_t<std::is_floating_point_v<T>, bool> is_dividable(T) { return true; }
template <typename T> constexpr std::enable_if_t<std::is_integral_v<T>, bool> is_dividable(T) { return false; }

//! 数域操作实现如下

// real(float), 实数, 以及整数
template <typename T> constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> ident(T) { return 1; }
template <typename T> constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> zero(T) { return 0; }
template <typename T> constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> conj(T x) { return x; }
template <typename T> constexpr std::enable_if_t<std::is_floating_point_v<T>, T> inv(T x) { return (T)(1.0 / x); }
template <typename T> constexpr std::enable_if_t<std::is_integral_v<T>, T> norm(T x) { return abs(x); }
template <typename T> constexpr std::enable_if_t<std::is_floating_point_v<T>, T> norm(T x) { return fabs(x); }
template <typename T> constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> norm2(T x) { return x * x; }
template <typename T> constexpr std::enable_if_t<std::is_scalar_v<T>, int> line(T) { return 1; }

// 字面量0被重载是有可能的, 建议调用前用构造器转换
inline void print(void* x, int) { printf("%p", x); }

inline void print(char x, int) { printf("%c", x); }
inline void print(wchar_t x, int) { printf("%C", x); }

inline void print(const char* x, int) { printf("%s", x); }
inline void print(const wchar_t* x, int) { printf("%S", x); }

inline void print(bool x, int) { printf("%s", x ? "true" : "false"); }
inline void print(signed char x, int) { printf("%hhd", x); }
inline void print(unsigned char x, int) { printf("%hhu", x); }
inline void print(short x, int) { printf("%hd", x); }
inline void print(unsigned short x, int) { printf("%hu", x); }
inline void print(int x, int) { printf("%d", x); }
inline void print(unsigned int x, int) { printf("%u", x); }
inline void print(long x, int) { printf("%ld", x); }
inline void print(unsigned long x, int) { printf("%lu", x); }
inline void print(long long x, int) { printf("%lld", x); }
inline void print(unsigned long long x, int) { printf("%llu", x); }

inline void print(float x, int) { printf("%g", x); }
inline void print(double x, int) { printf("%lg", x); }
inline void print(long double x, int) { printf("%Lg", x); }

template <typename T> inline std::enable_if_t<std::is_enum_v<T>, void> print(T x, int) { printf("%d", x); }

template <typename T>
void print(const T& x) {
    int s = line(T());
    for (int i = 0; i < s; i++)
        print(x, i);
}

template <typename T>
void println(const T& x) {
    int s = line(T());
    for (int i = 0; i < s; i++){
        print(x, i);
        putchar('\n');
    }
}

template <typename T, typename I> constexpr std::enable_if_t<std::is_integral_v<I>, T> pow(T x, I y){
    T r = ident(T());
    while (y){
        if (y & 1)
            r *= x;
        x *= x;
        y >>= 1;
    }
    return r;
}

int ffsi(int);
int ffsl(long long);
int ctzi(int);
int ctzl(long long);
int clzi(int);
int clzl(long long);
int popcnti(int);
int popcntl(long long);

template <typename T, typename E> constexpr bool equal(const T& x, const T& y, E epsilon) { return norm2(x - y) <= epsilon * epsilon; }

template <typename T1, typename T2> constexpr std::enable_if_t<std::is_scalar_v<T1> && std::is_scalar_v<T2>, void> convert(T1& x, const T2& y) { x = y; }

#endif /* MATH_BASE_H */
