#ifndef FRACTION_H
#define FRACTION_H

#include "math_base.h"
#include "number_theory.h"

template <typename T>
class TFraction {
public:
    T x;
    T y;

    constexpr TFraction() : x(zero(T())), y(ident(T())) {}
    constexpr TFraction(const TFraction& b) : x(b.x), y(b.y) {}
    constexpr TFraction &operator=(const TFraction& b) { x = b.x; y = b.y; return *this; }
    constexpr TFraction(const T& x) : x(x), y(ident(T())) {}
    constexpr TFraction(const T& x, const T& y) : x(x), y(y) {}
    constexpr TFraction &operator=(const T& x) { this->x = x; this->y = ident(T()); return *this; }
};

template <typename T> constexpr bool is_conjugate_identical(TFraction<T>) { return is_conjugate_identical(T()); }
template <typename T> constexpr bool is_commutative(TFraction<T>) { return is_commutative(T()); }
template <typename T> constexpr bool is_associative(TFraction<T>) { return is_associative(T()); }
template <typename T> constexpr bool is_alternative(TFraction<T>) { return is_alternative(T()); }

template <typename T> constexpr bool is_unital(TFraction<T>) { return is_unital(T()); }
template <typename T> constexpr bool is_dividable(TFraction<T>) { return true; }

template <typename T> constexpr TFraction<T> sim(const TFraction<T>& x) { T t = gcd(abs(x.x), x.y); return TFraction<T>(x.x / t, x.y / t); }

template <typename T> constexpr TFraction<T> operator+(const TFraction<T>& a, const TFraction<T>& b) { return TFraction<T>(a.x * b.y + b.x * a.y, a.y * b.y); }
template <typename T> constexpr TFraction<T> operator+(const TFraction<T>& a, const T& b) { return TFraction<T>(a.x + b * a.y, a.y); }
template <typename T> constexpr TFraction<T> operator-(const TFraction<T>& a, const TFraction<T>& b) { return TFraction<T>(a.x * b.y - b.x * a.y, a.y * b.y); }
template <typename T> constexpr TFraction<T> operator-(const TFraction<T>& a, const T& b) { return TFraction<T>(a.x - b * a.y, a.y); }
template <typename T> constexpr TFraction<T> operator*(const TFraction<T>& a, const TFraction<T>& b) { return TFraction<T>(a.x * b.x, a.y * b.y); }
template <typename T> constexpr TFraction<T> operator*(const TFraction<T>& a, const T& b) { return TFraction<T>(a.x * b, a.y); }
template <typename T> constexpr TFraction<T> operator/(const TFraction<T>& a, const TFraction<T>& b) { return inv(b) * a; }
template <typename T> constexpr TFraction<T> operator/(const TFraction<T>& a, const T& b) { return TFraction<T>(a.x, a.y * b); }
template <typename T> constexpr TFraction<T>& operator+=(TFraction<T>& a, const TFraction<T>& b) { return a = a + b; }
template <typename T> constexpr TFraction<T>& operator+=(TFraction<T>& a, const T& b) { return a = a + b; }
template <typename T> constexpr TFraction<T>& operator-=(TFraction<T>& a, const TFraction<T>& b) { return a = a - b; }
template <typename T> constexpr TFraction<T>& operator-=(TFraction<T>& a, const T& b) { return a = a - b; }
template <typename T> constexpr TFraction<T>& operator*=(TFraction<T>& a, const TFraction<T>& b) { return a = a * b; }
template <typename T> constexpr TFraction<T>& operator*=(TFraction<T>& a, const T& b) { return a = a * b; }
template <typename T> constexpr TFraction<T>& operator/=(TFraction<T>& a, const TFraction<T>& b) { return a = a / b; }
template <typename T> constexpr TFraction<T>& operator/=(TFraction<T>& a, const T& b) { return a = a / b; }
template <typename T> constexpr TFraction<T> operator+(const TFraction<T>& a) { return a; }
template <typename T> constexpr TFraction<T> operator-(const TFraction<T>& a) { return TFraction<T>(-a.x, a.y); }
template <typename T> constexpr TFraction<T> operator^(const TFraction<T>& a, int n) { return pow(a, n); }

template <typename T> constexpr bool operator==(const TFraction<T>& a, const TFraction<T>& b) { return a.x == b.x && a.y == b.y; }
template <typename T> constexpr bool operator==(const TFraction<T>& a, const T& b) { return a.x == b * a.y; }
template <typename T> constexpr bool operator!=(const TFraction<T>& a, const TFraction<T>& b) { return !(a == b); }
template <typename T> constexpr bool operator!=(const TFraction<T>& a, const T& b) { return !(a == b); }
template <typename T> constexpr bool operator<(const TFraction<T>& a, const TFraction<T>& b) { return a.x * b.y < b.x * a.y; }
template <typename T> constexpr bool operator<(const TFraction<T>& a, const T& b) { return a.x < b * a.y; }
template <typename T> constexpr bool operator<(const T& a, const TFraction<T>& b) { return a * b.y < b.x; }
template <typename T> constexpr bool operator>(const TFraction<T>& a, const TFraction<T>& b) { return a.x * b.y > b.x * a.y; }
template <typename T> constexpr bool operator>(const TFraction<T>& a, const T& b) { return a.x > b * a.y; }
template <typename T> constexpr bool operator>(const T& a, const TFraction<T>& b) { return a * b.y > b.x; }
template <typename T> constexpr bool operator<=(const TFraction<T>& a, const TFraction<T>& b) { return a.x * b.y <= b.x * a.y; }
template <typename T> constexpr bool operator<=(const TFraction<T>& a, const T& b) { return a.x <= b * a.y; }
template <typename T> constexpr bool operator<=(const T& a, const TFraction<T>& b) { return a * b.y <= b.x; }
template <typename T> constexpr bool operator>=(const TFraction<T>& a, const TFraction<T>& b) { return a.x * b.y >= b.x * a.y; }
template <typename T> constexpr bool operator>=(const TFraction<T>& a, const T& b) { return a.x >= b * a.y; }
template <typename T> constexpr bool operator>=(const T& a, const TFraction<T>& b) { return a * b.y >= b.x; }

template <typename T> constexpr TFraction<T> ident(TFraction<T>) { return TFraction<T>(ident(T()), ident(T())); }
template <typename T> constexpr TFraction<T> zero(TFraction<T>) { return TFraction<T>(zero(T()), ident(T())); }
template <typename T> constexpr TFraction<T> conj(const TFraction<T>& a) { return TFraction<T>(conj(a.x), conj(a.y)); }
template <typename T> constexpr TFraction<T> inv(const TFraction<T>& a) { return TFraction<T>(a.y, a.x); }
template <typename T> constexpr double norm(const TFraction<T>& a) { return (double)a.x / a.y; }
template <typename T> constexpr double norm2(const TFraction<T>& a) { return (double)(a.x * a.x) / (a.y * a.y); }
template <typename T> constexpr int line(TFraction<T>) { return line(T()); }

template <typename T> constexpr bool isnan(const TFraction<T>& a) { return a.y == zero(T()); }
template <typename T> constexpr TFraction<T> nan(TFraction<T>) { return TFraction<T>(zero(T()), zero(T())); }
template <typename T> constexpr TFraction<T> pow(const TFraction<T>& a, int n) { return TFraction<T>(pow(a.x, n), pow(a.y, n)); }

template <typename T> void print(const TFraction<T>& a, int l) { print(a.x, l); if (a.y != 1){ putchar(l == line(T()) / 2 ? '/' : ' '); print(a.y, l); } }
template <typename T> void print(const TFraction<T>& a) {
    int s = line(a.x);
    for (int i = 0; i < s; i++){
        print(a, i);
        putchar('\n');
    }
}

using Fraction = TFraction<default_int>;

constexpr Fraction operator ""_f(unsigned long long x) { return Fraction(x); }

#endif /* FRACTION_H */