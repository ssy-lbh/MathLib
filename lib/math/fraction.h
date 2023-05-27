#ifndef FRACTION_H
#define FRACTION_H

#include "math_base.h"
#include "number_theory.h"

template <typename T>
class Fraction {
public:
    T x;
    T y;

    constexpr Fraction() : x(zero(T())), y(ident(T())) {}
    constexpr Fraction(const Fraction& b) : x(b.x), y(b.y) {}
    constexpr Fraction &operator=(const Fraction& b) { x = b.x; y = b.y; return *this; }
    constexpr Fraction(const T& x) : x(x), y(ident(T())) {}
    constexpr Fraction(const T& x, const T& y) : x(x), y(y) {}
    constexpr Fraction &operator=(const T& x) { this->x = x; this->y = ident(T()); return *this; }
};

template <typename T> constexpr bool is_conjugate_identical(Fraction<T>) { return is_conjugate_identical(T()); }
template <typename T> constexpr bool is_commutative(Fraction<T>) { return is_commutative(T()); }
template <typename T> constexpr bool is_associative(Fraction<T>) { return is_associative(T()); }
template <typename T> constexpr bool is_alternative(Fraction<T>) { return is_alternative(T()); }

template <typename T> constexpr bool is_unital(Fraction<T>) { return is_unital(T()); }
template <typename T> constexpr bool is_dividable(Fraction<T>) { return true; }

template <typename T> constexpr Fraction<T> sim(const Fraction<T>& x) { T t = gcd(x.x, x.y); return Fraction<T>(x.x / t, x.y / t); }

template <typename T> constexpr Fraction<T> operator+(const Fraction<T>& a, const Fraction<T>& b) { return sim(Fraction<T>(a.x * b.y + b.x * a.y, a.y * b.y)); }
template <typename T> constexpr Fraction<T> operator-(const Fraction<T>& a, const Fraction<T>& b) { return sim(Fraction<T>(a.x * b.y - b.x * a.y, a.y * b.y)); }
template <typename T> constexpr Fraction<T> operator*(const Fraction<T>& a, const Fraction<T>& b) { return Fraction<T>(a.x * b.x, a.y * b.y); }
template <typename T> constexpr Fraction<T> operator/(const Fraction<T>& a, const Fraction<T>& b) { return inv(b) * a; }
template <typename T> constexpr Fraction<T>& operator+=(Fraction<T>& a, const Fraction<T>& b) { return a = a + b; }
template <typename T> constexpr Fraction<T>& operator-=(Fraction<T>& a, const Fraction<T>& b) { return a = a - b; }
template <typename T> constexpr Fraction<T>& operator*=(Fraction<T>& a, const Fraction<T>& b) { return a = a * b; }
template <typename T> constexpr Fraction<T>& operator/=(Fraction<T>& a, const Fraction<T>& b) { return a = a / b; }
template <typename T> constexpr Fraction<T> operator+(const Fraction<T>& a) { return a; }
template <typename T> constexpr Fraction<T> operator-(const Fraction<T>& a) { return Fraction<T>(-a.x, a.y); }

template <typename T> constexpr bool operator==(const Fraction<T>& a, const Fraction<T>& b) { return a.x == b.x && a.y == b.y; }
template <typename T> constexpr bool operator!=(const Fraction<T>& a, const Fraction<T>& b) { return !(a == b); }

template <typename T> constexpr Fraction<T> ident(Fraction<T>) { return Fraction<T>(ident(T()), ident(T())); }
template <typename T> constexpr Fraction<T> zero(Fraction<T>) { return Fraction<T>(zero(T()), ident(T())); }
template <typename T> constexpr Fraction<T> conj(const Fraction<T>& a) { return Fraction<T>(conj(a.x), conj(a.y)); }
template <typename T> constexpr Fraction<T> inv(const Fraction<T>& a) { return Fraction<T>(a.y, a.x); }
template <typename T> constexpr T norm(const Fraction<T>& a) { return a.x / a.y; }
template <typename T> constexpr T norm2(const Fraction<T>& a) { return (a.x * a.x) / (a.y * a.y); }
template <typename T> constexpr int line(Fraction<T>) { return line(T()); }

template <typename T> void print(const Fraction<T>& a, int l) { print(x, l); putchar(l == line(T()) / 2 ? '/' : ' '); print(y, l); }
template <typename T> void print(const Fraction<T>& a) {
    int s = line(T());
    for (int i = 0; i < s; i++){
        print(a, i);
        putchar('\n');
    }
}

#endif /* FRACTION_H */