#ifndef BOOLEAN_H
#define BOOLEAN_H

#include "math_base.h"

class Boolean {
public:
    bool n;

    constexpr Boolean() : n(false) {}
    constexpr Boolean(const Boolean& b) : n(b.n) {}
    constexpr Boolean &operator=(const Boolean& b) { n = b.n; return *this; }
    template <typename V> constexpr Boolean(const V& n) : n(n) {}
    template <typename V> constexpr Boolean &operator=(const V& n) { this->n = bool(n); return *this; }
    constexpr operator bool() const { return n; }
};

constexpr bool is_conjugate_identical(Boolean) { return true; }
constexpr bool is_commutative(Boolean) { return true; }
constexpr bool is_associative(Boolean) { return true; }
constexpr bool is_alternative(Boolean) { return true; }

constexpr bool is_unital(Boolean) { return true; }
constexpr bool is_dividable(Boolean) { return true; }

constexpr Boolean operator+(Boolean a, Boolean b) { return a.n || b.n; }
constexpr Boolean operator-(Boolean a, Boolean b) { return a.n && !b.n; }
constexpr Boolean operator*(Boolean a, Boolean b) { return a.n && b.n; }
constexpr Boolean operator/(Boolean a, Boolean b) { return inv(b) * a; }
constexpr Boolean& operator+=(Boolean& a, Boolean b) { return a = a + b; }
constexpr Boolean& operator-=(Boolean& a, Boolean b) { return a = a - b; }
constexpr Boolean& operator*=(Boolean& a, Boolean b) { return a = a * b; }
constexpr Boolean& operator/=(Boolean& a, Boolean b) { return a = a / b; }
constexpr Boolean operator+(Boolean a) { return a; }
constexpr Boolean operator-(Boolean a) { return !a; }
constexpr Boolean operator!(Boolean a) { return !a; }

constexpr bool operator==(Boolean a, Boolean b) { return a.n == b.n; }
constexpr bool operator!=(Boolean a, Boolean b) { return a.n != b.n; }

constexpr Boolean ident(Boolean) { return true; }
constexpr Boolean zero(Boolean) { return false; }
constexpr Boolean conj(Boolean a) { return a; }
constexpr Boolean inv(Boolean a) { assert(a); return a; }
constexpr int norm(Boolean a) { return a.n; }
constexpr int norm2(Boolean a) { return a.n; }
constexpr int line(Boolean) { return 1; }

void print(Boolean a, int) { putchar(a ? 'T' : 'F'); }
void print(Boolean a) { putchar(a ? 'T' : 'F'); }

#endif