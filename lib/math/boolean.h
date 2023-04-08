#ifndef BOOLEAN_H
#define BOOLEAN_H

#include "math_base.h"

class Boolean {
public:
    bool n;

    constexpr Boolean() : n(false) {}
    constexpr Boolean(const Boolean& b) : n(b.n) {}
    constexpr Boolean &operator=(const Boolean& b) { n = b.n; return *this; }
    constexpr Boolean(bool n = false) : n(n) {}
    constexpr Boolean &operator=(bool n) { this->n = n; return *this; }
    constexpr operator bool() const { return n; }
};

bool is_conjugate_identical(Boolean) { return true; }
bool is_commutative(Boolean) { return true; }
bool is_associative(Boolean) { return true; }
bool is_alternative(Boolean) { return true; }

bool is_unital(Boolean) { return true; }
bool is_dividable(Boolean) { return true; }

constexpr Boolean operator+(Boolean a, Boolean b) { return a.n || b.n; }
constexpr Boolean operator-(Boolean a, Boolean b) { return a.n && !b.n; }
constexpr Boolean operator*(Boolean a, Boolean b) { return a.n && b.n; }
constexpr Boolean operator/(Boolean a, Boolean b) { return inv(b) * a; }
constexpr Boolean& operator+=(Boolean& a, Boolean b) { return a = a + b; }
constexpr Boolean& operator-=(Boolean& a, Boolean b) { return a = a - b; }
constexpr Boolean& operator*=(Boolean& a, Boolean b) { return a = a * b; }
constexpr Boolean& operator/=(Boolean& a, Boolean b) { return a = a / b; }
constexpr Boolean operator+(Boolean a) { return a; }
constexpr Boolean operator!(Boolean a) { return !a; }

constexpr Boolean ident(Boolean) { return true; }
constexpr Boolean zero(Boolean) { return false; }
constexpr Boolean conj(Boolean a) { return a; }
constexpr Boolean inv(Boolean a) { assert(a); return a; }
constexpr bool norm(Boolean a) { return a; }
constexpr bool norm2(Boolean a) { return a; }
constexpr int line(Boolean) { return 1; }

void print(Boolean a, int) { putchar(a ? 'T' : 'F'); }
void print(Boolean a) { putchar(a ? 'T' : 'F'); }

#endif