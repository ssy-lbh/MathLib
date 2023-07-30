#ifndef ELLIPTIC_CURVE_H
#define ELLIPTIC_CURVE_H

#include "math_base.h"
#include "number_theory.h"
#include "tensor.h"

template <typename T>
using Point2 = TTensor<T, 2>;

template <typename T>
class TEllipticCurve {
public:
    // y^2 = x^3 + A * x + B
    T a, b;
    TEllipticCurve() : a(zero(T())), b(zero(T())) {}
    TEllipticCurve(T a, T b) : a(a), b(b) {}
    TEllipticCurve(const TEllipticCurve& c) : a(c.a), b(c.b) {}
    TEllipticCurve& operator=(const TEllipticCurve& c) { a = c.a; b = c.b; return *this; }
    bool operator==(const TEllipticCurve& c) const { return a == c.a && b == c.b; }
    bool operator!=(const TEllipticCurve& c) const { return a != c.a || b != c.b; }
    bool is_singular() const { return num(a, 4) * a * a * a + num(b, 27) * b * b == zero(a); }
    bool is_smooth() const { return !is_singular(); }
    bool is_nonsingular() const { return !is_singular(); }
    bool is_supersingular() const { return a == zero(a) && b == zero(b); }
    bool is_ordinary() const { return !is_supersingular(); }
    bool is_special() const { return is_supersingular(); }
    bool is_infinity(const Point2<T>& p) const { return p[0] == zero(p[0]) && p[1] == zero(p[1]); }
    bool is_on_curve(const Point2<T>& p) const { return p[1] * p[1] == (p[0] * p[0] + a) * p[0] + b; }

    Point2<T> add(const Point2<T>& p, const Point2<T>& q) const {
        if (is_infinity(p)) return q;
        if (is_infinity(q)) return p;
        if (p[0] == q[0] && p[1] == -q[1]) return zero(p);
        T s = zero(p[0]);
        if (p == q) s = (num(p[0], 3) * p[0] * p[0] + a) / (num(p[1], 2) * p[1]);
        else s = (q[1] - p[1]) / (q[0] - p[0]);
        T x = s * s - p[0] - q[0];
        T y = s * (p[0] - x) - p[1];
        return Point2<T>({x, y});
    }

    template <typename U>
    Point2<T> mul(U k, const Point2<T>& p) const {
        Point2<T> r = zero(p);
        Point2<T> q = p;
        while (k){
            if (k & 1) r = add(r, q);
            q = add(q, q);
            k >>= 1;
        }
        return r;
    }

    Point2<T> neg(const Point2<T>& p) const {
        return Point2<T>({p[0], -p[1]});
    }

    // 需要数据类型有模数
    Point2<T> random() const {
        T x = rand(a);
        T y_square = (x * x + a) * x + b;
        while (legendre(y_square) == -1){
            x = rand(a);
            y_square = (x * x + a) * x + b;
        }
        T y = sqrt(y_square);
        return Point2<T>({x, y});
    }
};

#endif