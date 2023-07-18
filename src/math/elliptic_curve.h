#ifndef ELLIPTIC_CURVE_H
#define ELLIPTIC_CURVE_H

#include "math_base.h"
#include "number_theory.h"
#include "tensor.h"

template <typename T>
class TEllipticCurve {
public:
    using Point = TTensor<T, 2>;

    // y^2 = x^3 + A * x + B
    T a, b;
    TEllipticCurve() : a(0), b(0) {}
    TEllipticCurve(T a, T b) : a(a), b(b) {}
    TEllipticCurve(const TEllipticCurve& c) : a(c.a), b(c.b) {}
    TEllipticCurve& operator=(const TEllipticCurve& c) { a = c.a; b = c.b; return *this; }
    bool operator==(const TEllipticCurve& c) const { return a == c.a && b == c.b; }
    bool operator!=(const TEllipticCurve& c) const { return a != c.a || b != c.b; }
    bool operator<(const TEllipticCurve& c) const { return a < c.a || (a == c.a && b < c.b); }
    bool operator>(const TEllipticCurve& c) const { return a > c.a || (a == c.a && b > c.b); }
    bool operator<=(const TEllipticCurve& c) const { return a <= c.a || (a == c.a && b <= c.b); }
    bool operator>=(const TEllipticCurve& c) const { return a >= c.a || (a == c.a && b >= c.b); }
    bool is_singular() const { return 4 * a * a * a + 27 * b * b == 0; }
    bool is_smooth() const { return !is_singular(); }
    bool is_nonsingular() const { return !is_singular(); }
    bool is_supersingular() const { return a == 0 && b == 0; }
    bool is_ordinary() const { return !is_supersingular(); }
    bool is_special() const { return is_supersingular(); }
    bool is_infinity(const Point& p) const { return p[0] == 0 && p[1] == 0; }
    bool is_on_curve(const Point& p) const { return p[1] * p[1] == (p[0] * p[0] + a) * p[0] + b; }

    Point add(const Point& p, const Point& q) const {
        if (is_infinity(p)) return q;
        if (is_infinity(q)) return p;
        if (p[0] == q[0] && p[1] == -q[1]) return Point({0, 0});
        T s;
        if (p[0] == q[0] && p[1] == q[1]) s = (3 * p[0] * p[0] + a) / (2 * p[1]);
        else s = (q[1] - p[1]) / (q[0] - p[0]);
        T x = s * s - p[0] - q[0];
        T y = s * (p[0] - x) - p[1];
        return Point({x, y});
    }

    Point mul(const Point& p, uint64_t k) const {
        Point r{0, 0};
        Point q = p;
        while (k){
            if (k & 1) r = add(r, q);
            q = add(q, q);
            k >>= 1;
        }
        return r;
    }

    Point neg(const Point& p) const {
        return Point({p.x, -p.y});
    }

    Point random(T mod) const {
        T x = rand(T());
        T y_square = (x * x + a) * x + b;
        while (legendre(y_square, mod) == -1){
            x = rand(T());
            y_square = (x * x + a) * x + b;
        }
        T y = sqrt(y_square);
        return Point({x, y});
    }
};

#endif