#ifndef ELLIPTIC_CURVE_H
#define ELLIPTIC_CURVE_H

#include "math_base.h"
#include "number_theory.h"
#include "tensor.h"

template <typename T>
class TEllipticCurve;

template <typename T>
using TPoint2 = TTensor<T, 2>;

template <typename T>
class TEllipticCurve {
public:
    // y^2 = x^3 + a * x + b
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
    bool is_infinity(const TPoint2<T>& p) const { return p[0] == 0 && p[1] == 0; }
    bool is_on_curve(const TPoint2<T>& p) const { return p[1] * p[1] == (p[0] * p[0] + a) * p[0] + b; }

    TPoint2<T> add(const TPoint2<T>& p, const TPoint2<T>& q) const {
        if (is_infinity(p)) return q;
        if (is_infinity(q)) return p;
        if (p.x == q.x && p.y == -q.y) return TPoint2<T>({0, 0});
        T s;
        if (p.x == q.x && p.y == q.y) s = (3 * p.x * p.x + a) / (2 * p.y);
        else s = (q.y - p.y) / (q.x - p.x);
        T x = s * s - p.x - q.x;
        T y = s * (p.x - x) - p.y;
        return TPoint2<T>({x, y});
    }

    TPoint2<T> mul(const TPoint2<T>& p, uint64_t k) const {
        TPoint2<T> r{0, 0};
        while (k > 0){
            if (k & 1) r = add(r, p);
            p = add(p, p);
            k >>= 1;
        }
    }

    TPoint2<T> neg(const TPoint2<T>& p) const {
        return TPoint2<T>(p.x, -p.y, this);
    }

    TPoint2<T> random() const {
        T x = rand(T());
        T y_square = x * x * x + a * x + b;
        while (legendre(y_square) == -1){
            x = rand(T());
            y_square = x * x * x + a * x + b;
        }
        T y = sqrt(y_square);
        return TPoint2<T>({x, y});
    }
};

using EllipticCurve = TEllipticCurve<default_int>;
using Point2 = TPoint2<default_int>;

void choose_discriminant(uint64_t n, uint64_t& d, uint64_t ms[], uint32_t& num);
bool factor_orders(uint64_t m, uint64_t n, uint64_t& k, uint64_t& q);
bool test_order(EllipticCurve& ec, uint64_t m);
void curve_parameters(uint64_t d, uint64_t n, std::vector<uint64_t>& a, std::vector<uint64_t>& b);

bool atkin_morain(uint64_t n, EllipticCurve& ec){
    /*
    Atkin-Morain ECPP Algorithm.
    Args:
        n: Probable Prime

    Returns:
        Certificate of primality, or False.
    */
   const uint64_t arbitrary_bound = 1e6;

    if (n < arbitrary_bound)
        return miller_rabin(n);

    uint64_t d = 0;
    bool m_found = false;
    uint64_t ms[10];
    uint64_t k, q;
    while (!m_found){
        uint32_t num;
        choose_discriminant(n, d, ms, num);
        for (uint32_t i = 0; i < num; i++){
            uint64_t m = ms[i];
            if (factor_orders(m, n, k, q)){
                std::vector<uint64_t> a, b;
                curve_parameters(d, n, a, b);
                // Test to see if the order of the curve is really m
                ec = EllipticCurve(a.back(), b.back());
                a.pop_back(); b.pop_back();
                while (!test_order(ec, m)){
                    ec = EllipticCurve(a.back(), b.back());
                    a.pop_back(); b.pop_back();
                }
                // #print n, a, b
                m_found = true;
                break;
            }
        }
    }
    // # if no proper m can be found. Go back to choose_discriminant()
    // If this step fails need to return false.

    // operate on point
    Point2 V;
    while (true){
        Point2 P = ec.random();
        Point2 U = ec.mul(P, k); // # U = [m/q]P
        if (!ec.is_infinity(U)) break;
        V = ec.mul(U, q);
    }

    if (!ec.is_infinity(V))
        return false;
    else {
        if (q > arbitrary_bound)
            return atkin_morain(q, ec);
        else {
            if (miller_rabin(q)) return true;
            else return false;
        }
    }
}

#endif