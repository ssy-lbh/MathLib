#include "misc.h"

#include <vector>

BigFrac sqrt_approx1(BigInt x, uint64_t qsize){
    BigInt a = sqrt(x);
    BigInt f0 = a, f1 = 1;
    while (sizeinbase(f1, 2) < qsize){
        // a = a + (x - a2) / (2 * a) = (a2 + x) / (2 * a)
        // p/q = (p^2/q^2 + x) / (2 * p/q) = (p^2 + x * q^2) / (2 * p * q)
        BigInt nf0 = f0 * f0 + x * (f1 * f1);
        BigInt nf1 = 2 * f0 * f1;
        f0 = nf0;
        f1 = nf1;
    }
    return BigFrac(f0, f1);
}

BigFrac sqrt_approx1(BigFrac x, uint64_t qsize){
    BigInt p, q;
    x.get_num_den(p, q);
    return sqrt_approx1(p * q, qsize) / q;
}

BigFrac cbrt_approx1(BigInt x, uint64_t qsize){
    BigInt a = root(x, 3);
    BigInt p = a, q = 1;
    while (sizeinbase(q, 2) < qsize){
        BigInt p2 = p * p, p3 = p2 * p;
        // a = a + (x - a3) / (3 * a2) = (2 * a3 + x) / (3 * a2)
        // p/q = (2 * p^3/q^3 + x) / (3 * p^2/q^2) = (2 * p^3 + x * q^3) / (3 * p^2 * q)
        BigInt np = 2 * p3 + x * (q * q) * q;
        BigInt nq = 3 * p2 * q;
        p = np;
        q = nq;
    }
    return BigFrac(p, q);
}

BigFrac cbrt_approx1(BigFrac x, uint64_t qsize){
    BigInt p, q;
    x.get_num_den(p, q);
    return cbrt_approx1(p * (q * q), qsize) / q;
}

BigFloat pi_borwein3(mpfr_prec_t prec){
    BigFloat a = BigFloat(1., prec) / 3.;
    BigFloat s = (sqrt(BigFloat(3., prec)) - 1) / 2.;
    BigFloat p3 = BigFloat(1., prec);
    mpfr_prec_t iter = log(prec) / log(3.) + 1;
    while (iter --> 0){
        BigFloat r = 3 / (1 + 2 * cbrt(1 - s * s * s));
        s = (r - 1) / 2;
        BigFloat r2 = r * r;
        a = r2 * a - p3 * (r2 - 1);
        p3 *= 3.;
    }
    return 1. / a;
}

BigFloat pi_borwein4(mpfr_prec_t prec){
    BigFloat y = sqrt(BigFloat(2., prec)) - 1;
    BigFloat a = square(y) << 1;
    BigInt p = 8;
    mpfr_prec_t iter = log2(prec) * 0.5 + 1;
    while (iter --> 0){
        BigFloat rt = root(1. - y, 4);
        y = (1. - rt) / (1. + rt);
        a = a * square(square(1. + y)) - p * y * (1. + y + square(y));
        p <<= 2;
    }
    return 1. / a;
}

BigFloat pi_borwein9(mpfr_prec_t prec){
    BigFloat a = BigFloat(1., prec) / 3.;
    BigFloat r = (sqrt(BigFloat(3., prec)) - 1) / 2;
    BigFloat s = cbrt(1 - r * r * r);
    BigFloat p3 = a;
    mpfr_prec_t iter = log(prec) / log(9.) + 1;
    while (iter --> 0){
        BigFloat t = 1 + 2 * r;
        BigFloat u = cbrt(9 * r * (1 + r + r * r));
        BigFloat v = t * t + t * u + u * u;
        BigFloat w = 27 * (1 + s + s * s) / v;
        a = w * a + p3 * (1 - w);
        s = pow(1 - r, 3) / ((t + 2 * u) * v);
        r = cbrt(1 - s * s * s);
        p3 *= 9.;
    }
    return 1. / a;
}

BigFloat pi_agm(mpfr_prec_t prec){
    BigFloat a(1.0, prec), b = sqrt(BigFloat(0.5, prec)), t = BigFloat(0.25, prec);
    BigInt p = 1;
    mpfr_prec_t iter = log2(prec) + 1;
    while (iter --> 0) {
        BigFloat na = (a + b) >> 1;
        b = sqrt(a * b);
        t = t - p * square(a - na);
        a = na;
        p <<= 1;
    }
    return square(a + b) / (4.0 * t);
}

uint32_t pi_bbp16(uint64_t n){
    uint32_t s = 0;
    uint32_t t = 1;
    uint32_t k = 0;
    while (k < n){
        uint32_t r = 8 * k + 1;
        uint32_t a = 4 * k + 1;
        uint32_t b = 2 * k + 1;
        uint32_t c = 1;
        uint32_t d = 0;
        while (r > 0){
            if (r & 1){
                uint32_t t2 = (a * t + b * s) % c;
                s = (a * s + d * t) % c;
                t = t2;
            }
            uint32_t a2 = (a * a + 5 * b * b) % c;
            b = (2 * a * b) % c;
            c = (c * c) % 0x10000;
            d = (d * d + 5 * c) % 0x10000;
            a = a2;
            r >>= 1;
        }
        s = (16 * s) % 0x10000;
        k++;
    }
    return s;
}

// q * q - t * p * p = 1
// @return p / q
BigFrac pell_equation1(BigInt n){
    BigInt m = sqrt(n);
    std::vector<BigInt> f;
    f.push_back(m);
    BigInt a = m, b = 1;
    while (true){
        b = (n - a * a) / b;
        if (b == 1) break;
        BigInt c = (m + a) / b;
        a = c * b - a;
        f.push_back(c);
    }
    BigInt p = 0, q = 1;
    for (int i = f.size() - 1; i >= 0; --i){
        BigInt t = p;
        p = q;
        q = f[i] * q + t;
    }
    return BigFrac(p, q);
}