#include "polynomial.h"

void poly_rev_bit(uint32_t poly[], uint32_t n){
    assert(is_power_of_2(n));
    uint32_t* r = new uint32_t[n];
    r[0] = 0;
    for (uint32_t i = 1; i < n; i++){
        r[i] = (r[i >> 1] >> 1) | ((i & 1) ? (n >> 1) : 0);
        if (i < r[i]){ uint32_t t; t = poly[i]; poly[i] = poly[r[i]]; poly[r[i]] = t; }
    }
    delete[] r;
}

void poly_ntt(uint32_t poly[], uint32_t n, uint32_t g, uint32_t mod, bool intt){
    poly_rev_bit(poly, n);
    assert((mod - 1) % n == 0);
    if (intt)
        g = pow(g, mod - 2, mod);
    for (uint32_t i = 1; (1u << i) <= n; i++){
        uint32_t h = (1 << (i - 1));
        for (uint32_t j = 0; j < (n >> i); j++){
            uint32_t m = 1; // 幺元
            uint32_t o = pow(g, (mod - 1) >> i, mod); // 生成元与逆元
            for (uint32_t k = (j << i); k < (j << i) + h; k++){
                uint32_t A = poly[k], B = (uint64_t)m * poly[k + h] % mod;
                poly[k] = A + B >= mod ? A + B - mod : A + B;
                poly[k + h] = A >= B ? A - B : A - B + mod;
                m = (uint64_t)m * o % mod;
            }
        }
    }
    if (intt){
        uint32_t invn = inv(n, mod);
        for (uint32_t i = 0; i < n; i++)
            poly[i] = (uint64_t)poly[i] * invn % mod;
    }
}

void poly_add(uint32_t poly1[], uint32_t poly2[], uint32_t n){
    for (uint32_t i = 0; i < n; i++)
        poly1[i] += poly2[i];
}

void poly_add(uint32_t poly1[], uint32_t poly2[], uint32_t n, uint32_t mod){
    for (uint32_t i = 0; i < n; i++)
        poly1[i] = (poly1[i] + poly2[i]) % mod;
}

void poly_sub(uint32_t poly1[], uint32_t poly2[], uint32_t n){
    for (uint32_t i = 0; i < n; i++)
        poly1[i] -= poly2[i];
}

void poly_sub(uint32_t poly1[], uint32_t poly2[], uint32_t n, uint32_t mod){
    for (uint32_t i = 0; i < n; i++)
        poly1[i] = (poly1[i] + mod - poly2[i]) % mod;
}

uint32_t poly_deg(uint32_t poly[], uint32_t n){
    for (uint32_t i = n - 1; i > 0; i--)
        if (poly[i])
            return i;
    return 0;
}

void poly_rev(uint32_t poly[], uint32_t n){
    for (uint32_t i = 0, j = n - 1; i < j; i++, j--){
        uint32_t t = poly[i];
        poly[i] = poly[j];
        poly[j] = t;
    }
}

void poly_multiply(uint32_t poly[], uint32_t n, uint32_t k, uint32_t mod){
    for (uint32_t i = 0; i < n; i++)
        poly[i] = (uint64_t)poly[i] * k % mod;
}

void poly_divide(uint32_t poly[], uint32_t n, uint32_t k, uint32_t mod){
    uint32_t invk = inv(k, mod);
    for (uint32_t i = 0; i < n; i++)
        poly[i] = (uint64_t)poly[i] * invk % mod;
}

void poly_fwt_or(uint32_t poly[], uint32_t n, bool inv){
    for (uint32_t i = 1; i < n; i <<= 1)
        for (uint32_t j = 0; j < n; j += i << 1)
            for (uint32_t k = 0; k < i; k++)
                if (inv)
                    poly[i + j + k] = poly[j + k] - poly[i + j + k];
                else
                    poly[i + j + k] = poly[j + k] + poly[i + j + k];
}

void poly_fwt_and(uint32_t poly[], uint32_t n, bool inv){
    for (uint32_t i = 1; i < n; i <<= 1)
        for (uint32_t j = 0; j < n; j += i << 1)
            for (uint32_t k = 0; k < i; k++)
                if (inv)
                    poly[j + k] = poly[i + j + k] - poly[j + k];
                else
                    poly[j + k] = poly[i + j + k] + poly[j + k];
}

void poly_fwt_xor(uint32_t poly[], uint32_t n, bool inv){
    for (uint32_t i = 1; i < n; i <<= 1)
        for (uint32_t j = 0; j < n; j += i << 1)
            for (uint32_t k = 0; k < i; k++){
                uint32_t x = poly[j + k], y = poly[i + j + k];
                if (inv){
                    poly[j + k] = (x + y) >> 1;
                    poly[i + j + k] = (x - y) >> 1;
                }
                else{
                    poly[j + k] = (x + y);
                    poly[i + j + k] = (x - y);
                }
            }
}

// 多项式乘法
// len(poly1) = 2n
// len(poly2) = 2n
// poly1 <= result
// poly2 <= ntt(poly2)
void poly_mul(uint32_t poly1[], uint32_t poly2[], uint32_t n, uint32_t g, uint32_t mod){
    assert(is_power_of_2(n));
    uint32_t len = n << 1;
    poly_ntt(poly1, len, g, mod, false);
    poly_ntt(poly2, len, g, mod, false);
    for (uint32_t i = 0; i < len; i++)
        poly1[i] = (uint64_t)poly1[i] * poly2[i] % mod;
    poly_ntt(poly1, len, g, mod, true);
}

// 多项式求逆
// len(poly) = 2n
// len(tmp) = 2n
// len(tmp2) = 2n
// poly <= inv(poly)
// tmp <= ntt(poly)
void poly_inv(uint32_t poly[], uint32_t tmp[], uint32_t n, uint32_t g, uint32_t mod){
    assert(is_power_of_2(n));
    if (n == 1){
        poly[0] = inv(poly[0], mod);
        return;
    }
    // G = A^(-1) (mod x^[n/2])
    // B = A^(-1) (mod x^n)
    // B - G = 0 (mod x^[n/2])
    // (B - G)^2 = 0 (mod x^n)
    // B^2 - 2BG + G^2 = 0 (mod x^n)
    // AB^2 - 2ABG + AG^2 = 0 (mod x^n)
    // B - 2G + AG^2 = 0 (mod x^n)
    // B = 2G - AG^2 (mod x^n) = G(2 - AG) (mod x^n)
    uint32_t len = n << 1;
    memcpy(tmp, poly, sizeof(uint32_t) * n);
    poly_inv(poly, tmp + n, n >> 1, g, mod);
    memset(tmp + n, 0, sizeof(uint32_t) * n);
    memset(poly + n, 0, sizeof(uint32_t) * n);
    poly_ntt(tmp, len, g, mod, false);
    poly_ntt(poly, len, g, mod, false);
    for (uint32_t i = 0; i < len; i++)
        tmp[i] = (2 + mod - (uint64_t)poly[i] * tmp[i] % mod) * poly[i] % mod;
    poly_ntt(tmp, len, g, mod, true);
    memcpy(poly, tmp, sizeof(uint32_t) * n);
    memset(poly + n, 0, sizeof(uint32_t) * n);
}

// 多项式除法
// len(dividend) = 2n
// len(divisor) = 2n
// len(quotient) = 2n
// dividend <= rev(dividend) * rev(divisor)^(-1) (mod x^n)
// divisor <= ntt(rev(divisor)^_01) (mod x^n)
// quotient <= result
void poly_div(uint32_t dividend[], uint32_t divisor[], uint32_t quotient[], uint32_t n, uint32_t g, uint32_t mod){
    assert(is_power_of_2(n));
    uint32_t len = n << 1;
    uint32_t deg_dividend = poly_deg(dividend, len);
    uint32_t deg_divisor = poly_deg(divisor, len);
    if (deg_dividend < deg_divisor){
        memset(quotient, 0, sizeof(uint32_t) * n);
        return;
    }
    uint32_t deg_quotient = deg_dividend - deg_divisor;
    poly_rev(dividend, deg_dividend);
    poly_rev(divisor, deg_divisor);
    poly_inv(divisor, quotient, n, g, mod);
    memset(quotient + n, 0, sizeof(uint32_t) * n);
    poly_mul(dividend, divisor, n, g, mod);
    memcpy(quotient, dividend, sizeof(uint32_t) * (deg_quotient + 1));
    poly_rev(quotient, deg_quotient + 1);
}

// 多项式开根
// len(poly) = 2n
// poly <= sqrt(poly)
void poly_sqrt(uint32_t poly[], uint32_t n, uint32_t g, uint32_t mod, uint32_t inv2){
    assert(is_power_of_2(n));
    if (n == 1){
        poly[0] = (uint32_t)cipolla(poly[0], mod);
        return;
    }
    if (!inv2)
        inv2 = inv(2, mod);
    poly_sqrt(poly, n >> 1, g, mod, inv2);
    uint32_t len = n << 1;
    uint32_t* tmp = new uint32_t[len];
    memcpy(tmp, poly, sizeof(uint32_t) * n);
    memset(tmp + n, 0, sizeof(uint32_t) * n);
    poly_inv(poly, tmp, n, g, mod);
    poly_ntt(tmp, len, g, mod, false);
    poly_ntt(poly, len, g, mod, false);
    for (uint32_t i = 0; i < len; i++)
        poly[i] = (poly[i] + (uint64_t)tmp[i] * poly[i] % mod) * inv2 % mod;
    poly_ntt(poly, len, g, mod, true);
    memset(poly + n, 0, sizeof(uint32_t) * n);
    delete[] tmp;
}

// 多项式微分
// len(poly) = n
// poly <= derivative(poly)
void poly_deriv(uint32_t poly[], uint32_t n, uint32_t mod){
    for (uint32_t i = 1; i < n; i++)
        poly[i - 1] = (uint64_t)poly[i] * i % mod;
    poly[n - 1] = 0;
}

// 多项式积分
// len(poly) = n
// invs = inv(0, 1, 2, ..., n - 1)
// poly <= int(poly)
void poly_int(uint32_t poly[], uint32_t invs[], uint32_t n, uint32_t mod){
    for (uint32_t i = n - 1; i > 0; i--)
        poly[i] = (uint64_t)poly[i - 1] * invs[i] % mod;
    poly[0] = 0;
}

void poly_int(uint32_t poly[], uint32_t n, uint32_t mod){
    uint32_t* invs = new uint32_t[n];
    inv(invs, n, mod);
    poly_int(poly, invs, n, mod);
    delete[] invs;
}

// 多项式对数
// len(poly) = 2n
// len(tmp_diff) = 2n
// len(tmp_inv) = 2n
// poly <= log(poly)
void poly_log(uint32_t poly[], uint32_t tmp_diff[], uint32_t tmp_inv[], uint32_t n, uint32_t g, uint32_t mod){
    assert(is_power_of_2(n));
    assert(poly[0] == 1);
    uint32_t len = n << 1;
    memcpy(tmp_diff, poly, sizeof(uint32_t) * n);
    poly_deriv(tmp_diff, n, mod);
    memset(tmp_diff + n, 0, sizeof(uint32_t) * n);
    poly_inv(poly, tmp_inv, n, g, mod);
    poly_ntt(poly, len, g, mod, false);
    poly_ntt(tmp_diff, len, g, mod, false);
    for (uint32_t i = 0; i < len; i++)
        poly[i] = (uint64_t)poly[i] * tmp_diff[i] % mod;
    poly_ntt(poly, len, g, mod, true);
    poly_int(poly, n, mod);
    memset(poly + n, 0, sizeof(uint32_t) * n);
}

void poly_log(uint32_t poly[], uint32_t tmp_diff[], uint32_t n, uint32_t g, uint32_t mod){
    uint32_t* tmp_inv = new uint32_t[n << 1];
    poly_log(poly, tmp_diff, tmp_inv, n, g, mod);
    delete[] tmp_inv;
}

void poly_log(uint32_t poly[], uint32_t n, uint32_t g, uint32_t mod){
    uint32_t* tmp_diff = new uint32_t[n << 1];
    uint32_t* tmp_inv = new uint32_t[n << 1];
    poly_log(poly, tmp_diff, tmp_inv, n, g, mod);
    delete[] tmp_diff;
    delete[] tmp_inv;
}

// 多项式指数
// len(poly) = 2n
// len(tmp) = 2n
// poly <= exp(poly)
void poly_exp(uint32_t poly[], uint32_t tmp[], uint32_t n, uint32_t g, uint32_t mod){
    assert(is_power_of_2(n));
    assert(poly[0] == 0);
    if (n == 1){
        poly[0] = 1;
        return;
    }
    poly_exp(poly, tmp, n >> 1, g, mod);
    uint32_t len = n << 1;
    poly_log(poly, tmp, n, g, mod);
    for (uint32_t i = 0; i < n; i++)
        tmp[i] = (mod + poly[i] - tmp[i]) % mod;
    tmp[0] = (tmp[0] + 1) % mod;
    memset(tmp + n, 0, sizeof(uint32_t) * n);
    poly_ntt(poly, len, g, mod, false);
    poly_ntt(tmp, len, g, mod, false);
    for (uint32_t i = 0; i < len; i++)
        poly[i] = (uint64_t)poly[i] * tmp[i] % mod;
    poly_ntt(poly, len, g, mod, true);
    memset(poly + n, 0, sizeof(uint32_t) * n);
}

// 多项式正余弦
// len(poly) = 2n
// len(sinp) = 2n
// len(cosp) = 2n
void poly_sincos(const uint32_t poly[], uint32_t sinp[], uint32_t cosp[], uint32_t n, uint32_t g, uint32_t mod){
    assert(is_power_of_2(n));
    assert(poly[0] == 0);
    uint32_t I = cipolla(-1, mod);
    assert(I != -1);
    uint32_t inv2 = inv(2, mod);
    uint32_t* t1 = new uint32_t[n << 1];
    memcpy(sinp, poly, sizeof(uint32_t) * n);
    poly_multiply(sinp, n, I, mod);
    poly_exp(sinp, cosp, n, g, mod);
    memcpy(t1, sinp, sizeof(uint32_t) * n);
    poly_inv(t1, cosp, n, g, mod);
    memcpy(cosp, sinp, sizeof(uint32_t) * n);
    poly_add(cosp, t1, n, mod);
    poly_multiply(cosp, n, inv2, mod);
    poly_sub(sinp, t1, n, mod);
    poly_multiply(sinp, n, (uint64_t)inv2 * (mod - I) % mod, mod);
    delete[] t1;
}

// 多项式幂
// len(poly) = 2n
// len(tmp) = 2n
// poly <= pow(poly, k)
void poly_pow(uint32_t poly[], uint32_t tmp[], uint32_t n, uint32_t k, uint32_t g, uint32_t mod){
    assert(is_power_of_2(n));
    if (k == 0){
        poly[0] = 1;
        memset(poly + 1, 0, sizeof(uint32_t) * (n - 1));
        return;
    }
    uint32_t len = n << 1;
    poly_log(poly, tmp, n, g, mod);
    for (uint32_t i = 0; i < n; i++)
        tmp[i] = (uint64_t)tmp[i] * k % mod;
    poly_exp(tmp, poly, n, g, mod);
}
