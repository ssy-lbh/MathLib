#include "polynomial.h"

constexpr bool is_power_of_2(uint32_t n){
    return n && ((n & (n - 1)) == 0);
}

void poly_rev(uint32_t poly[], uint32_t n){
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
    poly_rev(poly, n);
    assert((mod - 1) % n == 0);
    if (intt)
        g = pow(g, mod - 2, mod);
    for (uint32_t i = 1; (1u << i) <= n; i++){
        uint32_t h = (1 << (i - 1));
        for (uint32_t j = 0; j < (n >> i); j++){
            uint32_t m = 1; // ��Ԫ
            uint32_t o = pow(g, (mod - 1) >> i, mod); // ����Ԫ����Ԫ
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

void poly_sub(uint32_t poly1[], uint32_t poly2[], uint32_t n){
    for (uint32_t i = 0; i < n; i++)
        poly1[i] -= poly2[i];
}

void poly_add(uint32_t poly1[], uint32_t poly2[], uint32_t n){
    for (uint32_t i = 0; i < n; i++)
        poly1[i] += poly2[i];
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

// ����ʽ�˷�
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

// ����ʽ����
// len(poly) = 2n
// len(tmp) = 2n
// poly <= inv(poly)
// tmp <= ntt(poly)
void poly_inv(uint32_t poly[], uint32_t tmp[], uint32_t n, uint32_t g, uint32_t mod){
    assert(is_power_of_2(n));
    if (n == 1){
        poly[0] = inv(poly[0], mod);
        return;
    }
    poly_inv(poly, tmp, n >> 1, g, mod);
    uint32_t len = n << 1;
    memcpy(tmp, poly, sizeof(uint32_t) * n);
    memset(tmp + n, 0, sizeof(uint32_t) * n);
    poly_ntt(tmp, len, g, mod, false);
    poly_ntt(poly, len, g, mod, false);
    for (uint32_t i = 0; i < len; i++)
        poly[i] = (uint64_t)poly[i] * (mod + 2 - (uint64_t)tmp[i] * poly[i] % mod) % mod;
    poly_ntt(poly, len, g, mod, true);
    memset(poly + n, 0, sizeof(uint32_t) * n);
}

// ����ʽ����
// len(poly1) = 3n
// len(poly2) = 2n
// poly1 <= result
// poly2 <= ntt(poly2)
void poly_div(uint32_t poly1[], uint32_t poly2[], uint32_t n, uint32_t g, uint32_t mod){
    assert(is_power_of_2(n));
    uint32_t len = n << 1;
    poly_inv(poly2, poly1 + n, n, g, mod);
    poly_ntt(poly1, len, g, mod, false);
    poly_ntt(poly2, len, g, mod, false);
    for (uint32_t i = 0; i < len; i++)
        poly1[i] = (uint64_t)poly1[i] * poly2[i] % mod;
    poly_ntt(poly1, len, g, mod, true);
    memset(poly1 + n, 0, sizeof(uint32_t) * n);
}

// ����ʽ����
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

// ����ʽ΢��
// len(poly) = n
// poly <= diff(poly)
void poly_diff(uint32_t poly[], uint32_t n, uint32_t mod){
    for (uint32_t i = 1; i < n; i++)
        poly[i - 1] = (uint64_t)poly[i] * i % mod;
    poly[n - 1] = 0;
}

// ����ʽ����
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

// ����ʽ����
// len(poly) = 2n
// len(tmp_diff) = 2n
// len(tmp_inv) = n
// poly <= log(poly)
void poly_log(uint32_t poly[], uint32_t tmp_diff[], uint32_t tmp_inv[], uint32_t n, uint32_t g, uint32_t mod){
    assert(is_power_of_2(n));
    uint32_t len = n << 1;
    memcpy(tmp_diff, poly, sizeof(uint32_t) * n);
    poly_diff(tmp_diff, n, mod);
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
    uint32_t* tmp_inv = new uint32_t[n];
    poly_log(poly, tmp_diff, tmp_inv, n, g, mod);
    delete[] tmp_inv;
}

void poly_log(uint32_t poly[], uint32_t n, uint32_t g, uint32_t mod){
    uint32_t* tmp_diff = new uint32_t[n << 1];
    uint32_t* tmp_inv = new uint32_t[n];
    poly_log(poly, tmp_diff, tmp_inv, n, g, mod);
    delete[] tmp_diff;
    delete[] tmp_inv;
}

// ����ʽָ��
// len(poly) = 2n
// len(tmp) = 2n
// poly <= exp(poly)
void poly_exp(uint32_t poly[], uint32_t tmp[], uint32_t n, uint32_t g, uint32_t mod){
    assert(is_power_of_2(n));
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
    poly_ntt(poly, len, g, mod, false);
    poly_ntt(tmp, len, g, mod, false);
    for (uint32_t i = 0; i < len; i++)
        poly[i] = (uint64_t)poly[i] * tmp[i] % mod;
    poly_ntt(poly, len, g, mod, true);
    memset(poly + n, 0, sizeof(uint32_t) * n);
}

// ����ʽ��
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