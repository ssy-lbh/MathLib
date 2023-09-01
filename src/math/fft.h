#ifndef FFT_H
#define FFT_H

#include "math_base.h"
#include "tensor.h"

template <typename T, int N>
void fft_rev(TTensor<T, N>& s){
    static_assert(N == (N & -N), "N must be power of 2");
    int* r = new int[N];
    r[0] = 0;
    for (int i = 1; i < N; i++){
        r[i] = (r[i >> 1] >> 1) | ((i & 1) ? (N >> 1) : 0);
        if (i < r[i]){ T t; t = s[i]; s[i] = s[r[i]]; s[r[i]] = t; }
    }
    delete[] r;
}

template <typename T>
void fft_rev(T poly[], int n){
    assert(n == (n & -n) && "n must be power of 2");
    int* r = new int[n];
    r[0] = 0;
    for (int i = 1; i < n; i++){
        r[i] = (r[i >> 1] >> 1) | ((i & 1) ? (n >> 1) : 0);
        if (i < r[i]){ T t; t = poly[i]; poly[i] = poly[r[i]]; poly[r[i]] = t; }
    }
    delete[] r;
}

// 要求T类型是环，满足幂等律，存在单位元和乘法循环子群
// C - {0} = S_1 \times R^+ ==log=> S_1 \times R, [A, B] = 0
// H - {0} = S_3 \times R^+ ==log=> S_2 \times S_1 \times R, [A, B] = ?
template <typename T, int N>
void fft(TTensor<T, N>& poly, bool ifft = false){
    static_assert(N == (N & -N), "N must be power of 2");
    fft_rev(poly);
    for (int i = 1; (1 << i) <= N; i++){
        int h = (1 << (i - 1));
        for (int j = 0; j < (N >> i); j++){
            T o, m = ident(o); // 单位元
            o = ifft ? inv(gen(o, 1 << i)) : gen(o, 1 << i); // 生成元
            for (int k = (j << i); k < (j << i) + h; k++){
                T A = poly[k], B = poly[k + h] * m; // 系数在左，变量在右，右乘生成元的幂
                poly[k] = A + B; poly[k + h] = A - B;
                m *= o;
            }
        }
    }
    if (ifft)
        poly *= inv(T(N));
}

template <typename T>
void fft(T poly[], int n, bool ifft = false){
    assert(n == (n & -n) && "n must be power of 2");
    fft_rev(poly, n);
    for (int i = 1; (1 << i) <= n; i++){
        int h = (1 << (i - 1));
        for (int j = 0; j < (n >> i); j++){
            T o, m = ident(o); // 单位元
            o = ifft ? inv(gen(o, 1 << i)) : gen(o, 1 << i); // 生成元
            for (int k = (j << i); k < (j << i) + h; k++){
                T A = poly[k], B = poly[k + h] * m;
                poly[k] = A + B; poly[k + h] = A - B;
                m *= o;
            }
        }
    }
    if (ifft){
        T t = inv(T(n));
        for (int i = 0; i < n; i++)
            poly[i] *= t;
    }
}

template <typename T, int CL, int... L>
void fft(TTensor<T, CL, L...>& poly, bool ifft = false){
    static_assert(CL == (CL & -CL), "CL must be power of 2");
    for (int i = 0; i < CL; i++)
        fft(poly[i], ifft);
    fft_rev(poly.n, CL);
    for (int i = 1; (1 << i) <= CL; i++){
        int h = (1 << (i - 1));
        for (int j = 0; j < (CL >> i); j++){
            T o, m = ident(o); // 单位元
            o = ifft ? inv(gen(o, 1 << i)) : gen(o, 1 << i); // 生成元
            for (int k = (j << i); k < (j << i) + h; k++){
                TTensor<T, L...> A = poly[k], B = poly[k + h] * m; // 系数在左，变量在右，右乘生成元的幂
                poly[k] = A + B; poly[k + h] = A - B;
                m *= o;
            }
        }
    }
    if (ifft)
        poly *= inv(T(CL));
}

template <typename T, int N>
void fht(TTensor<T, N>& s, bool ifht = false){
    static_assert(N == (N & -N), "N must be power of 2");
    fft_rev(s);
    for (int i = 0; i < N; i += 2){
        T t = s[i];
        s[i] = t + s[i+1];
        s[i+1] = t - s[i+1];
    }
    for (int s1 = 4; s1 <= N; s1 <<= 1){
        const int s2 = s1 >> 1;
        const int s4 = s1 >> 2;
        for (int k = 0; k < N; k += s1){
            int i0 = k;
            int i2 = i0 + s4;
            int i1 = i2 + s4;
            int i3 = i1 + s4;
            T t0 = s[i0] + s[i1];
            T t1 = s[i0] - s[i1];
            T t2 = s[i2] + s[i3];
            T t3 = s[i2] - s[i3];
            s[i0] = t0;
            s[i1] = t1;
            s[i2] = t2;
            s[i3] = t3;
        }
        const float angle = 2.0f * (float)PI / (float)s1;
        const float sv = sin(ifht ? -angle : angle), cv = cos(angle);
        // v = gen(T(), s1);
        float sk = sv, ck = cv;
        for (int j = 1; j < s4; j++){
            for (int k = 0; k < N; k += s1){
                int i0 = k + j;
                int i2 = k + s2 - j;
                int i1 = k + s2 + j;
                int i3 = k + s1 - j;
                T m0 = ck * s[i1] + sk * s[i3];
                T m1 = ck * s[i3] - sk * s[i1];
                T t0 = s[i0] + m0;
                T t1 = s[i0] - m0;
                T t2 = s[i2] - m1;
                T t3 = s[i2] + m1;
                s[i0] = t0;
                s[i1] = t1;
                s[i2] = t2;
                s[i3] = t3;
            }
            float t0 = ck * cv - sk * sv;
            float t1 = ck * sv + sk * cv;
            ck = t0; sk = t1;
        }
    }
    if (ifht)
        s *= inv(T(N));
}

#endif /* FFT_H */
