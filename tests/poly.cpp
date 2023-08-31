#include "math/number_theory.h"
#include "math/polynomial.h"
#include "math/fft.h"

void test_poly_mul(){
    const uint32_t MOD = 998244353;
    const uint32_t G = 3;

    const uint32_t N = 16;
    uint32_t poly1[N << 1]{};
    uint32_t poly2[N << 1]{};
    uint32_t target[N << 1]{};

    for (uint32_t i = 0; i < 20; i++){
        for (uint32_t j = 0; j < N; j++){
            poly1[j] = rand() % MOD;
            poly2[j] = rand() % MOD;
        }
        memset(poly1 + N, 0, sizeof(uint32_t) * N);
        memset(poly2 + N, 0, sizeof(uint32_t) * N);
        memset(target, 0, sizeof(uint32_t) * (N << 1));
        for (uint32_t j = 0; j < N; j++){
            for (uint32_t k = 0; k < N; k++){
                target[j + k] = (target[j + k] + (uint64_t)poly1[j] * poly2[k]) % MOD;
            }
        }
        poly_mul(poly1, poly2, N, G, MOD);
        for (uint32_t j = 0; j < (N << 1); j++){
            assert(poly1[j] == target[j]);
        }
    }
}

void test_poly_inv(){
    const uint32_t MOD = 998244353;
    const uint32_t G = 3;

    const uint32_t N = 16;
    uint32_t poly1[N << 1]{};
    uint32_t poly2[N << 1]{};
    uint32_t tmp[N << 1]{};

    for (uint32_t i = 0; i < 20; i++){
        for (uint32_t j = 0; j < N; j++){
            poly1[j] = poly2[j] = rand() % MOD;
        }
        memset(poly1 + N, 0, sizeof(uint32_t) * N);
        memset(poly2 + N, 0, sizeof(uint32_t) * N);
        poly1[0] = poly2[0] = 1;
        poly_inv(poly1, tmp, N, G, MOD);
        poly_mul(poly1, poly2, N, G, MOD);
        for (uint32_t j = 0; j < N; j++){
            //printf("%u ", poly1[j]);
            assert(poly1[j] == (j == 0));
        }
    }
}

void test_high_order_fft(){
    using T = NMod<998244353>;
    constexpr uint32_t N = 4;
    TTensor<T, N, N, N> A{};
    TTensor<T, N, N, N> B{};
    A[1][0][2] = 1;
    B[0][3][1] = 1;
    fft(A); fft(B);
    TTensor<T, N, N, N> C = hadamard(A, B);
    fft(C, true);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                assert(C[i][j][k] == (i == 1 && j == 3 && k == 3));
}

int main(){
    test_poly_mul();
    test_poly_inv();
    test_high_order_fft();
    return 0;
}