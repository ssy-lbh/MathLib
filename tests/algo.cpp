#include <cstdlib>

#include "math/math_base.h"
#include "math/matrix.h"
#include "math/number_theory.h"
#include "math/fft.h"
#include "math/fraction.h"
#include "math/complex.h"
#include "math/polynomial.h"
#include "math/crt.h"

void algo_fftmul(){
    TTensor<NMod<998244353>, 8> a = {5, 6, 1, 9, 0, 0, 0, 0};
    TTensor<NMod<998244353>, 8> b = {7, 6, 3, 0, 0, 0, 0, 0};
    TTensor<NMod<998244353>, 8> c = {5, 5, 5, 3, 6, 3, 3, 0};
    fft(a); fft(b);
    for (int i = 0; i < 8; i++)
        a[i] *= b[i];
    fft(a, true);
    for (int i = 0; i < 7; i++){
        a[i+1].n += a[i].n / 10;
        a[i].n %= 10;
        if (a[i] != c[i]){
            print(a);
            assert(false);
        }
    }
}

void algo_matacc(){
    // 0 => 1
    // 1 => 1
    // 2 => 2
    // ...
    auto calc_fab = [](int n){
        NMod<998244353> a = 1, b = 1;
        for (int i = 2; i <= n; i++){
            NMod<998244353> c = a + b;
            a = b;
            b = c;
        }
        return b;
    };
    TMatrix<NMod<998244353>, 2, 2> A = {{1, 1}, {1, 0}};
    for (int i = 0; i < 20; i++){
        int n = randmod(1000u);
        auto B = pow(A, n);
        auto a = B[0][0]; // * [1(0), 0(-1)]
        auto b = calc_fab(n);
        if (a != b){
            printf("%u %u %u\n", n, a.n, b.n);
            print(B);
            assert(false);
        }
    }
}

void algo_eigen(){
    const default_type eps = 1e-6;

    Matrix<2, 2> A = {{3, 1}, {1, 2}};
    Tensor<2> E;
    Matrix<2, 2> V;

    jacobi_eigen(A, E, V);
    
    assert(norm(A * V[0] - E[0] * V[0]) < eps);
    assert(norm(A * V[1] - E[1] * V[1]) < eps);
    assert(norm(A - inv(V) * diag(E) * V) < eps);
    assert(norm(trace(A) - sum(E)) < eps);
    assert(norm(det(A) - prod(E)) < eps);
}

void algo_frac_inv(){
    TMatrix<Fraction, 2, 2> A = {{8, 5}, {1, 4}};
    assert(det(A) == 27_f);
    assert(inv(A) * A == ident(A));

    TMatrix<Fraction, 2, 2> L, U;
    lu_decom(A, L, U);
    assert(L * U == A);
}

void algo_eigen_poly(){
    constexpr int N = 8;
    constexpr double eps = 1e-6;
    TMatrix<TComplex<double>, N, N> A;
    for (int t = 0; t < 10; t++){
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++){
                A[i][j][0] = ((int)randmod(64u) - 32) / 32.0;
                A[i][j][1] = ((int)randmod(64u) - 32) / 32.0;
            }
        assert(norm2(subs<0>(eigen_poly(A), A)) < eps);
    }
}

void algo_crt(){
    constexpr uint32_t N = 8;
    uint64_t mods[N] = {2, 3, 5, 7, 11, 13, 17, 19};
    for (int t = 0; t < 10; t++){
        uint64_t rems[N];
        for (int i = 0; i < N; i++)
            rems[i] = randmod(mods[i]);
        uint64_t ans = crt_calc(rems, mods, N);
        for (int i = 0; i < N; i++)
            assert(ans % mods[i] == rems[i]);
    }
}

void algo_solve3x3(){
    TMatrix<NMod<4>, 3, 3> A = {{1, 1, 0}, {1, 1, 1}, {0, 1, 1}};
    TMatrix<NMod<4>, 3, 3> IA = gauss_commutative_ring_inv(A);
    TMatrix<NMod<4>, 9, 9> B = nest(IA);
    TTensor<NMod<4>, 9> b = {0, 1, 3, 0, 2, 3, 1, 0, 2};
    print(B * b);
}

void algo_solve4x4(){
    TMatrix<NMod<4>, 4, 4> A = {
        {1, 1, 0, 0},
        {1, 1, 1, 0},
        {0, 1, 1, 1},
        {0, 0, 1, 1}
    };
    TMatrix<NMod<4>, 4, 4> IA = gauss_commutative_ring_inv(A);
    TMatrix<NMod<4>, 16, 16> B = nest(IA);
    TTensor<NMod<4>, 16> b = {1, 3, 2, 3, 1, 2, 1, 0, 3, 2, 2, 2, 0, 2, 1, 3};
    print(B * b);
}

const int N = 4 + 5 + 6 + 7 + 6 + 5 + 4;

TMatrix<int, N, N> gen_hexagon_matrix(){
    const int row[7] = {4, 5, 6, 7, 6, 5, 4};
    const int start[7] = {0, 4, 9, 15, 22, 28, 33};
    TMatrix<int, N, N> A = zero(TMatrix<int, N, N>());
    for (int i = 0, idx = 0; i < 7; i++){
        for (int j = 0; j < row[i]; j++, idx++){
            // 影响此行前后两个点，以及上下两行的两个点
            A[idx][idx] = 1;
            if (j > 0)
                A[idx][idx - 1] = 1;
            if (j + 1 < row[i])
                A[idx][idx + 1] = 1;
            if (i == 0){
                // 下行j,j+1
                A[idx][start[1] + j] = 1;
                A[idx][start[1] + j + 1] = 1;
            } else if (i < 3){
                // 上行j-1,j 下行j,j+1
                if (j > 0) A[idx][start[i - 1] + j - 1] = 1;
                if (j < row[i - 1]) A[idx][start[i - 1] + j] = 1;
                A[idx][start[i + 1] + j] = 1;
                A[idx][start[i + 1] + j + 1] = 1;
            } else if (i == 3){
                // 上行j-1,j 下行j-1,j
                if (j > 0) A[idx][start[i - 1] + j - 1] = 1;
                if (j < row[i - 1]) A[idx][start[i - 1] + j] = 1;
                if (j > 0) A[idx][start[i + 1] + j - 1] = 1;
                if (j < row[i + 1]) A[idx][start[i + 1] + j] = 1;
            } else if (i < 6){
                // 上行j,j+1 下行j-1,j
                A[idx][start[i - 1] + j] = 1;
                A[idx][start[i - 1] + j + 1] = 1;
                if (j > 0) A[idx][start[i + 1] + j - 1] = 1;
                if (j < row[i + 1]) A[idx][start[i + 1] + j] = 1;
            } else {
                // 上行j,j+1
                A[idx][start[5] + j] = 1;
                A[idx][start[5] + j + 1] = 1;
            }
        }
    }
    return A;
}

void algo_solve_hexagon(){
    // 整个图形是一个六边形，边长为4
    // 一个点的操作会影响自身和相邻的6个点
    // 点坐标编码: 共计7行，每行点数4,5,6,7,6,5,4，点坐标(0,0)在左上角
    const int N = 4 + 5 + 6 + 7 + 6 + 5 + 4;
    TMatrix<NMod<2>, N, N> A;
    convert(A, gen_hexagon_matrix());
    print(A);
}

int main(){
    // algo_fftmul();
    // algo_matacc();
    // algo_eigen();
    // algo_frac_inv();
    // algo_eigen_poly();
    // algo_crt();

    //algo_solve3x3();
    //algo_solve4x4();
    algo_solve_hexagon();

    return 0;
}
