#include <cstdlib>

#include "math/math_base.h"
#include "math/matrix.h"
#include "math/number_theory.h"
#include "math/fft.h"
#include "math/fraction.h"
#include "math/complex.h"

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
        int n = rand(1000u);
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

int main(){
    algo_fftmul();
    algo_matacc();
    algo_eigen();
    algo_frac_inv();
    return 0;
}
