#include <cstdlib>

#include "math/math_base.h"
#include "math/big_num.h"
#include "math/matrix.h"
#include "math/number_theory.h"
#include "math/fft.h"
#include "math/fraction.h"
#include "math/complex.h"

void bignum_pow(){
    BigFloat a = 223;
    BigFloat b = 423;
    println(pow(a, b));
}

void bignum_eigen(){
    TMatrix<BigFloat, 2, 2> A = {{1, 1}, {1, 0}};
    TTensor<BigFloat, 2> E;
    TMatrix<BigFloat, 2, 2> V;
    jacobi_eigen(A, E, V);
    println(E);
    println(V);
}

int main(){
    bignum_pow();
    bignum_eigen();
    return 0;
}
