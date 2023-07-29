#include <cstdlib>

#include "math/math_base.h"
#include "math/big_num.h"
#include "math/matrix.h"
#include "math/number_theory.h"
#include "math/fraction.h"
#include "math/elliptic_curve.h"

void bignum_pow(){
    BigFloat a = 23;
    BigFloat b = 22;
    assert(pow(a, b) == "907846434775996175406740561329");
}

void bignum_arith(){
    BigInt A = "8866128975287528";
    BigInt B = "-8778405442862239";
    BigInt C = "-2736111468807040";
    assert(A * A * A + B * B * B + C * C * C == 33);
}

void bignum_eigen(){
    const double eps = 1e-20;
    
    TMatrix<BigFloat, 2, 2> A = {{1, 1}, {1, 0}};
    TTensor<BigFloat, 2> E;
    TMatrix<BigFloat, 2, 2> V;

    jacobi_eigen(A, E, V, eps, 1000);

    assert(norm(A * V[0] - E[0] * V[0]) < eps);
    assert(norm(A * V[1] - E[1] * V[1]) < eps);
    assert(norm(A - inv(V) * diag(E) * V) < eps);
    assert(norm(trace(A) - sum(E)) < eps);
    assert(norm(det(A) - prod(E)) < eps);
}

/* assert(riemann_hypothesis()); */
void bignum_prime(){
    BigInt A = "300000000000000000053";
    assert(miller_prime_proof(A));
    A = "900900900900900900900900900900990990990990990990990990990991";
    assert(miller_prime_proof(A));
    A = "1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111";
    assert(!miller_prime_proof(A));
    A = "153211620887015423991278431667808361439217294295901387715486473457925534859044796980526236853";
    assert(miller_prime_proof(A));
}

void bignum_test_ec(){
    TEllipticCurve<BigFrac> ec(1, 6);
    Point2<BigFrac> p = {3, 6}, q = {2, 4};
    assert(ec.is_on_curve(p));
    assert(ec.is_on_curve(q));
    assert(ec.is_on_curve(ec.add(p, q))); // [-1 2]
    assert(ec.is_on_curve(ec.add(q, q))); // [-87/64 747/512]
    assert(ec.is_on_curve(ec.add(p, p))); // [-5/9 62/27]

    // select a random point G on the curve
    // 1. G is a generator of the group
    // 2. G is a point on the curve
    // 3. G has a large prime order
    // generate a random number r
    // compute Q = r * G
    // Q is a public key
    // r is a private key
    // to encrypt a message M, compute k * G = (x1, y1)
    // compute M + k * Q = (x2, y2)
    // to decrypt, compute (x1, y1) - r * (x2, y2)

    TEllipticCurve<NMod<257>> ec2(3, 19);
    for (int i = 0; i < 10; i++)
        assert(ec2.is_on_curve(ec2.random()));
}

void bignum_pell_equation(){
    BigInt t = 39352;
    BigInt m = sqrt(t);
    std::vector<BigInt> f;
    f.push_back(m);
    BigInt a = m, b = 1;
    while (true){
        b = (t - a * a) / b;
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
    assert(q * q - t * p * p == 1);
    assert(q == "2662019309411216232806345449321879270495478346383");
    assert(p == "13419236180444562206941613278603693768762410762");
    // all solutions [[q, t * p], [p, q]]^n * [1 0]
}

int main(){
    bignum_pow();
    bignum_arith();
    bignum_eigen();
    bignum_prime();
    bignum_test_ec();
    bignum_pell_equation();
    return 0;
}
