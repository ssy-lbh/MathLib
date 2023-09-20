#include <cstdlib>

#include "math/math_base.h"
#include "math/big_num.h"
#include "math/matrix.h"
#include "math/number_theory.h"
#include "math/fraction.h"
#include "math/elliptic_curve.h"
#include "math/montgomery_mul.h"

void bignum_pow(){
    BigFloat a = 23.;
    BigFloat b = 22.;
    assert(pow(a, b) == "907846434775996175406740561329");
}

void bignum_arith(){
    BigInt A = "8866128975287528";
    BigInt B = "-8778405442862239";
    BigInt C = "-2736111468807040";
    assert(A * A * A + B * B * B + C * C * C == 33);
}

void bignum_eigen(){
    const BigFloat eps = 1e-20;
    
    TMatrix<BigFloat, 2, 2> A = {{1., 1.}, {1., 0.}};
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

    TEllipticCurve<NMod<257>> ec2(3, 19);
    for (int i = 0; i < 10; i++)
        assert(ec2.is_on_curve(ec2.random()));
}

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
void bignum_test_ecc_enc(){
    const BigInt MOD = "6277101735386680763835789423207666416083908700390324961279";
    assert(miller_prime_proof(MOD));
    TEllipticCurve<TMod<BigInt>> EC({1, MOD}, {6, MOD});
    Point2<TMod<BigInt>> G = {{2, MOD}, {4, MOD}}; // generator
    assert(EC.is_on_curve(G));
    BigInt R = randmod(MOD); // private key
    Point2<TMod<BigInt>> Q = EC.mul(R, G); // public key
    Point2<TMod<BigInt>> M = EC.random(); // message
    BigInt K = randmod(MOD);
    Point2<TMod<BigInt>> KG = EC.mul(K, G);
    Point2<TMod<BigInt>> KQ = EC.mul(K, Q);
    Point2<TMod<BigInt>> C = EC.add(M, KQ); // encrypted message
    Point2<TMod<BigInt>> D = EC.add(C, EC.neg(EC.mul(R, KG))); // decrypted message
    // EC: y^2 = x^3 + x + 6 (mod 6277101735386680763835789423207666416083908700390324961279)
    // private key: 2257266665275028717368858634151775288157095944676107451977
    // generator: [2 4]
    // encrypted point pair: ([364562710969197180833286231649197863335006589095609232292 4064846720840770304970884966639104197186368667790994813115], [4345212677055395448941375689444562490478462710431726579977 154092500694242968000986641518415931559619348920495233257]) (kG, k(rG) + M)
    // message: [842479473113610285124554190380454783966523983876323970880 5424199499805154532289260646877346460171942119518172451807]
    assert(D == M);
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

void bignum_factorize(){
    BigInt n = "5079161913028937957193211717202415680711159149463252";
    BigInt prime[100];
    uint64_t exp[100];
    uint32_t len = factorize(n, prime, exp, 100);
    for (uint32_t i = 0; i < len; i++){
        print(prime[i]);
        print("^");
        print(exp[i]);
        print(" ");
    }
}

void bignum_montgomery_mul(){
    BigInt mod = 123456789;
    MontgomeryMul<BigInt> mm(mod);
    for (int i = 0; i < 100; i++){
        BigInt a = randmod(mod);
        BigInt b = randmod(mod);
        BigInt c = mm.modmul(a, b);
        assert(c == a * b % mod);
    }
}

void bignum_discrete_log(){
    BigInt p = "11269160459";
    // BigInt p = "3087913911303507899";

    BigInt prime[64];
    uint64_t exp[64];
    uint64_t cnt = factorize(p - 1, prime, exp, 64);
    BigInt g = find_root(p, prime, cnt);

    BigInt b = randmod(p - 1) + 1;

    BigInt x = index_calculus_log(g, b, p);
    print("ind_");
    print(g);
    print("[");
    print(b);
    print("] = ");
    print(x);
    print("\n");
    assert(pow(g, x, p) == b);
    // ind_2[2049958631902929026] = 3041110249776312847 (mod 3087913911303507899)
}

int main(){
    bignum_pow();
    bignum_arith();
    bignum_eigen();
    bignum_prime();
    bignum_test_ec();
    bignum_pell_equation();
    //bignum_factorize();
    bignum_test_ecc_enc();
    bignum_montgomery_mul();
    //bignum_discrete_log();
    return 0;
}
