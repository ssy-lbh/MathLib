#include "elliptic_curve.h"
#include "number_theory.h"
#include "big_num.h"
#include "sieves.h"

#include <cstdlib>

static TTensor<uint64_t, 2> compute_bounds(BigInt n){
    uint64_t log_n = sizeinbase(n, 10);
    if (log_n <= 30)
        return {2000, 147396};
    else if (log_n <= 40)
        return {11000, 1873422};
    else if (log_n <= 50)
        return {50000, 12746592};
    else if (log_n <= 60)
        return {250000, 128992510};
    else if (log_n <= 70)
        return {1000000, 1045563762};
    else if (log_n <= 80)
        return {3000000, 5706890290};
    return {430000000, 20000000000};
}
 
static Point2<BigInt> point_add(Point2<BigInt> p, Point2<BigInt> q, Point2<BigInt> r, BigInt n){
    BigInt u = (p[0] - p[1]) * (q[0] + q[1]);
    BigInt v = (p[0] + p[1]) * (q[0] - q[1]);
    BigInt upv = u + v;
    BigInt umv = u - v;
    BigInt x = mul(r[1], upv * upv, n);
    BigInt z = mul(r[0], umv * umv, n);
    return {x, z};
}

static Point2<BigInt> point_double(Point2<BigInt> p, BigInt n, BigInt a24){
    BigInt u = p[0] + p[1];
    BigInt v = p[0] - p[1];
    BigInt u2 = u * u;
    BigInt v2 = v * v;
    BigInt t = (u2 >= v2 ? u2 - v2 : u2 + n - v2);
    BigInt x = mul(u2, v2, n);
    BigInt z = t * (v2 + a24 * t) % n;
    return {x, z};
}

static Point2<BigInt> scalar_multiply(BigInt k, Point2<BigInt> p, BigInt n, BigInt a24){
    uint64_t lk = sizeinbase(k, 2);
    Point2<BigInt> q = p;
    Point2<BigInt> r = point_double(p, n, a24);
 
    for (uint64_t i = 1; i < lk; i++){
        if (k.testbit(lk - i - 1)){
            q = point_add(r, q, p, n);
            r = point_double(r, n, a24);
        } else {
            r = point_add(q, r, p, n);
            q = point_double(q, n, a24);
        }
    }

    return q;
}

static uint64_t pi_limit(uint64_t x){
    double log_x = log(x);
    double log_log_x = log(log_x);
    return (uint64_t)(x / log_x * (1 + 1.2762 / log_x + 1.2762 * log_log_x / log_x / log_x));
}

const uint64_t MAX_CURVES_ECM = 10000;
const uint64_t MAX_RND_ECM = 0x8000000000000000ULL;

BigInt ecm_factorize(BigInt n){
    if (n <= 1 || miller_prime_proof(n))
        return n;

    auto bound = compute_bounds(n);
    uint64_t D = sqrt_floor(bound[1]);

    BigInt* beta = new BigInt[D];
    Point2<BigInt>* S = new Point2<BigInt>[D];

    uint64_t curves = 0;

    uint64_t lim = pi_limit(bound[1]);
    uint64_t* primes = new uint64_t[lim];
    bool* tag = new bool[bound[1]];
    uint64_t num_primes = egypt_sieve(bound[1], tag, primes);
    delete[] tag;
    assert(num_primes <= lim);
    uint64_t idx_B0 = std::lower_bound(primes, primes + num_primes, bound[0]) - primes;

    BigInt k = 1;
    double log_B0 = log(bound[0]);
    for (uint64_t i = 0; i < idx_B0; i++){
        uint64_t p = primes[i];
        k *= pow(BigInt(p), (uint64_t)(log_B0 / log(p)));
    }
 
    BigInt g = 1;
    while ((g == 1 || g == n) && curves <= MAX_CURVES_ECM){
        curves++;
        BigInt sigma = rand(6, MAX_RND_ECM);

        // Generate a new random curve in Montgomery form with Suyama's parametrization
        BigInt u = ((sigma * sigma) - 5) % n;
        BigInt v = mul(4, sigma, n);
        BigInt vmu = v - u;
        BigInt A = ((vmu * vmu * vmu) * (3 * u + v) / (4 * u * u * u * v) - 2) % n;
        if (A < 0) A += n;
        BigInt a24 = (A + 2) >> 2;
 
        // ----- Stage 1 -----
        Point2<BigInt> p = {((u * u * u) / (v * v * v)) % n, 1};
        Point2<BigInt> q = scalar_multiply(k, p, n, a24);
        g = gcd(n, q[1]);
 
        // If stage 1 is successful, return a non-trivial factor else
        // move on to stage 2
        if (g != 1 && g != n){
            delete[] beta;
            delete[] S;
            delete[] primes;
            return g;
        }
 
        // ----- Stage 2 -----
        S[0] = point_double(q, n, a24);
        S[1] = point_double(S[0], n, a24);
        beta[0] = mul(S[0][0], S[0][1], n);
        beta[1] = mul(S[1][0], S[1][1], n);
        for (uint64_t d = 2; d < D; d++){
            S[d] = point_add(S[d - 1], S[0], S[d - 2], n);
            beta[d] = mul(S[d][0], S[d][1], n);
        }

        BigInt g = 1;
        uint64_t B = bound[0] - 1;
 
        auto r = scalar_multiply(B, q, n, a24);
        auto t = scalar_multiply(B - (D << 1), q, n, a24);
        uint64_t idx = idx_B0;
        uint64_t step = D << 1;
 
        for (uint64_t ridx = B; ridx < bound[1]; ridx += step){
            auto alpha = mul(r[0], r[1], n);
            uint64_t limit = ridx + step;
            while (idx < num_primes && primes[idx] <= limit){
                uint64_t d = ((primes[idx] - ridx) >> 1) - 1;
                auto f = (r[0] - S[d][0]) * (r[1] + S[d][1]) - alpha + beta[d];
                g = mul(g, f, n);
                idx++;
            }
            auto tr = r;
            r = point_add(r, S[D - 1], t, n);
            t = tr;
        }
        g = gcd(n, g);
    }
    delete[] beta;
    delete[] S;
    delete[] primes;
    if (curves > MAX_CURVES_ECM)
        return -1;
    else
        return g;
}