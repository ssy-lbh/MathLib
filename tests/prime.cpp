#include "math/number_theory.h"

#include "math/sieves.h"
#include "math/prime_count.h"

#include <algorithm>

constexpr uint32_t N = 20;

void test_pollard_rho(){
    uint64_t n = 673984; // 10531 * 2^6

    uint64_t prime[N];
    uint32_t exp[N];
    uint32_t cnt = pollard_rho(n, prime, exp, N);

    assert(cnt == 2);
    if (prime[0] > prime[1]){
        std::swap(prime[0], prime[1]);
        std::swap(exp[0], exp[1]);
    }
    assert(prime[0] == 2);
    assert(prime[1] == 10531);
    assert(exp[0] == 6);
    assert(exp[1] == 1);
}

void test_pohlig_hellman_log(){
    uint64_t g = 2;
    uint64_t b = 3;
    uint64_t p = 5;

    uint64_t prime[N];
    uint32_t exp[N];
    uint32_t cnt = pollard_rho(p - 1, prime, exp, N);

    assert(check_root(g, p, prime, cnt));
    uint64_t x = pohlig_hellman_log(g, b, p);
    assert(x == 3);
}

void test_index_calculus_log(){
    uint64_t g = 3;
    uint64_t b = 524234;
    uint64_t p = 998244353;

    uint64_t prime[N];
    uint32_t exp[N];
    uint32_t cnt = pollard_rho(p - 1, prime, exp, N);

    assert(check_root(g, p, prime, cnt));
    uint64_t x = index_calculus_log(b, g, p);
    printf("%llu\n", x);
}

void test_dujiao_sieve(){
    int64_t n = 4729389;
    // mobius function
    int64_t res = dujiao_sieve(n, [](int64_t x) -> int64_t { return 1; });
    assert(res == 61);
    // unit function
    res = dujiao_sieve(n, [](int64_t x) -> int64_t { return x; });
    assert(res == 1);
    // euler function
    res = dujiao_sieve(n, [](int64_t x) -> int64_t { return x * (x + 1) / 2; });
    assert(res == 6798790501158ll);
}

int main(){
    test_pollard_rho();
    test_pohlig_hellman_log();
    //test_index_calculus_log();
    test_dujiao_sieve();
    return 0;
}