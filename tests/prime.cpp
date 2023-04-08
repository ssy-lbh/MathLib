#include "math/number_theory.h"

constexpr uint32_t N = 16;

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

void test_discrete_log(){
    uint64_t g = 2;
    uint64_t b = 3;
    uint64_t p = 5;

    uint64_t prime[N];
    uint32_t exp[N];
    uint32_t cnt = pollard_rho(p - 1, prime, exp, N);

    assert(check_root(g, p, prime, cnt));
    uint64_t x = pohlig_hellman_log(g, b, p);
    printf("%llu\n", x);
    assert(x == 3);
}

int main(){
    test_pollard_rho();
    test_discrete_log();
    return 0;
}