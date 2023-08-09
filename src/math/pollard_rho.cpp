#include "number_theory.h"
#include "big_num.h"
#include "ecm_factorize.h"
#include "sieves.h"

#include <queue>

uint64_t pollard_rho(uint64_t n){
    if (n == 4)
        return 2;
    if (miller_rabin(n))
        return n;
    while (true){
        uint64_t c = rand(1, n); // 生成随机的c
        auto f = [=](uint64_t x) {
            uint64_t y = (uint64_t)mul(x, x, n) + c;
            return (y >= n ? y - n : y);
        }; // lll表示__int128，防溢出
        uint64_t t = 0, r = 0, p = 1, q;
        do{
            for (int i = 0; i < 128; i++){ // 令固定距离C=128
                t = f(t), r = f(f(r));
                if (t == r || (q = mul(p, abs((int64_t)t - (int64_t)r), n)) == 0) // 如果发现环，或者积即将为0，退出
                    break;
                p = q;
            }
            uint64_t d = gcd(p, n);
            if (d > 1)
                return d;
        } while (t != r);
    }
}

uint32_t pollard_rho(uint32_t n){
    if (n == 4)
        return 2;
    if (miller_rabin(n))
        return n;
    while (true){
        uint32_t c = rand(1, n); // 生成随机的c
        auto f = [=](uint32_t x) {
            uint32_t y = (uint64_t)x * x % n + c;
            return (y >= n ? y - n : y);
        };
        uint32_t t = 0, r = 0, p = 1, q;
        do{
            for (int i = 0; i < 16; i++){ // 令固定距离C=16
                t = f(t), r = f(f(r));
                if (t == r || (q = (uint64_t)p * (uint32_t)abs((int32_t)t - (int32_t)r) % n) == 0) // 如果发现环，或者积即将为0，退出
                    break;
                p = q;
            }
            uint32_t d = gcd(p, n);
            if (d > 1)
                return d;
        } while (t != r);
    }
}

uint32_t factorize(uint64_t n, uint64_t prime[], uint32_t exp[], uint32_t len){
    if (n <= 1)
        return 0;
    uint32_t cnt = 0;
    std::queue<uint64_t> q;
    q.push(n);
    while (!q.empty()){
        uint64_t x = q.front();
        q.pop();
        if (x <= 1)
            continue;
        uint64_t d;
        for (uint32_t i = 0; i < cnt; i++)
            while (x % prime[i] == 0){
                exp[i]++;
                x /= prime[i];
                if (x <= 1)
                    goto next;
            }
        d = (x <= 0xFFFFFFFFull ? pollard_rho((uint32_t)x) : pollard_rho(x));
        if (d == x){
            for (uint32_t i = 0; i < cnt; i++)
                if (prime[i] == x){
                    exp[i]++;
                    goto next;
                }
            if (cnt == len)
                return cnt;
            prime[cnt] = x;
            exp[cnt++] = 1;
            continue;
        }
        q.push(d);
        q.push(x / d);
        next:;
    }
    return cnt;
}

BigInt pollard_rho(BigInt n){
    if (n == 4)
        return 2;
    if (miller_prime_proof(n))
        return n;
    while (true){
        BigInt c = randmod(n) + 1; // 生成随机的c
        auto f = [=](BigInt x) {
            BigInt y = (BigInt)mul(x, x, n) + c;
            return (y >= n ? y - n : y);
        };
        BigInt t = 0, r = 0, p = 1, q;
        do {
            for (int i = 0; i < 128; i++){ // 令固定距离C=128
                t = f(t), r = f(f(r));
                if (t == r || (q = mul(p, abs(t - r), n)) == 0) // 如果发现环，或者积即将为0，退出
                    break;
                p = q;
            }
            BigInt d = gcd(p, n);
            if (d > 1)
                return d;
        } while (t != r);
    }
}

const uint64_t PRECOMPUTE_LIMIT = 25000ULL;
const BigInt POLLARD_RHO_MAX = "100000000000000000000";

uint64_t factorize_hard(BigInt n, BigInt prime[], uint64_t exp[], uint64_t len, uint64_t cnt){
    if (n <= 1)
        return 0;
    std::queue<BigInt> q;
    q.push(n);
    while (!q.empty()){
        BigInt x = q.front();
        q.pop();
        if (x <= 1)
            continue;
        BigInt d;
        if (x > POLLARD_RHO_MAX){
            d = ecm_factorize(x);
            if (d == -1){
                d = pollard_rho(x);
            }
        } else {
            d = pollard_rho(x);
        }
        if (d == x){
            for (uint64_t i = 0; i < cnt; i++)
                if (prime[i] == x){
                    exp[i]++;
                    goto next;
                }
            if (cnt == len)
                return cnt;
            prime[cnt] = x;
            exp[cnt++] = 1;
            continue;
        }
        q.push(d);
        q.push(x / d);
        next:;
    }
    return cnt;
}

uint64_t factorize(BigInt n, BigInt prime[], uint64_t exp[], uint64_t len, uint64_t filter){
    uint64_t cnt = 0;

    if (filter) {
        uint64_t lim = pi_limit(filter);
        uint64_t* primes = new uint64_t[lim];
        bool* tag = new bool[filter];
        memset(tag, false, sizeof(bool) * filter);
        uint64_t num_primes = egypt_sieve(filter, tag, primes);
        delete[] tag;
        assert(num_primes <= lim);

        for (uint64_t i = 0; i < num_primes; i++){
            if (n % primes[i] == 0){
                prime[cnt] = primes[i];
                exp[cnt] = 1;
                n /= primes[i];
                while (n % primes[i] == 0){
                    exp[cnt]++;
                    n /= primes[i];
                }
                cnt++;
                if (cnt == len)
                    return cnt;
            }
        }
        delete[] primes;
    }

    return factorize_hard(n, prime, exp, len, cnt);
}

uint64_t factorize(BigInt n, BigInt prime[], uint64_t exp[], uint64_t len){
    uint64_t bit = sizeinbase(n, 2) / 5;
    uint64_t filter = (bit >= 64 ? PRECOMPUTE_LIMIT : ((1ULL << bit) > PRECOMPUTE_LIMIT ? PRECOMPUTE_LIMIT : (1ULL << bit)));
    if (filter < 8)
        filter = 8;
    return factorize(n, prime, exp, len, filter);
}
