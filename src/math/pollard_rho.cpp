#include "number_theory.h"

#include <queue>

uint64_t pollard_rho(uint64_t n){
    if (n == 4)
        return 2;
    if (miller_rabin(n))
        return n;
    while (true){
        uint64_t c = rand(1, n - 1); // 生成随机的c
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
        uint32_t c = rand(1, n - 1); // 生成随机的c
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

uint32_t pollard_rho(uint64_t n, uint64_t prime[], uint32_t exp[], uint32_t len){
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