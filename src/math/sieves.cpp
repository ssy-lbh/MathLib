#include "number_theory.h"

#include <unordered_map>

uint32_t euler_sieve(uint32_t n, uint32_t fac[], uint32_t prime[]){
    uint32_t cnt = 0;
    for (uint32_t i = 2; i < n; i++){
        if (!fac[i]){
            fac[i] = i;
            prime[cnt++] = i;
        }
        for (uint32_t j = 0; j < cnt; j++){
            if ((uint64_t)i * prime[j] >= n || prime[j] > fac[i]) break;
            fac[i * prime[j]] = prime[j];
        }
    }
    return cnt;
}

uint64_t euler_sieve(uint64_t n, uint64_t fac[], uint64_t prime[]){
    uint64_t cnt = 0;
    for (uint64_t i = 2; i < n; i++){
        if (!fac[i]){
            fac[i] = i;
            prime[cnt++] = i;
        }
        for (uint64_t j = 0; j < cnt; j++){
            if (i * prime[j] >= n || prime[j] > fac[i]) break;
            fac[i * prime[j]] = prime[j];
        }
    }
    return cnt;
}

uint32_t egypt_sieve(uint32_t n, bool tag[], uint32_t prime[]){
    uint32_t cnt = 0;
    for (uint32_t i = 2; i < n; i++){
        if (!tag[i]){
            prime[cnt++] = i;
            for (uint64_t j = i * i; j < n; j += i)
                tag[j] = true;
        }
    }
    return cnt;
}

uint64_t egypt_sieve(uint64_t n, bool tag[], uint64_t prime[]){
    uint64_t cnt = 0;
    for (uint64_t i = 2; i < n; i++){
        if (!tag[i]){
            prime[cnt++] = i;
            if (i >= (n + i - 1) / i)
                continue;
            for (uint64_t j = i * i; j < n; j += i)
                tag[j] = true;
        }
    }
    return cnt;
}

// min25筛 O(n^{2/3})

// 杜教筛
// 求解 \sum_{i=1}^n f(i)
// 输入 h = \sum{i=1}^n f * g
// 默认 g = 1

// \phi * 1 = id
// \mu * 1 = \epsilon
// \phi * \mu = id

int64_t dujiao_sum(int64_t x, std::unordered_map<int64_t, int64_t>& cache, int64_t(*h)(int64_t)){
    if (cache.count(x)) return cache[x];
    int64_t res = h(x);
    for (int64_t l = 2, r; l <= x; l = r + 1){
        int64_t y = x / l;
        // r = x / (x / l), [x / n] == y 求出的n区间
        r = x / y;
        res -= (r - l + 1) * dujiao_sum(y, cache, h);
    }
    return cache[x] = res;
}

int64_t dujiao_sieve(int64_t n, int64_t(*h)(int64_t)){
    std::unordered_map<int64_t, int64_t> cache;
    return dujiao_sum(n, cache, h);
}
