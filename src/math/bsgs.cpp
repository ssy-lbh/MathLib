#include "number_theory.h"
#include "big_num.h"

#include <unordered_map>

// a^x = b (mod p)
uint64_t bsgs(uint64_t a, uint64_t b, uint64_t p){
    uint64_t m = ceil(sqrt(p));
    uint64_t am = pow(a, m, p);
    uint64_t inv = pow(am, p - 2, p);
    uint64_t cur = b;
    std::unordered_map<uint64_t, uint64_t> mp;
    for (uint64_t i = 0; i < m; i++){
        mp[cur] = i;
        cur = mul(cur, inv, p);
    }
    cur = 1;
    for (uint64_t i = 0; i < m; i++){
        if (mp.count(cur)){
            return i + mp[cur] * m;
        }
        cur = mul(cur, a, p);
    }
    return -1;
}

struct BigIntHash {
    size_t operator()(const BigInt& p) const {
        return std::hash<std::string>()(to_string(p, 16));
    }
};

// a^x = b (mod p)
BigInt bsgs(const BigInt& a, const BigInt& b, const BigInt& p){
    BigInt m = sqrt(p) + 1;
    BigInt am = pow(a, m, p);
    BigInt inv = pow(am, p - 2, p);
    BigInt cur = b;
    std::unordered_map<BigInt, BigInt, BigIntHash> mp;
    for (BigInt i = 0; i < m; i++){
        mp[cur] = i;
        cur = mul(cur, inv, p);
    }
    cur = 1;
    for (BigInt i = 0; i < m; i++){
        if (mp.count(cur)){
            return i + mp[cur] * m;
        }
        cur = mul(cur, a, p);
    }
    return -1;
}