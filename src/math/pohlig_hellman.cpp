#include "number_theory.h"

#include "crt.h"

constexpr uint32_t N = 16;

// g^x === b (mod p)
uint64_t pohlig_hellman(uint64_t g, uint64_t b, uint64_t p, uint64_t pm1_prime[], uint32_t exp[], uint32_t cnt){
    uint32_t numc = 0;
    uint64_t crta[N];
    uint64_t crtb[N];
    uint64_t pow_prime[64];
    uint64_t factor[64];

    uint64_t phi = p - 1;

    for (uint32_t i = 0; i < cnt; i++) { // 枚举每个质因子
        memset(factor, 0, sizeof(factor));
        pow_prime[0] = 1;
        for (uint32_t j = 0; j <= exp[i]; j++){
            pow_prime[j+1] = pow_prime[j] * pm1_prime[i]; // 预处理
        }
        uint64_t sum = 1;
        uint64_t rem = 0, mod = pow_prime[exp[i]];
        for (uint32_t j = 1; j <= exp[i]; j++){ //求出每个系数
            uint64_t a0 = pow(g, phi / pm1_prime[i], p);
            uint64_t b0 = pow(b, phi / pow_prime[j], p);
            uint64_t ap = 1;
            for (uint64_t x = 0; x <= pm1_prime[i] - 1; x++){ // 遍历找出系数x
                if (mul(sum, ap, p) == b0){
                    factor[j] = x;
                    uint64_t xs = 0; // 已求出的系数按权求和
                    for (uint64_t k = 0; k < j; k++)
                        xs = (xs + (pow_prime[k] * factor[k + 1] % phi)) % phi;
                    sum = pow(g, mul(xs, phi / pow_prime[j + 1], phi), p);
                    break;
                }
                ap = mul(ap, a0, p);
            }
            rem += pow_prime[j - 1] * factor[j];
        }
        if (exp[i]){
            crta[numc] = rem;
            crtb[numc] = mod;
            numc++;
        }
    }
    return crt_calc(crta, crtb, numc, phi);
}

// g^x === b (mod p)
uint64_t pohlig_hellman(uint64_t g, uint64_t b, uint64_t p){
    uint64_t phi = p - 1;

    uint64_t pm1_prime[N];
    uint32_t exp[N];

    uint32_t cnt = pollard_rho(phi, pm1_prime, exp, N);

    return pohlig_hellman(g, b, p, pm1_prime, exp, cnt);
}

// a^x === b (mod p)
uint64_t pohlig_hellman_log(uint64_t a, uint64_t b, uint64_t p){
    uint64_t phi = p - 1;

    uint64_t pm1_prime[N];
    uint32_t exp[N];

    uint32_t cnt = pollard_rho(phi, pm1_prime, exp, N);

    if (check_root(a, p, pm1_prime, cnt))
        return pohlig_hellman(a, b, p, pm1_prime, exp, cnt);

    if (!has_root(pm1_prime, exp, cnt))
        return -1;
    uint64_t g = find_root(p, pm1_prime, cnt);

    uint64_t lga = pohlig_hellman(g, a, p, pm1_prime, exp, cnt);
    uint64_t lgb = pohlig_hellman(g, a, p, pm1_prime, exp, cnt);

    int64_t x, y;
    uint64_t fac = exgcd(lga, phi, x, y);

    if (lgb % fac != 0)
        return -1;
    
    uint64_t px = x < 0 ? x + phi : x;
    return mul(px, lgb / fac, phi / fac);
}
