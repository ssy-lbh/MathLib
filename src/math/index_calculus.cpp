#include "number_theory.h"

#include "sieves.h"

#include <algorithm>

/* 
* Index Calculus
* 求解离散对数问题, 即 a^x === b (mod p)
* 1. 预处理
*    1.1. 筛出一定数量的素数, 确保 b 能被这些素数完全分解, 素数的上界一般取 e^{\frac{\sqrt(\log mod * \log \log mod)}{2} + 1}
*    1.2. 这些素数依次编号为 p_i, 后文中令 e_i 为 p_i 的指数, a^{x_i} === p_i (mod p)
* 2. 列举方程
*    2.1. 令 t = randmod(0, p-1), 计算 a^t (mod p)
*    2.2. 以所有的 p 为基底, 对 a^t (mod p) 执行质因数分解, 得到 a^t (mod p) = \prod p_i^{e_i}, 如果无法分解为这些素数的乘积则舍去
*    2.3. 所得方程 t = \sum e_i * x_i
* 3. 解方程组
*    3.1. 采用模意义下的高斯消元法, 解出对于每个 p_i 的离散对数 x_i
* 4. 计算答案
*    4.1. 直接将 b 分解, 离散对数 x = \sum e_i * x_i
*/

struct IndexCalculusContext {
    uint64_t** a;
    uint64_t n;
    uint64_t m;

    uint64_t* prime;
    uint64_t* inv_prime;
    uint64_t* ind;
};

static uint32_t get_opti_value(uint64_t mod){
    double lgmod = log(mod);
    return (uint32_t)exp(sqrt(lgmod * log(lgmod)) * 0.5 + 1.0);
}

uint32_t index_calculus_init1(IndexCalculusContext& ctx, uint64_t mod){
    uint64_t limit = get_opti_value(mod);

    uint64_t* tag = new uint64_t[limit];
    uint64_t* prime = new uint64_t[pi_limit(limit)];
    memset(tag, 0, sizeof(uint64_t) * limit);
    uint64_t cnt = (uint64_t)euler_sieve(limit, tag, prime);
    delete[] tag;

    ctx.prime = new uint64_t[cnt];
    ctx.inv_prime = new uint64_t[cnt];
    ctx.ind = new uint64_t[cnt];
    memcpy(ctx.prime, prime, sizeof(uint64_t) * cnt);
    delete[] prime;

    ctx.n = cnt;
    ctx.m = cnt * 3;
    ctx.a = new uint64_t*[ctx.m];
    for (uint64_t i = 0; i < ctx.m; i++)
        ctx.a[i] = new uint64_t[cnt + 1];
    
    for (uint64_t i = 0; i < cnt; i++)
        ctx.inv_prime[i] = inv(ctx.prime[i], mod);

    return cnt;
}

bool gauss(IndexCalculusContext& ctx, uint64_t mod){
    for (uint32_t i = 0; i < ctx.n; i++){
        if (gcd(ctx.a[i][i], mod) != 1){
            uint32_t t = i;
            for (uint32_t j = i + 1; j < ctx.m; j++){
                if (gcd(ctx.a[j][i], mod) == 1){
                    t = j;
                    break;
                }
            }
            if (i == t)
                return false;
            std::swap(ctx.a[i], ctx.a[t]);
        }
        uint64_t cur_inv = inv(ctx.a[i][i], mod);
        for (uint32_t j = i + 1; j < ctx.m; j++){
            if (ctx.a[j][i]){
                uint64_t t = mul(cur_inv, ctx.a[j][i], mod);
                for (uint32_t k = i; k <= ctx.n; k++){
                    ctx.a[j][k] = (ctx.a[j][k] - mul(t, ctx.a[i][k], mod) + mod) % mod;
                }
            }
        }
    }
    for (int32_t i = ctx.n - 1; i >= 0; i--){
        for (uint32_t j = i + 1; j < ctx.n; j++){
            ctx.a[i][ctx.n] = (ctx.a[i][ctx.n] - mul(ctx.a[i][j], ctx.a[j][ctx.n], mod) + mod) % mod;
        }
        ctx.a[i][ctx.n] = mul(ctx.a[i][ctx.n], inv(ctx.a[i][i], mod), mod);
    }
    return true;
}

void index_calculus_init2(IndexCalculusContext& ctx, uint64_t g, uint64_t p){
    uint64_t phi = p - 1;
    do {
        uint32_t i = 0;
        for (uint64_t j = randmod(phi), k = pow(g, j, p); i < ctx.m; j = randmod(phi), k = pow(g, j, p)){
            for (uint64_t l = k, x = 1; x < phi; l = mul(l, k, p), x++){
                uint32_t t1 = ctzl(l);
                uint64_t t2 = l >> t1, t3;
                ctx.a[i][0] = t1;
                uint32_t y;
                for (y = 1; y < ctx.n && t2 > 1; y++){
                    if ((y == 9 && t2 > 1e15) || (y == 29 && t2 > 1e12))
                        break;
                    ctx.a[i][y] = 0;
                    while ((t3 = mul(t2, ctx.inv_prime[y], p)) < t2){
                        t2 = t3;
                        ctx.a[i][y]++;
                    }
                }
                memset(ctx.a[i] + y, 0, sizeof(uint64_t) * (ctx.n - y));
                if (t2 == 1){
                    ctx.a[i][ctx.n] = mul(j, x, phi);
                    i++;
                    break;
                }
            }
        }
    } while (!gauss(ctx, phi));
    for (uint32_t i = 0; i < ctx.n; i++)
        ctx.ind[i] = ctx.a[i][ctx.n];
}

inline uint64_t index_calculus(IndexCalculusContext& ctx, uint64_t a, uint64_t g, uint64_t p){
    if (a <= 1)
        return a - 1;
    uint64_t phi = p - 1;
    for (uint64_t i = randmod(phi), j = pow(g, i, p); ; i = randmod(phi), j = pow(g, i, p)){
        for (uint64_t k = j, l = 1; l < phi; k = mul(k, j, p), l++){
            uint32_t t1;
            uint64_t t2 = mul(a, k, p), t3;
            t1 = ctzl(t2);
            t2 >>= t1;
            uint64_t ans = (mul((uint64_t)t1, ctx.ind[0], phi) - mul(i, l, phi) + phi) % phi;
            for (uint32_t x = 1; x < ctx.n && t2 > 1; x++){
                if ((x == 9 && t2 > 1e15) || (x == 29 && t2 > 1e12))
                    break;
                while ((t3 = mul(t2, ctx.inv_prime[x], p)) < t2){
                    t2 = t3;
                    ans = (ans + ctx.ind[x]) % phi;
                }
            }
            if (t2 == 1){
                return ans;
            }
        }
    }
}

void index_calculus_free(IndexCalculusContext& ctx){
    for (uint64_t i = 0; i < ctx.m; i++)
        delete[] ctx.a[i];
    delete[] ctx.a;
    delete[] ctx.prime;
    delete[] ctx.inv_prime;
    delete[] ctx.ind;
}

// a^x === b (mod p)
uint64_t index_calculus_log(uint64_t a, uint64_t b, uint64_t g, uint64_t p){
    IndexCalculusContext ctx;
    index_calculus_init1(ctx, p);
    index_calculus_init2(ctx, g, p);
    if (b == 1)
        return 0;
    uint64_t phi = p - 1;
    uint64_t x = index_calculus(ctx, a, g, p);
    uint64_t d = gcd(x, phi);
    uint64_t y = index_calculus(ctx, b, g, p);
    index_calculus_free(ctx);
    if (y % d != 0)
        return -1;
    phi /= d;
    return mul(y / d, inv(x / d, phi), phi);
}

// g^x === a (mod p)
uint64_t index_calculus_log(uint64_t a, uint64_t g, uint64_t p){
    IndexCalculusContext ctx;
    index_calculus_init1(ctx, p);
    index_calculus_init2(ctx, g, p);
    uint64_t res = index_calculus(ctx, a, g, p);
    index_calculus_free(ctx);
    return res;
}
