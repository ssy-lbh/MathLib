#include "prime_count.h"

#include "sieves.h"

constexpr uint32_t M = 2;
constexpr uint32_t PM = 2 * 3 * 5;

static uint32_t N;

static bool* np;
static uint32_t* prime;
static uint32_t* pi;
static uint32_t phi[PM + 1][M + 1];
static uint32_t sz[M + 1];

class PrimeCount {
    #define MAXN 100
    #define MAXM 10001
    #define MAXP 40000
    #define MAX 400000
    #define clr(ar) memset(ar, 0, sizeof(ar))
    #define dbg(x) cout << #x << " = " << x << endl
    #define chkbit(ar, i) (((ar[(i) >> 6]) & (1 << (((i) >> 1) & 31))))
    #define setbit(ar, i) (((ar[(i) >> 6]) |= (1 << (((i) >> 1) & 31))))
    #define isprime(x) (( (x) && ((x)&1) && (!chkbit(ar, (x)))) || ((x) == 2))

    uint64_t dp[MAXN][MAXM];
    uint32_t ar[(MAX >> 6) + 5] = {0};
    int len = 0, primes[MAXP], counter[MAX];

    void sieve(){
        setbit(ar, 0), setbit(ar, 1);
        for (int i = 3; (i * i) < MAX; i++, i++){
            if (!chkbit(ar, i)){
                int k = i << 1;
                for (int j = (i * i); j < MAX; j += k) setbit(ar, j);
            }
        }

        for (int i = 1; i < MAX; i++){
            counter[i] = counter[i - 1];
            if (isprime(i)) primes[len++] = i, counter[i]++;
        }
    }

    void init(){
        sieve();
        for (int n = 0; n < MAXN; n++){
            for (int m = 0; m < MAXM; m++){
                if (!n) dp[n][m] = m;
                else dp[n][m] = dp[n - 1][m] - dp[n - 1][m / primes[n - 1]];
            }
        }
    }

    uint64_t phi(uint64_t m, int n){
        if (n == 0) return m;
        if (primes[n - 1] >= m) return 1;
        if (m < MAXM && n < MAXN) return dp[n][m];
        return phi(m, n - 1) - phi(m / primes[n - 1], n - 1);
    }

    uint64_t lehmer(uint64_t m){
        if (m < MAX) return counter[m];

        uint64_t w, res = 0;
        int i, a, s, c, x, y;
        s = sqrt(0.9 + m), y = c = cbrt(0.9 + m);
        a = counter[y], res = phi(m, a) + a - 1;
        for (i = a; primes[i] <= s; i++) res = res - lehmer(m / primes[i]) + lehmer(primes[i]) - 1;
        return res;
    }
};

static void alloc(uint32_t N){
    np = new bool[N];
    prime = new uint32_t[N];
    pi = new uint32_t[N];
    memset(np, 0, sizeof(bool) * N);
}

static void dealloc(){
    delete[] np;
    delete[] prime;
    delete[] pi;
}

static void init(){
    egypt_sieve(N, np, prime);
    pi[0] = pi[1] = 0;
    for (uint32_t i = 2; i < N; i++)
        pi[i] = pi[i - 1] + !np[i];
    
    sz[0] = 1;
    for (uint32_t i = 0; i <= PM; i++) phi[i][0] = i;
    for (uint32_t i = 1; i <= M; i++) {
        sz[i] = prime[i] * sz[i - 1];
        for (uint32_t j = 1; j <= PM; j++)
            phi[j][i] = phi[j][i - 1] - phi[j / prime[i]][i - 1];
    }
}

static uint64_t getphi(uint64_t x, uint32_t s) {
    if (s == 0) return x;
    if (s <= M) return phi[x % sz[s]][s] + (x / sz[s]) * phi[sz[s]][s];
    if (x <= prime[s] * prime[s]) return pi[x] - s + 1;
    if (x <= prime[s] * prime[s] * prime[s] && x < N) {
        uint32_t s2x = pi[sqrt_ceil(x)];
        uint64_t ans  = pi[x] - (s2x + s - 2) * (s2x - s + 1) / 2;
        for (uint32_t i = s + 1; i <= s2x; ++i) ans += pi[x / prime[i]];
        return ans;
    }
    return getphi(x, s - 1) - getphi(x / prime[s], s - 1);
}

static uint64_t getpi(uint64_t x) {
    if (x < N) return pi[x];
    uint64_t ans = getphi(x, pi[cbrt_ceil(x)]) + pi[cbrt_ceil(x)] - 1;
    for (uint32_t i = pi[cbrt_ceil(x)] + 1, ed = pi[sqrt_ceil(x)]; i <= ed; ++i)
        ans -= getpi(x / prime[i]) - i + 1;
    return ans;
}

// 小于等于n的素数有多少个
static uint64_t lehmer_pi_calc(uint64_t x){
    if (x < N) return pi[x];
    uint32_t a = (uint32_t)lehmer_pi_calc(sqrt_ceil(sqrt_ceil(x)));
    uint32_t b = (uint32_t)lehmer_pi_calc(sqrt_ceil(x));
    uint32_t c = (uint32_t)lehmer_pi_calc(cbrt_ceil(x));
    uint64_t sum = getphi(x, a) + (uint64_t)(b + a - 2) * (b - a + 1) / 2;
    for (uint32_t i = a + 1; i <= b; i++) {
        uint64_t w = x / prime[i];
        sum -= lehmer_pi_calc(w);
        if (i > c) continue;
        uint64_t lim = lehmer_pi_calc(sqrt_ceil(w));
        for (uint32_t j = i; j <= lim; j++)
            sum -= lehmer_pi_calc(w / prime[j]) - (j - 1);
    }
    return sum;
}

uint64_t lehmer_pi(uint64_t x){
    uint32_t N = (uint32_t)cbrt_ceil(x) + 1;
    alloc(N);
    init();
    uint64_t ans = lehmer_pi_calc(x);
    dealloc();
    return ans;
}
