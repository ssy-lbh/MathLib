#include "number_theory.h"

// �����ѵõ������ĳ˻�
static uint64_t totient_step(uint64_t x, uint64_t& n){
    // ���зֽ�
    uint64_t r = pollard_rho(x);
    // x Ϊ����, ��������������
    if (r == x){
        n = n / r * (r - 1);
        return r;
    }
    // x ��������, r1����Ϊ����, �ݹ鴦��
    uint64_t s = totient_step(r, n);
    uint64_t g = gcd(x, s);
    while (g > 1) {
        x /= g;
        g = gcd(x, s);
    }
    if (x == 1)
        return s;
    uint64_t t = totient_step(x, n);
    return s * t;
}

uint64_t totient(uint64_t n){
    if (n == 1)
        return 1;
    totient_step(n, n);
    return n;
}
