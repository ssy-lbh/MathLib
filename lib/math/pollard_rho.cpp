#include "number_theory.h"

uint64_t pollard_rho(uint64_t n){
    if (n == 4)
        return 2;
    if (miller_rabin(n))
        return n;
    while (1){
        uint64_t c = rand(1, n - 1); // ���������c
        auto f = [=](uint64_t x) {
            uint64_t y = (uint64_t)mul(x, x, n) + c;
            return (y >= n ? y - n : y);
        }; // lll��ʾ__int128�������
        int64_t t = 0, r = 0, p = 1, q;
        do{
            for (int i = 0; i < 128; i++){ // ��̶�����C=128
                t = f(t), r = f(f(r));
                if (t == r || (q = mul(p, abs(t - r), n)) == 0) // ������ֻ������߻�����Ϊ0���˳�
                    break;
                p = q;
            }
            uint64_t d = gcd(p, n);
            if (d > 1)
                return d;
        } while (t != r);
    }
}