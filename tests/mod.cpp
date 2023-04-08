#include <cstdlib>

#include "math/number_theory.h"

void mod_base(){
    for (int i = 0; i < 1000; i++){
        unsigned int a = rand();
        assert(TMod<19>(a) == a % 19);
        assert(TMod<12289>(a) == a % 12289);
        assert(TMod<7340033>(a) == a % 7340033);
        assert(TMod<998244353>(a) == a % 998244353);
    }
}

void mod_inv(){
    for (int i = 0; i < 1000; i++){
        unsigned int a = rand();
        assert(a % 19 == 0 || inv(TMod<19>(a)) * a == 1);
        assert(a % 7681 == 0 || inv(TMod<7681>(a)) * a == 1);
        assert(a % 12289 == 0 || inv(TMod<12289>(a)) * a == 1);
        assert(a % 7340033 == 0 || inv(TMod<7340033>(a)) * a == 1);
        assert(a % 998244353 == 0 || inv(TMod<998244353>(a)) * a == 1);
    }
}

void mod_sqrt(){
    for (int i = 0; i < 1000; i++){
        unsigned int a = rand();
        TMod<19> r1 = sqrt(TMod<19>(a));
        assert((a % 19 == 0 && r1 == 0) || isnan(r1) || r1 * r1 == a % 19);
        TMod<7681> r2 = sqrt(TMod<7681>(a));
        assert((a % 7681 == 0 && r2 == 0) || isnan(r2) || r2 * r2 == a % 7681);
    }
}

void mod_fibonacci(){
    TMod<19> a = 5, rt;
    int n = 6;
    rt = sqrt(a);
    assert(!isnan(rt));
    assert((pow((rt + 1) / 2, n) - pow((-rt + 1) / 2, n)) / rt == 8);
}

int main(){
    mod_base();
    mod_inv();
    mod_sqrt();
    mod_fibonacci();
    return 0;
}