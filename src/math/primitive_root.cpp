#include "number_theory.h"
#include "big_num.h"

bool check_root(uint64_t x, uint64_t p, uint64_t pm1_prime[], uint32_t cnt){
    uint64_t phi = p - 1;
    for (uint32_t i = 0; i < cnt; i++){
        if (pow(x, phi / pm1_prime[i], p) == 1)
            return false;
    }
    return true;
}

bool check_root(BigInt x, BigInt p, BigInt pm1_prime[], uint64_t cnt){
    BigInt phi = p - 1;
    for (uint64_t i = 0; i < cnt; i++){
        if (pow(x, phi / pm1_prime[i], p) == 1)
            return false;
    }
    return true;
}

uint64_t find_root(uint64_t p, uint64_t pm1_prime[], uint32_t cnt){
    for (uint64_t i = 2; i < p - 1; i++){
        if (check_root(i, p, pm1_prime, cnt))
            return i;
    }
    if (p == 1 || p == 2)
        return 1;
    if (p == 3)
        return 2;
    if (p == 4)
        return 3;
    return 0;
}

BigInt find_root(BigInt p, BigInt pm1_prime[], uint64_t cnt){
    for (BigInt i = 2; i < p - 1; i++){
        if (check_root(i, p, pm1_prime, cnt))
            return i;
    }
    if (p == 1 || p == 2)
        return 1;
    if (p == 3)
        return 2;
    if (p == 4)
        return 3;
    return 0;
}

// 1, 2, 4, p^k, 2p^k
bool has_root(uint64_t pm1_prime[], uint32_t exp[], uint32_t cnt){
    if (cnt > 2)
        return false;
    if (cnt == 2){
        if (pm1_prime[0] == 2){
            if (exp[0] != 1)
                return false;
        }
        else if (pm1_prime[1] == 2){
            if (exp[1] != 1)
                return false;
        }
        else
            return false;
        return true;
    }
    if (cnt == 1){
        if (pm1_prime[0] == 2)
            return exp[0] <= 2;
        return true;
    }
    return true;
}

bool has_root(BigInt pm1_prime[], uint64_t exp[], uint64_t cnt){
    if (cnt > 2)
        return false;
    if (cnt == 2){
        if (pm1_prime[0] == 2){
            if (exp[0] != 1)
                return false;
        }
        else if (pm1_prime[1] == 2){
            if (exp[1] != 1)
                return false;
        }
        else
            return false;
        return true;
    }
    if (cnt == 1){
        if (pm1_prime[0] == 2)
            return exp[0] <= 2;
        return true;
    }
    return true;
}