#include "math_base.h"
#include "number_theory.h"
#include "big_num.h"

static const uint32_t t_list[] = {
    2,
    12,
    60,
    180,
    840,
    1260,
    1680,
    2520,
    5040,
    15120,
    55440,
    110880,
    720720,
    1441440,
    4324320,
    24504480,
    73513440
};

bool apr_cl(BigInt x){
    if (x <= 3)
        return x >= 2;
    
    return true;
}