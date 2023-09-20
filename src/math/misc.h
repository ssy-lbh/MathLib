#ifndef MISC_H
#define MISC_H

#include "math_base.h"
#include "big_num.h"

BigFrac sqrt_approx1(BigInt x, uint64_t qsize = 64);
BigFrac sqrt_approx1(BigFrac x, uint64_t qsize = 64);
BigFrac cbrt_approx1(BigInt x, uint64_t qsize = 64);
BigFrac cbrt_approx1(BigFrac x, uint64_t qsize = 64);
BigFloat pi_borwein3(mpfr_prec_t prec = 64);
BigFloat pi_borwein4(mpfr_prec_t prec = 64);
BigFloat pi_borwein9(mpfr_prec_t prec = 64);
BigFloat pi_agm(mpfr_prec_t prec = 64);
// 未测试
uint32_t pi_bbp16(uint64_t n);
BigFrac pell_equation1(BigInt n);

#endif