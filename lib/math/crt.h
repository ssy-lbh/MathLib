#ifndef _CRT_H
#define _CRT_H

#include <stdint.h>

uint64_t crt_calc(uint64_t remainders[], uint64_t mods[], uint32_t nlen);
uint64_t crt_calc(uint64_t remainders[], uint64_t mods[], uint32_t nlen, uint64_t mod);

#endif /* _CRT_H */