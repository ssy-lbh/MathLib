#ifndef _KMP_H
#define _KMP_H

void kmp_init(const char* pattern, int* prefix, int nlen);
int kmp_find(const char* text, const char* pattern, int* prefix, int nlen);
// 下次开始位置 text += nlen
int kmp_next(const char* text, const char* pattern, int* prefix, int nlen, int* state);

#endif /* _KMP_H */
