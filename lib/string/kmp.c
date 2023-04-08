#include "kmp.h"

void kmp_init(const char* pattern, int* prefix, int nlen){
    prefix[0] = -1;
    for (int i = 1; i < nlen; i++){
        int j = prefix[i - 1];
        while (j >= 0 && pattern[j + 1] != pattern[i])
            j = prefix[j];
        prefix[i] = pattern[j + 1] == pattern[i] ? j + 1 : -1;
    }
}

int kmp_find(const char* text, const char* pattern, int* prefix, int nlen){
    int state = -1;
    return kmp_next(text, pattern, prefix, nlen, &state);
}

int kmp_next(const char* text, const char* pattern, int* prefix, int nlen, int* state){
    int j = *state;
    for (int i = 0; text[i]; i++){
        while (j >= 0 && pattern[j + 1] != text[i])
            j = prefix[j];
        if (pattern[j + 1] == text[i])
            j++;
        if (j == nlen - 1){
            *state = prefix[j];
            return i - nlen + 1;
        }
    }
    *state = j;
    return -1;
}