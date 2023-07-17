#ifndef _MALLOC_H
#define _MALLOC_H
// 以后方便抽象

#include <malloc.h>

inline void* global_malloc(size_t _Size){
    return malloc(_Size);
}

inline void global_free(void* _Block){
    free(_Block);
}

inline void* global_realloc(void* _Block, size_t _Size){
    return realloc(_Block, _Size);
}

#endif // _MALLOC_H