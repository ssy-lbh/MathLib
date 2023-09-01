#ifndef RADIX_SORT_H
#define RADIX_SORT_H

#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <assert.h>

template <typename T>
void radix_sort(T* start, T* end, T* temp, uint32_t key_size = sizeof(T)) {
    assert(key_size % 2 == 0 && "key_size must be even");
    uint32_t len = end - start;
    uint32_t bucket[1 << 8];
    uint32_t* bucket_end = bucket + (1 << 8);
    T* temp_end = temp + len;
    T* i;
    uint8_t* j;
    for (uint32_t k = 0; k < key_size; k += 2){
        memset(bucket, 0, sizeof(bucket));
        j = (uint8_t*)start + k; while (j < (uint8_t*)end) { bucket[*j]++; j += sizeof(T); }
        i = bucket_end - 2; while (i >= bucket) { *i += *(i+1); i--; }
        i = start; while (i < end) { temp[--bucket[*(uint8_t*)i + k]] = *i; i++; }
        memset(bucket, 0, sizeof(bucket));
        j = (uint8_t*)temp + k + 1; while (j < (uint8_t*)temp_end) { bucket[*j]++; j += sizeof(T); }
        i = bucket + 1; while (i < bucket_end) { *i += *(i-1); i++; }
        i = temp; while (i < temp_end) { start[--bucket[*((uint8_t*)i + k + 1)]] = *i; i++; }
    }
}

template <typename T>
void radix_sort(T* start, T* end, uint32_t key_size = sizeof(T)) {
    uint32_t len = end - start;
    T* temp = new T[len];
    radix_sort(start, end, temp, key_size);
    delete[] temp;
}

template <typename T>
void radix_sort2(T* start, T* end, T* temp, uint32_t key_size = sizeof(T)) {
    static_assert(sizeof(T) % 2 == 0, "sizeof(T) must be even");
    assert(key_size % 4 == 0 && "key_size must be multiple of 4");
    key_size >>= 1;
    uint32_t len = end - start;
    uint32_t bucket[1 << 16];
    uint32_t* bucket_end = bucket + (1 << 16);
    T* temp_end = temp + len;
    T* i;
    uint16_t* j;
    for (uint32_t k = 0; k < key_size; k += 2){
        memset(bucket, 0, sizeof(bucket));
        j = (uint16_t*)start + k; while (j < (uint16_t*)end) { bucket[*j]++; j += (sizeof(T) >> 1); }
        i = bucket_end - 2; while (i >= bucket) { *i += *(i+1); i--; }
        i = start; while (i < end) { temp[--bucket[*(uint16_t*)i + k]] = *i; i++; }
        memset(bucket, 0, sizeof(bucket));
        j = (uint16_t*)temp + k + 1; while (j < (uint16_t*)temp_end) { bucket[*j]++; j += (sizeof(T) >> 1); }
        i = bucket + 1; while (i < bucket_end) { *i += *(i-1); i++; }
        i = temp; while (i < temp_end) { start[--bucket[*((uint16_t*)i + k + 1)]] = *i; i++; }
    }
}

template <typename T>
void radix_sort2(T* start, T* end, uint32_t key_size = sizeof(T)) {
    uint32_t len = end - start;
    T* temp = new T[len];
    radix_sort2(start, end, temp, key_size);
    delete[] temp;
}

// @params keys => keys to sort, order from low to high
template <typename T>
void radix_sort(T* start, T* end, T* temp, uint32_t* keys, uint32_t num) {
    assert(key_size % 2 == 0 && "key_size must be even");
    uint32_t len = end - start;
    uint32_t bucket[1 << 8];
    uint32_t* bucket_end = bucket + (1 << 8);
    T* temp_end = temp + len;
    T* i;
    uint8_t* j;
    for (uint32_t k = 0; k < num; k++){
        uint32_t key = keys[k];
        if (k & 1){
            memset(bucket, 0, sizeof(bucket));
            j = (uint8_t*)temp + key; while (j < (uint8_t*)temp_end) { bucket[*j]++; j += sizeof(T); }
            i = bucket + 1; while (i < bucket_end) { *i += *(i-1); i++; }
            i = temp; while (i < temp_end) { start[--bucket[*((uint8_t*)i + key)]] = *i; i++; }
        } else {
            memset(bucket, 0, sizeof(bucket));
            j = (uint8_t*)start + key; while (j < (uint8_t*)end) { bucket[*j]++; j += sizeof(T); }
            i = bucket_end - 2; while (i >= bucket) { *i += *(i+1); i--; }
            i = start; while (i < end) { temp[--bucket[*(uint8_t*)i + key]] = *i; i++; }
        }
    }
    if (num & 1){
        memcpy(start, temp, sizeof(T) * len);
    }
}

// @params keys => keys to sort, order from low to high
template <typename T>
void radix_sort(T* start, T* end, uint32_t* keys, uint32_t num) {
    uint32_t len = end - start;
    T* temp = new T[len];
    radix_sort(start, end, temp, keys, num);
    delete[] temp;
}

// @param offset => offsetof(type, member)
template <typename T>
void float_to_rank32(T* start, T* end, uint32_t offset){
    while (start < end){
        int* p = (int*)((uint8_t*)start + offset);
        int x = *p >> 31;
        *p = (*p ^ x) | (~x << 31);
        start++;
    }
}

// @param offset => offsetof(type, member)
template <typename T>
void rank32_to_float(T* start, T* end, uint32_t offset){
    while (start < end){
        int* p = (int*)((uint8_t*)start + offset);
        int x = *p >> 31;
        *p = (*p ^ ~x) | (x << 31);
        start++;
    }
}

// @param offset => offsetof(type, member)
template <typename T>
void float_to_rank64(T* start, T* end, uint32_t offset){
    while (start < end){
        long long* p = (long long*)((uint8_t*)start + offset);
        long long x = *p >> 63;
        *p = (*p ^ x) | (~x << 63);
        start++;
    }
}

// @param offset => offsetof(type, member)
template <typename T>
void rank64_to_float(T* start, T* end, uint32_t offset){
    while (start < end){
        long long* p = (long long*)((uint8_t*)start + offset);
        long long x = *p >> 63;
        *p = (*p ^ ~x) | (x << 63);
        start++;
    }
}

// @param offset => offsetof(type, member)
template <typename T>
void int_to_rank32(T* start, T* end, uint32_t offset){
    while (start < end){
        int* p = (int*)((uint8_t*)start + offset);
        *p = *p ^ (1 << 31);
        start++;
    }
}

// @param offset => offsetof(type, member)
template <typename T>
void rank32_to_int(T* start, T* end, uint32_t offset){
    while (start < end){
        int* p = (int*)((uint8_t*)start + offset);
        *p = *p ^ (1 << 31);
        start++;
    }
}

// @param offset => offsetof(type, member)
template <typename T>
void int_to_rank64(T* start, T* end, uint32_t offset){
    while (start < end){
        long long* p = (long long*)((uint8_t*)start + offset);
        *p = *p ^ (1ll << 63);
        start++;
    }
}

// @param offset => offsetof(type, member)
template <typename T>
void rank64_to_int(T* start, T* end, uint32_t offset){
    while (start < end){
        long long* p = (long long*)((uint8_t*)start + offset);
        *p = *p ^ (1ll << 63);
        start++;
    }
}

#endif // RADIX_SORT_H