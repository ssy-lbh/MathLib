
#ifndef SHA1_H
#define SHA1_H
/* ================ sha1.h ================ */
/*
SHA-1 in C
By Steve Reid <steve@edmweb.com>
100% Public Domain
*/

typedef struct {
    uint32_t state[5];
    uint32_t count[2];
    unsigned char buffer[64];
} SHA1_CTX;

void sha1_transform(uint32_t state[5], const unsigned char buffer[64]);
void sha1_init(SHA1_CTX* context);
void sha1_update(SHA1_CTX* context, const unsigned char* data, uint32_t len);
void sha1_final(unsigned char digest[20], SHA1_CTX* context);

#endif
