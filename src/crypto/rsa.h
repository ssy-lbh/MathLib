#ifndef RSA_H
#define RSA_H

#include <cstdint>

#include "math/big_num.h"
#include "math/number_theory.h"

class RSAPublicKey;
class RSAPrivateKey;

class RSAInfo {
public:
    const BigInt mod;
    const BigInt p;
    const BigInt q;

public:
    RSAInfo(const BigInt& p, const BigInt& q);
    ~RSAInfo();

    void generateKeys(RSAPublicKey& pub, RSAPrivateKey& priv);
};

struct RSAPlainData {
    BigInt M;
};

struct RSACipherData {
    BigInt C;
};

class RSAPublicKey {
private:
    BigInt mod;
    BigInt key;

    friend class RSAInfo;

public:
    RSAPublicKey();
    ~RSAPublicKey();

    RSACipherData encrypt(RSAPlainData& data) const;

    bool operator==(const RSAPublicKey& other) const;
    bool operator!=(const RSAPublicKey& other) const;
};

class RSAPrivateKey {
private:
    BigInt mod;
    BigInt key;

    friend class RSAInfo;

public:
    RSAPrivateKey();
    ~RSAPrivateKey();

    RSAPlainData decrypt(const RSACipherData& message) const;

    bool operator==(const RSAPrivateKey& other) const;
    bool operator!=(const RSAPrivateKey& other) const;
};

#endif // RSA_H