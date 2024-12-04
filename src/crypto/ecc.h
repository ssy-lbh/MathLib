#ifndef ECC_H
#define ECC_H

#include <cstdint>

#include "math/math_base.h"
#include "math/big_num.h"
#include "math/elliptic_curve.h"

class ECCPublicKey;
class ECCPrivateKey;

class ECCInfo {
public:
    const BigInt mod;
    const BigInt a;
    const BigInt b;

    Point2<TMod<BigInt>> g;

    const TEllipticCurve<TMod<BigInt>> ec;

public:
    ECCInfo(const BigInt& mod, const BigInt& a, const BigInt& b);
    ECCInfo(const BigInt& mod, const BigInt& a, const BigInt& b, const Point2<TMod<BigInt>>& g);
    ~ECCInfo();

    void generateKeys(ECCPublicKey& pub, ECCPrivateKey& priv);
};

struct ECCPlainData {
    Point2<TMod<BigInt>> M;
};

struct ECCCipherData {
    Point2<TMod<BigInt>> RG;
    Point2<TMod<BigInt>> S;
};

class ECCPublicKey {
private:
    const ECCInfo* info;
    Point2<TMod<BigInt>> point;

    friend class ECCInfo;

public:
    ECCPublicKey();
    ~ECCPublicKey();

    ECCCipherData encrypt(ECCPlainData& data) const;

    bool operator==(const ECCPublicKey& other) const;
    bool operator!=(const ECCPublicKey& other) const;
};

class ECCPrivateKey {
private:
    const ECCInfo* info;
    BigInt key;

    friend class ECCInfo;

public:
    ECCPrivateKey();
    ~ECCPrivateKey();

    ECCPlainData decrypt(const ECCCipherData& data) const;

    bool operator==(const ECCPrivateKey& other) const;
    bool operator!=(const ECCPrivateKey& other) const;
};

#endif // ECC_H