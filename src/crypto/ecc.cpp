#include "ecc.h"

ECCInfo::ECCInfo(const BigInt& mod, const BigInt& a, const BigInt& b) :
    mod(mod), a(a), b(b), ec({a, mod}, {b, mod}) {
    g = ec.random();
}

ECCInfo::ECCInfo(const BigInt& mod, const BigInt& a, const BigInt& b, const Point2<TMod<BigInt>>& g) :
    mod(mod), a(a), b(b), g(g), ec({a, mod}, {b, mod}) {}

ECCInfo::~ECCInfo() {}

void ECCInfo::generateKeys(ECCPublicKey& pub, ECCPrivateKey& priv){
    priv.info = this;
    pub.info = this;

    BigInt key = randmod(mod);
    priv.key = key;
    pub.point = ec.mul(key, g);
}

ECCPublicKey::ECCPublicKey() : info(nullptr) {}

ECCPublicKey::~ECCPublicKey() {}

ECCCipherData ECCPublicKey::encrypt(ECCPlainData& data) const {
    data.M = info->ec.random();
    BigInt r = randmod(info->mod);
    Point2<TMod<BigInt>> S = info->ec.add(data.M, info->ec.mul(r, point));
    Point2<TMod<BigInt>> RG = info->ec.mul(r, info->g);
    return {RG, S};
}

bool ECCPublicKey::operator==(const ECCPublicKey& other) const {
    return info == other.info && point == other.point;
}

bool ECCPublicKey::operator!=(const ECCPublicKey& other) const {
    return !(*this == other);
}

ECCPrivateKey::ECCPrivateKey() : info(nullptr) {}

ECCPrivateKey::~ECCPrivateKey() {}

ECCPlainData ECCPrivateKey::decrypt(const ECCCipherData& message) const {
    Point2<TMod<BigInt>> M = info->ec.add(message.S, info->ec.neg(info->ec.mul(key, message.RG)));
    return {M};
}

bool ECCPrivateKey::operator==(const ECCPrivateKey& other) const {
    return info == other.info && key == other.key;
}

bool ECCPrivateKey::operator!=(const ECCPrivateKey& other) const {
    return !(*this == other);
}