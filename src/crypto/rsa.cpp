#include "rsa.h"

RSAInfo::RSAInfo(const BigInt& p, const BigInt& q) :
    mod(p * q), p(p), q(q) {}

RSAInfo::~RSAInfo() {}

void RSAInfo::generateKeys(RSAPublicKey& pub, RSAPrivateKey& priv){
    priv.mod = mod;
    pub.mod = mod;

    BigInt phi = (p - 1) * (q - 1);
    BigInt key = randmod(phi);
    while (gcd(key, phi) != 1){
        key = randmod(phi);
    }
    priv.key = key;
    pub.key = inv(key, phi);
}

RSAPublicKey::RSAPublicKey() : mod(0), key(0) {}

RSAPublicKey::~RSAPublicKey() {}

RSACipherData RSAPublicKey::encrypt(RSAPlainData& data) const {
    data.M = randmod(mod);
    BigInt C = pow(data.M, key, mod);
    return {C};
}

RSAPlainData RSAPublicKey::verify(const RSACipherData& message) const {
    BigInt M = pow(message.C, key, mod);
    return {M};
}

bool RSAPublicKey::operator==(const RSAPublicKey& other) const {
    return mod == other.mod && key == other.key;
}

bool RSAPublicKey::operator!=(const RSAPublicKey& other) const {
    return !(*this == other);
}

RSAPrivateKey::RSAPrivateKey() : mod(0), key(0) {}

RSAPrivateKey::~RSAPrivateKey() {}

RSAPlainData RSAPrivateKey::decrypt(const RSACipherData& message) const {
    BigInt M = pow(message.C, key, mod);
    return {M};
}

RSACipherData RSAPrivateKey::sign(const RSAPlainData& data) const {
    BigInt C = pow(data.M, key, mod);
    return {C};
}

bool RSAPrivateKey::operator==(const RSAPrivateKey& other) const {
    return mod == other.mod && key == other.key;
}

bool RSAPrivateKey::operator!=(const RSAPrivateKey& other) const {
    return !(*this == other);
}