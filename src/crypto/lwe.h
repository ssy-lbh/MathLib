#ifndef LWE_H
#define LWE_H

#include <cstdint>

#include "math/big_num.h"
#include "math/number_theory.h"
#include "math/matrix.h"

template <int N, int M>
class LWEPublicKey;
template <int N, int M>
class LWEPrivateKey;

// LWE: Learning With Errors
// https://en.wikipedia.org/wiki/Learning_with_errors
// @require: M * B < mod / 4
// @reference: (N, M, MOD, \chi) = (256, 640, 4093, 3.3)
template <int N, int M>
class LWEInfo {
public:
    const BigInt mod;
    const BigInt B;

public:
    LWEInfo(const BigInt& mod, const BigInt& B);
    ~LWEInfo();

    void generateKeys(LWEPublicKey<N, M>& pub, LWEPrivateKey<N, M>& priv);
};

template <int N, int M>
struct LWEPlainData {
    TMatrix<bool, M, 1> MSG;
};

template <int N, int M>
struct LWECipherData {
    TMatrix<BigInt, N, 1> U;
    BigInt V;
};

template <int N, int M>
class LWEPublicKey {
private:
    const LWEInfo<N, M>* info;
    const TMatrix<BigInt, M, N> A;
    const TMatrix<BigInt, M, 1> b;

    friend class LWEInfo<N, M>;

public:
    LWEPublicKey();
    ~LWEPublicKey();

    LWECipherData<N, M> encrypt(LWEPlainData<N, M>& data) const;

    bool operator==(const LWEPublicKey<N, M>& other) const;
    bool operator!=(const LWEPublicKey<N, M>& other) const;
};

template <int N, int M>
class LWEPrivateKey {
private:
    const LWEInfo<N, M>* info;
    TMatrix<BigInt, N, 1> s;

    friend class LWEInfo<N, M>;

public:
    LWEPrivateKey();
    ~LWEPrivateKey();

    LWEPlainData<N, M> decrypt(const LWECipherData<N, M>& data) const;

    bool operator==(const LWEPrivateKey<N, M>& other) const;
    bool operator!=(const LWEPrivateKey<N, M>& other) const;
};

#endif // LWE_H