#include "lwe.h"

template <int N, int M>
LWEInfo<N, M>::LWEInfo(const BigInt& mod, const BigInt& B) : mod(mod), B(B) {}

template <int N, int M>
LWEInfo<N, M>::~LWEInfo() {}

template <int N, int M>
void LWEInfo<N, M>::generateKeys(LWEPublicKey<N, M>& pub, LWEPrivateKey<N, M>& priv){
    priv.info = this;
    pub.info = this;

    TMatrix<BigInt, M, 1> s;
    for(int i = 0; i < M; i++){
        s[i][0] = randmod(mod);
    }
    priv.s = s;

    TMatrix<BigInt, M, N> A;
    TMatrix<BigInt, N, 1> e;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            A[i][j] = randmod(mod);
        }
        e[i][0] = randmod(B);
    }
    pub.A = A;
    pub.b = (A * s + e) % mod;
}

template <int N, int M>
LWEPublicKey<N, M>::LWEPublicKey() : info(nullptr) {}

template <int N, int M>
LWEPublicKey<N, M>::~LWEPublicKey() {}

template <int N, int M>
LWECipherData<N, M> LWEPublicKey<N, M>::encrypt(LWEPlainData<N, M>& data) const {
    const BigInt& mod = info->mod;
    apply_ten(data.MSG, [](bool x){ return randmod(2u) == 1; });
    TMatrix<BigInt, N, 1> U = (A * data.M) % mod;
    BigInt V = (b * data.M)[0][0] % mod;
    return {U, V};
}

template <int N, int M>
bool LWEPublicKey<N, M>::operator==(const LWEPublicKey<N, M>& other) const {
    return info == other.info && A == other.A && b == other.b;
}

template <int N, int M>
bool LWEPublicKey<N, M>::operator!=(const LWEPublicKey<N, M>& other) const {
    return !(*this == other);
}

template <int N, int M>
LWEPrivateKey<N, M>::LWEPrivateKey() : info(nullptr) {}

template <int N, int M>
LWEPrivateKey<N, M>::~LWEPrivateKey() {}

template <int N, int M>
LWEPlainData<N, M> LWEPrivateKey<N, M>::decrypt(const LWECipherData<N, M>& data) const {
    const BigInt& mod = info->mod;
    TMatrix<BigInt, N, 1> e = (data.U * s - data.V) % mod;
    TMatrix<bool, M, 1> m;
    for(int i = 0; i < M; i++){
        m[i][0] = e[i][0] > mod / 2;
    }
    return {m};
}

template <int N, int M>
bool LWEPrivateKey<N, M>::operator==(const LWEPrivateKey<N, M>& other) const {
    return info == other.info && s == other.s;
}

template <int N, int M>
bool LWEPrivateKey<N, M>::operator!=(const LWEPrivateKey<N, M>& other) const {
    return !(*this == other);
}