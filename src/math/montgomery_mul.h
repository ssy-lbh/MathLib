#ifndef MONTGOMERY_MUL_H
#define MONTGOMERY_MUL_H

#include "math_base.h"
#include "number_theory.h"

#include <cstdint>

template <typename T>
class MontgomeryMul {
public:
    T N, R;
    T N_inv, N_inv_neg, R2;
    T R_mask;
    uint64_t logR;

    /// @brief N must be odd
    /// @param N modulus
    MontgomeryMul(const T& N) : N(N) {
        logR = binsize(N);
        R = T(1) << logR;
        N_inv = inv(N, R);
        assert(N_inv != -1 && "modular inverse does not exist");
        N_inv_neg = R - N_inv;
        R2 = (T(1) << (logR << 1)) % N;
        R_mask = R - T(1);
    }

    MontgomeryMul(const T& N, uint64_t logR) : N(N), logR(logR) {
        R = T(1) << logR;
        N_inv = inv(N, R);
        assert(N_inv != -1 && "modular inverse does not exist");
        N_inv_neg = R - N_inv;
        R2 = (T(1) << (logR << 1)) % N;
        R_mask = R - T(1);
    }

    T reduce(const T& X) const {
        T m = ((X & R_mask) * N_inv_neg) & R_mask; // m = (X % R * N_inv_neg) % R        
        T t = (X + m * N) >> logR; // t = (T + m * N) / R
        return (t >= N ? t - N : t);
    }

    T modmul(const T& a, const T& b) const {
        assert(a < N && b < N && "input integer must be smaller than the modulus N");
        T aR = reduce(a * R2); // convert a to Montgomery form
        T bR = reduce(b * R2); // convert b to Montgomery form
        T X = aR * bR; // standard multiplication
        T abR = reduce(X); // Montgomery reduction
        return reduce(abR); // covnert abR to normal ab
    }
};

#endif