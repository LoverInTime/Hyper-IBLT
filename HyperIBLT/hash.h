#ifndef HASH_H
#define HASH_H

#include <random>
#include <iostream>

#include <limits.h>
#include <stdint.h>
#include "murmur3.h"

namespace MurmurHash {
    template<typename T>
    inline uint32_t hash(const T& data, uint32_t seed = 0);
    inline uint32_t randomGenerator();

    static std::random_device rd;
    static std::mt19937 rng(rd());
    static std::uniform_real_distribution<double> dis(0, 1);

    inline uint32_t randomGenerator(){
        return rng();
    }

    template<typename T>
    inline uint32_t Hash(const T& data, uint32_t seed){
        uint32_t output;
        MurmurHash3_x86_32(&data, sizeof(T), seed, &output);
        return output;
    }
}

namespace RandomHash {
    int h(int i, int x);
}
namespace SpatialCoupling {
    int base_h0(int x);
    int h(int i, int x);
}

#endif
