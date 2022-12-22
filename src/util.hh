/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef UTIL_H
#define UTIL_H

#include <algorithm>
#include <bitset>
#include <vector>

#include "global.hh"

namespace util
{
    inline int log2(uint64_t a)
    {
        return 63 - __builtin_clzl(a);
    }

    /* given an adjacency list for an undirected graph,
     * directs it such that edges are made one way
     * with direction chosen uniformly at random. */
    void direct_undirected(std::vector<std::vector<int>> &adj);

    uint64_t irred_poly(int deg);
    bool gcd1(int i, std::bitset<64> p);
}
#endif
