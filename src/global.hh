/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef GLOBAL_H
#define GLOBAL_H

#include "xorshift.hh"

class GF2_n;
class GR4_n;

namespace util
{
    class rand64bit
    {
    private:
        Xorshift gen;
    public:
        rand64bit() {}
        void init(uint64_t seed) { this->gen.init(seed); }
        uint64_t operator()() { return this->gen.next(); }
    };
}

namespace global
{
    /* these are defined in main.cc. no writing, just reading. */
    extern util::rand64bit randgen;
    extern GF2_n *F;
    extern GR4_n *E;
    extern bool output;
}


#endif
