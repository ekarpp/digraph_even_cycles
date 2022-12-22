/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef GF_H
#define GF_H

#include <stdint.h>
#include <bitset>
#include <iostream>
#include <immintrin.h>
#include <set>

#include "global.hh"
#include "util.hh"

/* forward declare */
class GF_element;
class GR_element;

/* GF(2^n) */
class GF2_n
{
private:
    /* largest possible element in the field  + 1 */
    uint64_t mask;
    const int n;
    const uint64_t mod;

    uint64_t q_plus;
    uint64_t mod_ast;

    /* returns q s.t. for some r,
     * a = q*b + r is the division relation
     */
    uint64_t quo(uint64_t a, uint64_t b) const;

public:
    GF2_n(const int &e, const uint64_t &g);

    uint64_t ext_euclid(uint64_t a) const;

    /* carryless multiplication of a and b, polynomial multiplicatoin that is
     * done with Intel CLMUL
     */
    uint64_t clmul(uint64_t a, uint64_t b) const;

    virtual uint64_t rem(uint64_t a) const;
    // :)
    virtual __m256i wide_mul(__m256i a, __m256i b)
    { return _mm256_or_si256(a, b); }

    inline int get_n() const { return this->n; }
    inline uint64_t get_mod() const { return this->mod; }
    inline uint64_t get_mask() const { return this->mask; }
};

class GF2_16 : public GF2_n
{
public:
    using GF2_n::GF2_n;

    __m256i wide_mul(__m256i a, __m256i b);

    uint64_t rem(uint64_t a) const;
};

class GF2_32 : public GF2_n
{
public:
    using GF2_n::GF2_n;

    uint64_t rem(uint64_t a) const;
};

class GF_element
{
private:
    uint64_t repr;

public:
    GF_element() { }

    GF_element(const uint64_t n)
    {
        this->repr = n;
    }

    GF_element(const GF_element &e)
    {
        this->repr = e.get_repr();
    }

    inline GF_element operator+(const GF_element &other) const
    {
        return GF_element(this->repr ^ other.get_repr());
    }

    inline GF_element &operator+=(const GF_element &other)
    {
        this->repr ^= other.get_repr();
        return *this;
    }

    inline GF_element operator*(const GF_element &other) const
    {
        const uint64_t prod = global::F->clmul(
            this->repr,
            other.get_repr()
            );

        return GF_element(
            global::F->rem(prod)
            );
    }

    inline GF_element &operator*=(const GF_element &other)
    {
        const uint64_t prod = global::F->clmul(
            this->repr,
            other.get_repr()
            );

        this->repr = global::F->rem(prod);

        return *this;
    }

    inline GF_element inv() const
    {
        return GF_element(global::F->ext_euclid(this->repr));
    }

    inline void inv_in_place()
    {
        this->repr = global::F->ext_euclid(this->repr);
    }

    inline GF_element operator/(const GF_element &other) const
    {
        return *this * other.inv();
    }

    inline GF_element &operator/=(const GF_element &other)
    {
        const uint64_t inv = global::F->ext_euclid(other.get_repr());
        const uint64_t prod = global::F->clmul(
            this->repr,
            inv
            );

        this->repr = global::F->rem(prod);
        return *this;
    }

    inline bool operator==(const GF_element &other) const
    {
        return this->repr == other.get_repr();
    }

    inline uint64_t get_repr() const { return this->repr; }

    inline GF_element operator-(const GF_element &other) const
    {
        return *this + other;
    }

    inline GF_element &operator-=(const GF_element &other)
    {
        this->repr ^= other.get_repr();
        return *this;
    }

    inline GF_element &operator=(const GF_element &other)
    {
        this->repr = other.get_repr();
        return *this;
    }

    inline bool operator!=(const GF_element &other) const
    {
        return !(*this == other);
    }

    inline bool operator>(const GF_element &other) const
    {
        return this->repr > other.get_repr();
    }

    GR_element lift() const;

    void print() const
    {
        std::cout << std::bitset<8>(this->repr) << std::endl;
    }
};

namespace util
{
    std::vector<GF_element> distinct_elements(int n);

    inline GF_element GF_zero()
    {
        return GF_element(0);
    }

    inline GF_element GF_one()
    {
        return GF_element(1);
    }

    /* this can create zero, is it a problem? */
    inline GF_element GF_random()
    {
        return GF_element(global::randgen() & global::F->get_mask());
    }
}

#endif
