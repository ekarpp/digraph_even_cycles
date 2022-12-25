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
    const int n;
    const uint64_t mod;

    /* largest possible element in the field  + 1 */
    uint64_t mask;

    /* for rem in GF2_n */
    uint64_t q_plus;
    uint64_t mod_ast;

    /* returns q s.t. for some r,
     * a = q*b + r is the division relation
     */
    uint64_t quo(uint64_t a, const uint64_t b) const;

public:
    GF2_n(const int &e, const uint64_t &g);

    uint64_t ext_euclid(const uint64_t a) const;

    /* carryless multiplication of a and b, polynomial multiplicatoin that is
     * done with Intel CLMUL
     */
    inline uint64_t clmul(const uint64_t a, const uint64_t b) const
    {
        const __m128i prod = _mm_clmulepi64_si128(
            _mm_set_epi64x(0, a),
            _mm_set_epi64x(0, b),
            0x0
            );

        /* discard hi, only support up to 32 bit */
        return _mm_extract_epi64(prod, 0x0);
    }

    virtual uint64_t rem(const uint64_t a) const;

    /* multiply 8 bitsliced GF2_16 elements in one 256-bit vector.
     * define here for inlining and avoiding overhead from virtual func. */
    inline __m256i wide_mul(const __m256i &a, const __m256i &b) const
    {
        const __m128i prodlo = _mm_blend_epi32(
            _mm_shuffle_epi32(
                _mm_clmulepi64_si128(
                    _mm256_extractf128_si256(a, 0),
                    _mm256_extractf128_si256(b, 0),
                    0x11
                    ),
                0x8D
                ),
            _mm_shuffle_epi32(
                _mm_clmulepi64_si128(
                    _mm256_extractf128_si256(a, 0),
                    _mm256_extractf128_si256(b, 0),
                    0x00
                    ),
                0xD8
                ),
            0x3
            );

        const __m128i prodhi = _mm_blend_epi32(
            _mm_shuffle_epi32(
                _mm_clmulepi64_si128(
                    _mm256_extractf128_si256(a, 1),
                    _mm256_extractf128_si256(b, 1),
                    0x11
                    ),
                0x8D
                ),
            _mm_shuffle_epi32(
                _mm_clmulepi64_si128(
                    _mm256_extractf128_si256(a, 1),
                    _mm256_extractf128_si256(b, 1),
                    0x00
                    ),
                0xD8
                ),
            0x3
            );

        const __m256i prod = _mm256_set_m128i(prodhi, prodlo);

        const __m256i lomask = _mm256_set1_epi32(0xFFFF);

        const __m256i lo = _mm256_and_si256(
            prod,
            lomask
            );
        const __m256i hi = _mm256_srli_epi32(
            prod,
            16 // GF2_bits
            );

        const __m256i tmp = _mm256_xor_si256(
            hi,
            _mm256_xor_si256(
                _mm256_srli_epi16(hi, 14),
                _mm256_xor_si256(
                    _mm256_srli_epi16(hi, 13),
                    _mm256_srli_epi16(hi, 11)
                    )
                )
            );

        const __m256i rem_hi = _mm256_xor_si256(
            tmp,
            _mm256_xor_si256(
                _mm256_slli_epi16(tmp, 2),
                _mm256_xor_si256(
                    _mm256_slli_epi16(tmp, 3),
                    _mm256_slli_epi16(tmp, 5)
                    )
                )
            );

        return _mm256_xor_si256(rem_hi, lo);
    }

    inline int get_n() const { return this->n; }
    inline uint64_t get_mod() const { return this->mod; }
    inline uint64_t get_mask() const { return this->mask; }
};

class GF2_16 : public GF2_n
{
public:
    using GF2_n::GF2_n;

    uint64_t rem(const uint64_t a) const override;
};

class GF2_32 : public GF2_n
{
public:
    using GF2_n::GF2_n;

    uint64_t rem(const uint64_t a) const override;
};

class GF_element
{
private:
    uint64_t repr;

public:
    GF_element(): repr(0) { }

    explicit GF_element(const uint64_t n): repr(n) { }

    GF_element(const GF_element &e): repr(e.get_repr()) { }

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
        const uint64_t inv_repr = global::F->ext_euclid(other.get_repr());
        const uint64_t prod = global::F->clmul(
            this->repr,
            inv_repr
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
        std::cout << std::bitset<16>(this->repr) << " ";
    }
};

namespace util
{
    std::vector<GF_element> distinct_elements(const int n);

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
