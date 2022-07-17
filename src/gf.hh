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
class Extension_element;

/* GF(2^n) */
class GF2n
{
private:
    /* largest possible element in the field */
    uint64_t mask;
    int n;
    uint64_t mod;

    uint64_t q_plus;
    uint64_t mod_ast;

    /* returns q s.t. for some r,
     * a = q*b + r is the division relation
     */
    uint64_t quo(uint64_t a, uint64_t b) const
    {
        uint64_t q = 0b0;
        int degb = util::log2(b);
        while (a >= b)
        {
            int shift = util::log2(a);
            shift -= degb;
            /* shift = deg(a) - deg(b) */
            q ^= (1ll << shift);
            a ^= (b << shift);
        }
        return q;
    }

public:
    GF2n() {}
#if GF2_bits == 0
    void init(const int n, const uint64_t mod)
#else
    void init()
#endif
    {
#if GF2_bits == 16
        this->n = GF2_bits;
        /* x^16 + x^5 + x^3 + x^2 +  1 */
        this->mod = 0x1002D;
        this->mask = 0xFFFF;
#elif GF2_bits == 32
        this->n = GF2_bits;
        /* x^32 + x^7 + x^3 + x^2 + 1 */
        this->mod = 0x10000008D;
        this->mask = 0xFFFFFFFF;
#else
        this->n = n;
        this->mod = mod;
        this->mask = (1ll << this->n) - 1;
        this->q_plus = this->quo(1ull << (2*this->n), mod);
        this->mod_ast = this->mask & mod;
#endif
        if (global::output)
        {
            std::cout << "initialized GF(2^" << this->n << ") with modulus: ";
            for (int i = n; i >= 0; i--)
            {
                if ((this->mod >> i) & 1)
                    std::cout << "1";
                else
                    std::cout << "0";
            }
            std::cout << std::endl;
        }
    }

    GF_element zero() const;
    GF_element one() const;
    GF_element random() const;

    /* returns r s.t. for some q,
     * a = q*field.mod + r is the division relation (in Z(2^n))
     */
    uint64_t rem(uint64_t a) const
    {
        uint64_t lo = a & this->mask;
        uint64_t hi = a >> this->n;

#if GF2_bits == 16
        uint64_t r = hi ^ (hi >> 14) ^ (hi >> 13) ^ (hi >> 11);
        r ^= (r << 2) ^ (r << 3) ^ (r << 5);
#elif GF2_bits == 32
        uint64_t r = hi ^ (hi >> 30) ^ (hi >> 29) ^ (hi >> 25);
        r ^= (r << 2) ^ (r << 3) ^ (r << 7);
#else
        uint64_t r = this->clmul(hi, this->q_plus);
        r >>= this->n;
        r = this->clmul(r, this->mod_ast);
#endif
        r &= this->mask;
        return r ^ lo;
    }

#if GF2_bits == 16
    __m256i wide_mul(__m256i a, __m256i b)
    {
        /* al/bl might not be needed, just use a/b */
        __m256i lomask = _mm256_set_epi64x(
            0xFFFF,
            0xFFFF,
            0xFFFF,
            0xFFFF
        );

        const __m256i al = _mm256_and_si256(
            a,
            lomask
        );
        const __m256i bl = _mm256_and_si256(
            b,
            lomask
        );
        const __m256i ah = _mm256_srli_epi64(
            a,
            32
        );
        const __m256i bh = _mm256_srli_epi64(
            b,
            32
        );

        const __m128i prodhi = _mm_or_si128(
            _mm_unpacklo_epi64(
                _mm_clmulepi64_si128(
                    _mm256_extractf128_si256(al, 1),
                    _mm256_extractf128_si256(bl, 1),
                    0x00
                    ),
                _mm_clmulepi64_si128(
                    _mm256_extractf128_si256(al, 1),
                    _mm256_extractf128_si256(bl, 1),
                    0x11
                    )
                ),
            _mm_slli_epi64(
                _mm_unpacklo_epi64(
                    _mm_clmulepi64_si128(
                        _mm256_extractf128_si256(ah, 1),
                        _mm256_extractf128_si256(bh, 1),
                        0x00
                        ),
                    _mm_clmulepi64_si128(
                        _mm256_extractf128_si256(ah, 1),
                        _mm256_extractf128_si256(bh, 1),
                        0x11
                        )
                    ),
                32
                )
            );

        const __m128i prodlo = _mm_or_si128(
            _mm_unpacklo_epi64(
                _mm_clmulepi64_si128(
                    _mm256_extractf128_si256(al, 0),
                    _mm256_extractf128_si256(bl, 0),
                    0x00
                    ),
                _mm_clmulepi64_si128(
                    _mm256_extractf128_si256(al, 0),
                    _mm256_extractf128_si256(bl, 0),
                    0x11
                    )
                ),
            _mm_slli_epi64(
                _mm_unpacklo_epi64(
                    _mm_clmulepi64_si128(
                        _mm256_extractf128_si256(ah, 0),
                        _mm256_extractf128_si256(bh, 0),
                        0x00
                        ),
                    _mm_clmulepi64_si128(
                        _mm256_extractf128_si256(ah, 0),
                        _mm256_extractf128_si256(bh, 0),
                        0x11
                        )
                    ),
                32
                )
            );

        __m256i prod = _mm256_set_m128i(prodhi, prodlo);

        lomask = _mm256_set_epi32(
            0xFFFF,
            0xFFFF,
            0xFFFF,
            0xFFFF,
            0xFFFF,
            0xFFFF,
            0xFFFF,
            0xFFFF
        );
        __m256i lo = _mm256_and_si256(
            prod,
            lomask
        );
        __m256i hi = _mm256_srli_epi32(
            prod,
            16 // GF2_bits
        );

        __m256i tmp = _mm256_xor_si256(
            hi,
            _mm256_xor_si256(
                _mm256_srli_epi16(hi, 14),
                _mm256_xor_si256(
                    _mm256_srli_epi16(hi, 13),
                    _mm256_srli_epi16(hi, 11)
                )
            )
        );

        __m256i rem = _mm256_xor_si256(
            tmp,
            _mm256_xor_si256(
                _mm256_slli_epi16(tmp, 2),
                _mm256_xor_si256(
                    _mm256_slli_epi16(tmp, 3),
                    _mm256_slli_epi16(tmp, 5)
                )
            )
        );

        return _mm256_xor_si256(rem, lo);
    }

    uint64_t packed_rem(uint64_t a) const
    {
        uint64_t pack_mask = this->mask | (this->mask << 32);
        uint64_t lo = a & pack_mask;
        uint64_t hi = (a >> this->n) & pack_mask;

        uint64_t r = hi;

        uint64_t m = 0b11 | (0b11ull << 32);
        r ^= (hi >> 14) & m;
        m = 0b111 | (0b111ull << 32);
        r ^= (hi >> 13) & m;
        m = 0b11111 | (0b11111ull << 32);
        r ^= (hi >> 11) & m;
        r &= pack_mask;
        r ^= (r << 2) ^ (r << 3) ^ (r << 5);
        r &= pack_mask;

        return r ^ lo;
    }
#endif

    /* returns s s.t. for some t: s*a + t*field.mod = gcd(field.mod, a)
     * <=> s*a + t*field.mod = 1 taking mod field.mod we get
     * s*a = 1 mod field.mod and thus a^-1 = s mod field.mod*/
    uint64_t ext_euclid(uint64_t a) const
    {
        // assert(a != 0)
        uint64_t s = 0x1;
        uint64_t s_next = 0x0;
        uint64_t r = a;
        uint64_t r_next = this->mod;
        uint64_t tmp;

        while (r_next != 0x0)
        {
            uint64_t q = this->quo(r, r_next);
            tmp = r ^ this->clmul(q, r_next);
            r = r_next;
            r_next = tmp;

            tmp = s ^ this->clmul(q, s_next);
            s = s_next;
            s_next = tmp;
        }

        return s;
    }

    /* carryless multiplication of a and b, polynomial multiplicatoin that is
     * done with Intel CLMUL
     */
    uint64_t clmul(uint64_t a, uint64_t b) const
    {
        const __m128i prod = _mm_clmulepi64_si128(
            _mm_set_epi64x(0, a),
            _mm_set_epi64x(0, b),
            0x0
        );

        uint64_t lo = _mm_extract_epi64(prod, 0x0);
        /* discard hi, only support up to 32 bit */
        return lo;
    }

#if GF2_bits == 16
    uint64_t packed_clmul(uint64_t a, uint64_t b) const
    {
        __m128i aa = _mm_set_epi64x(a >> 32, a & this->mask);
        __m128i bb = _mm_set_epi64x(b >> 32, b & this->mask);
        __m128i prod[2];
        prod[0] = _mm_clmulepi64_si128(aa, bb, 0x00);
        prod[1] = _mm_clmulepi64_si128(aa, bb, 0x11);

        uint64_t res = _mm_extract_epi64(prod[0], 0x0);
        res |= _mm_extract_epi64(prod[1], 0x0) << 32;

        return res;
    }
#endif

    int get_n() const { return this->n; }
    uint64_t get_mod() const { return this->mod; }
    uint64_t get_mask() const { return this->mask; }
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

    GF_element operator+(const GF_element &other) const
    {
        return GF_element(this->repr ^ other.get_repr());
    }

    GF_element &operator+=(const GF_element &other)
    {
        this->repr ^= other.get_repr();
        return *this;
    }

    GF_element operator*(const GF_element &other) const
    {
        const uint64_t prod = global::F.clmul(
            this->repr,
            other.get_repr()
        );

        return GF_element(
            global::F.rem(prod)
        );
    }

    GF_element &operator*=(const GF_element &other)
    {
        const uint64_t prod = global::F.clmul(
            this->repr,
            other.get_repr()
        );

        this->repr = global::F.rem(prod);

        return *this;
    }

    GF_element inv() const
    {
        return GF_element(global::F.ext_euclid(this->repr));
    }

    void inv_in_place()
    {
        this->repr = global::F.ext_euclid(this->repr);
    }

    GF_element operator/(const GF_element &other) const
    {
        return *this * other.inv();
    }

    GF_element &operator/=(const GF_element &other)
    {
        const uint64_t inv = global::F.ext_euclid(other.get_repr());
        const uint64_t prod = global::F.clmul(
            this->repr,
            inv
            );

        this->repr = global::F.rem(prod);
        return *this;
    }

    bool operator==(const GF_element &other) const
    {
        return this->repr == other.get_repr();
    }

    uint64_t get_repr() const { return this->repr; }

    GF_element operator-(const GF_element &other) const
    {
        return *this + other;
    }

    GF_element &operator-=(const GF_element &other)
    {
        this->repr ^= other.get_repr();
        return *this;
    }

    GF_element &operator=(const GF_element &other)
    {
        this->repr = other.get_repr();
        return *this;
    }

    Extension_element lift() const;

    bool operator!=(const GF_element &other) const
    {
        return !(*this == other);
    }

    bool operator>(const GF_element &other) const
    {
        return this->repr > other.get_repr();
    }

    void print() const
    {
        std::cout << std::bitset<8>(this->repr) << std::endl;
    }
};

namespace util
{
    /* returns n distinct random elements from
     * global::F. (use LSFR?) */
    inline std::vector<GF_element> distinct_elements(int n)
    {
        std::vector<GF_element> vec(n);
        std::set<uint64_t> have;
        for (int i = 0; i < n; i++)
        {
            GF_element e = global::F.random();
            while (have.count(e.get_repr()) == 1)
                e = global::F.random();
            vec[i] = e;
            have.insert(e.get_repr());
        }
        return vec;
    }
}

#endif
