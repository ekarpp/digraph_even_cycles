/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <stdint.h>
#include <iostream>

#include "gf.hh"
#include "extension.hh"
#include "global.hh"

using namespace std;

GF2_n::GF2_n(const int &e, const uint64_t &g): n(e), mod(g)
{
    this->mask = (1ll << this->n) - 1;
    if (n != 16 || n != 32)
    {
        this->q_plus = this->quo(1ull << (2*this->n), mod);
        this->mod_ast = this->mask & mod;
    }

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


/* returns q s.t. for some r,
 * a = q*b + r is the division relation
 */
uint64_t GF2_n::quo(uint64_t a, uint64_t b) const
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

/* returns s s.t. for some t: s*a + t*field.mod = gcd(field.mod, a)
 * <=> s*a + t*field.mod = 1 taking mod field.mod we get
 * s*a = 1 mod field.mod and thus a^-1 = s mod field.mod*/
uint64_t GF2_n::ext_euclid(uint64_t a) const
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
uint64_t GF2_n::clmul(uint64_t a, uint64_t b) const
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

__m256i GF2_16::wide_mul(__m256i a, __m256i b)
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

    const __m256i rem = _mm256_xor_si256(
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

/* returns r s.t. for some q,
 * a = q*field.mod + r is the division relation (in Z(2^n))
 */
uint64_t GF2_n::rem(uint64_t a) const
{
    uint64_t lo = a & this->mask;
    uint64_t hi = a >> this->n;

    uint64_t r = this->clmul(hi, this->q_plus);
    r >>= this->n;
    r = this->clmul(r, this->mod_ast);
    r &= this->mask;
    return r ^ lo;
}

uint64_t GF2_16::rem(uint64_t a) const
{
    uint64_t lo = a & 0xFFFF;
    uint64_t hi = a >> 16;

    uint64_t r = hi ^ (hi >> 14) ^ (hi >> 13) ^ (hi >> 11);
    r ^= (r << 2) ^ (r << 3) ^ (r << 5);
    r &= 0xFFFF;
    return r ^ lo;
}

uint64_t GF2_32::rem(uint64_t a) const
{
    uint64_t lo = a & 0xFFFFFFFF;
    uint64_t hi = a >> 32;

    uint64_t r = hi ^ (hi >> 30) ^ (hi >> 29) ^ (hi >> 25);
    r ^= (r << 2) ^ (r << 3) ^ (r << 7);
    r &= 0xFFFFFFFF;
    return r ^ lo;
}

GR_element GF_element::lift() const
{
    return GR_element(this->repr, 0b0);
}

namespace util
{
    /* returns n distinct random elements from
     * global::F-> (use LSFR?) */
    std::vector<GF_element> distinct_elements(int n)
    {
        std::vector<GF_element> vec(n);
        std::set<uint64_t> have;
        for (int i = 0; i < n; i++)
        {
            GF_element e = util::GF_random();
            while (have.count(e.get_repr()) == 1)
                e = util::GF_random();
            vec[i] = e;
            have.insert(e.get_repr());
        }
        return vec;
    }


}
