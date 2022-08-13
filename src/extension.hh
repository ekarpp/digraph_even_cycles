/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef EXTENSION_H
#define EXTENSION_H

#include <iostream>
#include <stdint.h>
#include <immintrin.h>

#include "gf.hh"
#include "util.hh"
#include "global.hh"
#include "bitvectors.hh"

/* forward declare */
class GF_element;
class Extension_element;

#define LONGLONG_ONES 0xffffFFFFffffFFFFull

/* representation for elements of E(4^n)
 * each bit in lo is the low bit of the mod 4 coefficient.
 * similarly for hi
 */
class extension_repr
{
public:
    __m128i repr;

    extension_repr() { }

    extension_repr(uint64_t hi, uint64_t lo)
    {
        this->repr = _mm_set_epi64x(hi, lo);
    }

    extension_repr(__m128i repr)
    {
        this->repr = repr;
    }

    void print() const
    {
        uint64_t lo = _mm_extract_epi64(this->repr, 0x0);
        uint64_t hi = _mm_extract_epi64(this->repr, 0x1);
        std::cout << hi << lo << std::endl;
    }

    extension_repr &operator>>=(int n)
    {
        this->repr = _mm_srli_epi64(this->repr, n);
        return *this;
    }

    extension_repr &operator<<=(int n)
    {
        this->repr = _mm_slli_epi64(this->repr, n);

        return *this;
    }

    extension_repr &operator&=(uint64_t m)
    {
        this->repr = _mm_and_si128(this->repr, _mm_set_epi64x(m, m));

        return *this;
    }

    extension_repr operator>>(int n) const
    {
        return _mm_srli_epi64(this->repr, n);
    }

    extension_repr operator<<(int n) const
    {
        return _mm_slli_epi64(this->repr, n);
    }

    extension_repr operator&(uint64_t m) const
    {
        return _mm_and_si128(this->repr, _mm_set_epi64x(m, m));
    }

    extension_repr shiftr_and(int n, uint64_t m) const
    {
        return _mm_and_si128(
            _mm_srli_epi64(this->repr, n),
            _mm_set_epi64x(m, m)
        );
    }

    bool is_even() const
    {
        char ZF = _mm_testz_si128(
            this->repr,
            _mm_set_epi64x(0x0, LONGLONG_ONES)
        );

        return ZF == 1;
    }

    extension_repr div2() const
    {
        return _mm_and_si128(this->repr, _mm_set_epi64x(LONGLONG_ONES, 0x0));
    }

    bool operator==(extension_repr other) const
    {
        __m128i cmp = _mm_cmpeq_epi64(this->repr, other.repr);
        int res = _mm_movemask_epi8(cmp);
        return res == 0xFFFF;
    }

    uint64_t get_lo() const
    {
        return _mm_extract_epi64(this->repr, 0x0);
    }

    uint64_t get_hi() const
    {
        return _mm_extract_epi64(this->repr, 0x1);
    }

    uint64_t log2() const
    {
        return std::max(util::log2(this->get_lo()), util::log2(this->get_hi()));
    }

    extension_repr mul_const(extension_repr c) const
    {
        __m128i cand =
            _mm_and_si128(
                this->repr,
                _mm_shuffle_epi32(
                    c.repr,
                    0x44
                )
            );
        __m128i cxor =
            _mm_and_si128(
                _mm_set_epi64x(LONGLONG_ONES, 0x0),
                _mm_and_si128(
                    c.repr,
                    _mm_shuffle_epi32(
                        this->repr,
                        0x44
                    )
                )
            );
        return _mm_xor_si128(cand, cxor);
    }

    extension_repr negate() const
    {
        return _mm_xor_si128(
            this->repr,
            _mm_and_si128(
                _mm_set_epi64x(LONGLONG_ONES, 0x0),
                _mm_shuffle_epi32(
                    this->repr,
                    0x44
                )
            )
        );
    }

    extension_repr add(extension_repr other) const
    {
        __m128i carry =
            _mm_and_si128(
                _mm_set_epi64x(LONGLONG_ONES, 0x0),
                _mm_shuffle_epi32(
                    _mm_and_si128(this->repr, other.repr),
                    0x44
                )
            );
        return _mm_xor_si128(
            this->repr,
            _mm_xor_si128(
                carry,
                other.repr
            )
        );
    }
};

#if GF2_bits == 16
typedef uint128_t kronecker_form;
#else
struct kronecker_form
{
/* MSB */
    uint256_t big;
/* 32bit LSB */
    uint64_t small;
};
#endif

/* Extension of GF(2^n) to the ring E(4^n).
 * If GF(2^n) = Z2 / <g_2> for irreducible polynomial
 * g_2 of degree n, then if g_4 is g_2 but coefficients
 * projected to Z4 we have E(4^n) = Z4 / <g_4>. */
class Extension
{
private:
    int n;
    uint64_t mod;
    uint64_t mask;

#if GF2_bits == 0
    std::vector<uint64_t> mod_ast;
    std::vector<uint64_t> q_plus;
#else
    extension_repr mod_ast;
    extension_repr q_plus;
#endif


    extension_repr n_prime;
    extension_repr r_squared;

    /* euclidean division, only used once during initialization.
     * b has to be monic for this to work */
    extension_repr quo(extension_repr a, extension_repr b) const
    {
        extension_repr q = { 0, 0 };
        int dega = a.log2();
        int degb = b.log2();

        while (dega >= degb)
        {

            extension_repr s = a & (1ll << dega);
            s >>= degb;

            q = this->add(q, s);
            a = this->subtract(a, this->fast_mul(s, b));

            dega = a.log2();
        }

        return q;
    }


#define DEGA                                                            \
    {                                                                   \
        if (a[1].get_hi() || a[1].get_lo())                             \
            dega = 64 + a[1].log2();                                    \
        else                                                            \
            dega = a[0].log2();                                         \
    }


    /* returns quotient form division of r_prime * x^(2n) by mod.
     * r_prime * x^(2n) has degree up to 3n-1 (n<32), hence this method. */
    extension_repr dumb_quo(extension_repr r_prime) const
    {
        extension_repr a[2];
        a[1] = r_prime >> (64 - 2*this->n);
        a[0] = r_prime << (2*this->n);

        extension_repr q = { 0, 0 };
        int degb = util::log2(this->mod);
        int dega;

        DEGA;
        while (dega >= degb)
        {
            int lohi = (dega >= 64) ? 1 : 0;
            uint64_t lo = a[lohi].get_lo() & (1ll << (dega - lohi*64));
            uint64_t hi = a[lohi].get_hi() & (1ll << (dega - lohi*64));
            extension_repr s = a[lohi] & (1ll << (dega - lohi*64));

            if (lohi)
                s <<= 64 - degb;
            else
                s >>= degb;

            q = this->add(q, s);

            /* prod of sb, computed with 2 bit multiplier.
             * same method as with the reference multiplication. */
            int degc = dega - degb;

            uint64_t tmp_hi = 0;
            uint64_t tmp_lo = 0;
            if (hi)
                tmp_hi = 0xFFFFFFFFFll;
            if (lo)
                tmp_lo = 0xFFFFFFFFFll;
            extension_repr tmp(tmp_hi, tmp_lo);

            extension_repr sg[2];

            /* tmp * g = sg[0] */
            sg[0] = tmp & this->mod;
            if (dega >= 64)
                sg[1] = sg[0] >> (64 - degc);
            else
                sg[1] = { 0, 0 };

            sg[0] <<= degc;

            a[1] = this->subtract(a[1], sg[1]);
            a[0] = this->subtract(a[0], sg[0]);

            DEGA;
        }

        return q;
    }

public:
    Extension() {}

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

        // N' = x^15 + x^14 + 3x^12 + x^7 + 3x^5 + 3x^4 + x^3 + x^2 + 3
        this->n_prime = { 0x1031, 0xD0BD };
        this->r_squared = { 0x018C, 0x0451 };
#elif GF2_bits == 32
        this->n = GF2_bits;

        /* x^32 + x^7 + x^3 + x^2 + 1 */
        this->mod = 0x10000008D;
        this->mask = 0xFFFFFFFF;

        // N' = x^30 + 3x^29 + x^28 + x^27 + x^26 + x^25 + x^23 + 3x^22 + x^21 + 2x^20 + 2x^19 + 3x^18 + x^15 + 2x^14 + 3x^13 + 2x^12 + x^10 + 3x^9 + 2x^8 + 2x^5 + 3x^4 + x^3 + x^2 + 3
        this->n_prime = { 0x205C7331, 0x7EE4A61D };
        this->r_squared = { 0x000006AC, 0x00004051 };
#else
        this->n = n;
        this->mod = mod;
        this->mask = (1ll << this->n) - 1;

        extension_repr q_plus_repr =
            this->quo({0, 1ull << (2*this->n)} , { 0, this->mod });
        for (int i = 0; i < this->n + 1; i++)
        {
            if (((this->mod >> i) & 1) && i < this->n)
                this->mod_ast.push_back(i);

            char qphi = (q_plus_repr.get_hi() >> i) & 1;
            char qplo = (q_plus_repr.get_lo() >> i) & 1;
            if (qphi || qplo)
            {
                this->q_plus.push_back(i);
                this->q_plus.push_back((qphi)
                                       ? 0xFFFFFFFFFFFFFFFFull
                                       : 0x0);
                this->q_plus.push_back((qplo)
                                       ? 0xFFFFFFFFFFFFFFFFull
                                       : 0x0);
            }
        }

        this->r_squared = {
            0,
            1ull << (this->n * 2)
        };
        this->r_squared = this->rem(this->r_squared);

        // deg == 2*n
        extension_repr r = { 0x0, 1ull << (this->n*2) };
        int N = 1 << this->n;
        N -= 1;
        N *= 2;

        // deg <= n-1
        extension_repr r_rem = this->rem(r);
        extension_repr r_prime = r_rem;

        N--;
        long idx = 1ll << (util::log2(N) - 1);
        while (idx > 1)
        {
            r_prime = this->rem(this->mul(r_prime, r_prime));
            if (N & idx)
                r_prime = this->rem(this->mul(r_prime, r_rem));
            idx >>= 1;
        }

        // deg <= 3n - 1, overflow when 3n - 1 > 64 <=> n > 65 / 4 ~ 21.667
        this->n_prime = this->dumb_quo(r_prime);
#endif

        if (global::output)
        {
            std::cout << "initialized E(4^" << this->n << ") with modulus: ";
            for (int i = this->n; i >= 0; i--)
            {
                if ((this->mod >> i) & 1)
                    std::cout << "1";
                else
                    std::cout << "0";
            }
            std::cout << std::endl;
        }
    }

    Extension_element zero() const;
    Extension_element one() const;
    Extension_element random() const;

    extension_repr rem(extension_repr a) const
    {
        return intel_rem(a);
    }

    extension_repr euclid_rem(extension_repr a) const
    {
        while (a.get_lo() > this->mask || a.get_hi() > this->mask)
        {
            int shift = a.log2();
            shift -= this->n;
            /* mod has coefficients modulo 2, thus its negation
             * is just it applied to hi and lo (see negate function)*/
            a = this->add(a, { this->mod << shift, this->mod << shift });
        }
        return a;
    }

    /* https://dl.acm.org/doi/10.1016/j.ipl.2010.04.011 */
    extension_repr intel_rem(extension_repr a) const
    {
        extension_repr hi = a >> this->n;
        extension_repr lo = a & this->mask;

#if GF2_bits == 16
        extension_repr tmp = hi >> 14;
        tmp = this->add(tmp, hi >> 13);
        tmp = this->add(tmp, hi >> 11);
        tmp = this->subtract(hi, tmp);

        extension_repr r = this->add(tmp, tmp << 2);
        r = this->add(r, tmp << 3);
        r = this->add(r, tmp << 5);
#elif GF2_bits == 32
        extension_repr tmp = hi >> 30;
        tmp = this->add(tmp, hi >> 29);
        tmp = this->add(tmp, hi >> 25);
        tmp = this->subtract(hi, tmp);

        extension_repr r = this->add(tmp, tmp << 2);
        r = this->add(r, tmp << 3);
        r = this->add(r, tmp << 7);
#else
        /* deg n-2 * deg n*/
        extension_repr rem = { 0x0, 0x0 };
        for (uint i = 0; i < this->q_plus.size() / 3; i++)
        {
            extension_repr tmp = this->mul_const(
                hi,
                { this->q_plus[3*i + 1], this->q_plus[3*i + 2] }
            );
            tmp <<= this->q_plus[3*i + 0];
            rem = this->add(rem, tmp);
        }
        rem >>= this->n;
        /* deg n-1 * deg n - 2*/
        extension_repr r = { 0x0, 0x0 };
        for (uint i = 0; i < this->mod_ast.size(); i++)
        {
            extension_repr tmp = rem << this->mod_ast[i];
            r = this->add(r, tmp);
        }
#endif
        r &= this->mask;
        return this->subtract(lo, r);
    }

#if GF2_bits == 16
    extension_repr packed_intel_rem(extension_repr a) const
    {
        uint64_t pack_mask = (this->mask | (this->mask << 32));

        extension_repr hi = a.shiftr_and(this->n, pack_mask);
        extension_repr lo = a & pack_mask;

        uint64_t m = 0b11 | (0b11ull << 32);
        extension_repr tmp = hi.shiftr_and(14, m);

        m = 0b111 | (0b111ull << 32);
        tmp = this->add(
            tmp,
            hi.shiftr_and(13, m)
        );
        m = 0b11111 | (0b11111ull << 32);
        tmp = this->add(
            tmp,
            hi.shiftr_and(11, m)
        );
        tmp = this->subtract(hi, tmp);

        extension_repr r = this->add(tmp, tmp << 2);
        r = this->add(r, tmp << 3);
        r = this->add(r, tmp << 5);

        r &= pack_mask;
        return this->subtract(lo, r);
    }

    extension_repr packed_fast_mul(extension_repr a, extension_repr b) const
    {
        __m128i ahi = _mm_srli_epi64(a.repr, 32);
        __m128i alo =
            _mm_and_si128(
                a.repr,
                _mm_set_epi64x(0xFFFF, 0xFFFF)
            );

        __m128i bhi = _mm_srli_epi64(b.repr, 32);
        __m128i blo =
            _mm_and_si128(
                b.repr,
                _mm_set_epi64x(0xFFFF, 0xFFFF)
            );

        __m128i prodlo =
            _mm_and_si128(
                _mm_set_epi64x(LONGLONG_ONES, 0x0),
                _mm_shuffle_epi32(
                    _mm_xor_si128(
                        _mm_clmulepi64_si128(alo, blo, 0x01),
                        _mm_clmulepi64_si128(alo, blo, 0x10)
                    ),
                    0x44
                )
            );

        __m128i prodhi =
            _mm_slli_epi64(
                _mm_and_si128(
                    _mm_set_epi64x(LONGLONG_ONES, 0x0),
                    _mm_shuffle_epi32(
                        _mm_xor_si128(
                            _mm_clmulepi64_si128(ahi, bhi, 0x01),
                            _mm_clmulepi64_si128(ahi, bhi, 0x10)
                        ),
                        0x44
                    )
                ),
                32
            );

        __m128i hilo = _mm_or_si128(prodhi, prodlo);

        uint64_t hi = 0;
        uint64_t lo = 0;

        uint64_t b_lo = b.get_lo();
        uint64_t a_lo = a.get_lo();
        /* handle product of lo and lo */
        #pragma GCC unroll 32
        for (int i = 0; i < GF2_bits; i++)
        {
            uint64_t msk = 0;
            if ((b_lo >> i) & 1)
                msk |= 0xFFFFFFFFull;
            if ((b_lo >> (i+32)) & 1)
                msk |= 0xFFFFFFFFull << 32;
            hi ^= (a_lo << i) & lo & msk;
            lo ^= (a_lo << i) & msk;
        }

        __m128i lolo = _mm_set_epi64x(hi, lo);
        return _mm_xor_si128(lolo, hilo);
    }
#endif

    extension_repr mont_rem(extension_repr a) const
    {
        /* d-1 deg * d-1 deg */
        extension_repr u = this->mul(a, this->n_prime);

        u &= this->mask;
        /* d deg * d-1 deg */
        extension_repr c = this->add(a, this->fast_mul(u, {0, this->mod}));

        c >>= this->n;
        return c;
    }

    extension_repr mont_form(extension_repr a) const
    {
        return this->mont_rem(this->mul(a, this->r_squared));
    }

    extension_repr mont_reduce(extension_repr a) const
    {
        return this->mont_rem(this->mul(a, { 0, 1 }));
    }

    extension_repr add(extension_repr a, extension_repr b) const
    {
        return a.add(b);
    }

    extension_repr negate(extension_repr a) const
    {
        return a.negate();
    }

    extension_repr subtract(extension_repr a, extension_repr b) const
    {
        return this->add(a, this->negate(b));
    }

    extension_repr mul(extension_repr a, extension_repr b) const
    {
        return kronecker_mul(a, b);
    }

    extension_repr mul_const(extension_repr a, extension_repr c) const
    {
        return a.mul_const(c);
    }

    extension_repr ref_mul(extension_repr a, extension_repr b) const
    {
        extension_repr c(0,0);
        uint64_t a_hi = a.get_hi();
        uint64_t a_lo = a.get_lo();

        for (int i = 0; i <= global::E.get_n(); i++)
        {
            uint64_t tmp_hi = 0;
            uint64_t tmp_lo = 0;
            if ((a_hi >> i) & 1)
                tmp_hi = this->mask;
            if ((a_lo >> i) & 1)
                tmp_lo = this->mask;

            extension_repr tmp(tmp_hi, tmp_lo);
            /* 2 bit carryless multiplier */
            extension_repr aib = this->mul_const(b, tmp);

            aib <<= i;
            c = global::E.add(c, aib);
        }

        return c;
    }

    extension_repr fast_mul(extension_repr a, extension_repr b) const
    {
        __m128i hilo = _mm_and_si128(
            _mm_set_epi64x(LONGLONG_ONES, 0x0),
            _mm_shuffle_epi32(
                _mm_xor_si128(
                    _mm_clmulepi64_si128(a.repr, b.repr, 0x01),
                    _mm_clmulepi64_si128(a.repr, b.repr, 0x10)
                ),
                0x44
            )
        );

        uint64_t hi = 0;
        uint64_t lo = 0;

        /* handle product of lo and lo */
        uint64_t b_lo = b.get_lo();
        uint64_t a_lo = a.get_lo();
        #pragma GCC unroll 32
#if GF2_bits == 16
        for (int i = 0; i <= GF2_bits; i++)
#else
        for (int i = 0; i <= 32; i++)
#endif
        {
            if ((b_lo >> i)&1)
            {
                hi ^= (a_lo << i) & lo;
                lo ^= (a_lo << i);
            }
        }
        __m128i lolo = _mm_set_epi64x(hi, lo);
        return _mm_xor_si128(lolo, hilo);
    }

    /* only works if deg <= 15 for a AND b */
    extension_repr kronecker_mul(extension_repr a, extension_repr b) const
    {
        /* we use different representation of polynomials than before here.
         * each bit string can be split to sets of 2 bits where each set
         * corresponds to a coefficient modulo 4. */
        kronecker_form aa = this->kronecker_substitution(a);
        kronecker_form bb = this->kronecker_substitution(b);
#if GF2_bits == 16
        uint256_t prod = bit::mul_128bit(aa, bb);

        /* first store the interesting bits to a uint64_t,
         * that is the first two bits of each 8 bit limb.
         * it fits, as we have deg <= 15+15 and each coefficient
         * uses two bits. */
        uint64_t extmask = 0x0303030303030303ull;
        uint64_t tmp = 0;
        for (int i = 0; i < 4; i++)
            tmp |= _pext_u64(prod.words[i], extmask) << (16*i);

        /* extract the usual hi/lo representation */
        uint64_t hiextmask = 0xAAAAAAAAAAAAAAAAull;
        uint64_t loextmask = 0x5555555555555555ull;
        extension_repr ret(
            _pext_u64(tmp, hiextmask),
            _pext_u64(tmp, loextmask)
        );
        return ret;
#else
        uint512_t ahbh = bit::mul_256bit(aa.big, bb.big);
        uint512_t ahbl = bit::mul_256bit_64bit(aa.big, bb.small);
        uint512_t albh = bit::mul_256bit_64bit(bb.big, aa.small);
        uint64_t albl = aa.small * bb.small;

        uint576_t prod = bit::add_576bit(
            bit::widen_512bits(ahbh),
            bit::pad_words(ahbl, 4)
        );

        prod = bit::add_576bit(
            prod,
            bit::pad_words(albh, 4)
        );

        /* can't overflow */
        prod.words[8] += albl;

        /* extract */

        /* append zero bit to MSB of each word
         * to make each word contain exactly 7 coefficients:
         * 7*9 + 1 = 63 + 1 = 64 */

        for (int i = 8; i > 0; i--)
        {
            prod.words[i] <<= i;
            prod.words[i] |= prod.words[i-1] >> (64 - i);
        }

        uint64_t tmp[3];
        tmp[0] = 0; tmp[1] = 0; tmp[2] = 0;

        uint64_t extmask = 0x00C06030180C0603ull;
        for (int i = 0; i < 4; i++)
            tmp[0] |= _pext_u64(prod.words[i], extmask) << (14*i);

        for (int i = 4; i < 8; i++)
            tmp[1] |= _pext_u64(prod.words[i], extmask) << (14*(i-4));

        tmp[2] = _pext_u64(prod.words[8], extmask);

        uint64_t hiextmask = 0xAAAAAAAAAAAAAAAAull;
        uint64_t loextmask = 0x5555555555555555ull;

        uint64_t ret_hi = 0;
        uint64_t ret_lo = 0;
        for (int i = 0; i < 3; i++)
        {
            ret_hi |= _pext_u64(tmp[i], hiextmask) << (28*i);
            ret_lo |= _pext_u64(tmp[i], loextmask) << (28*i);
        }
        extension_repr ret(ret_hi, ret_lo);
        return ret;
#endif
    }

    kronecker_form kronecker_substitution(extension_repr x) const
    {
        /* combine lo and hi to single uint64_t
         * where 2 bits represent single coefficient.
         * the "more traditional" bit representation for polynomials */
        uint64_t extmask = 0x5555555555555555ull;
        /* TODO */
        uint64_t comb = _pdep_u64(x.get_lo(), extmask);
        comb |= _pdep_u64(x.get_hi(), extmask << 1);

        kronecker_form kron;

#if GF2_bits == 16
        /* contains the "polynomial" after kronecker substitution.
         * for us it is sufficient that each coefficient has 8 bits,
         * (see details in thesis) thus we need 16*8 = 128 bits
         * for the polynomial after substitution. */

        extmask = 0x0303030303030303ull;
        kron.words[0] = _pdep_u64(comb & 0xFFFF, extmask);
        kron.words[1] = _pdep_u64((comb >> 16) & 0xFFFF, extmask);
#else
        /* each coefficients takes 9 bits.
         * we have <= 32 coefficients. */

        /* mask has 2x ones 7x zeros repeating */
        extmask = 0x00C06030180C0603ull;
        for (int i = 0; i < 4; i++)
            kron.big.words[i] = _pdep_u64((comb >> (i*14)) & 0x3FFF, extmask);

        /* remove the MSB zero from each word
         * to make the bitstring continuous */
        for (int i = 0; i < 3; i++)
        {
            kron.big.words[i] >>= i;
            kron.big.words[i] |= kron.big.words[i+1] << (63 - i);
        }
        kron.small = _pdep_u64((comb >> 56) & 0xFF, extmask);
        kron.big.words[3] >>= 3;
        kron.big.words[3] |= kron.small << 60;
        kron.small >>= 4;
#endif

        return kron;
    }

    int get_n() const { return this->n; }
    uint64_t get_mod() const { return this->mod; }
    uint64_t get_mask() const { return this->mask; }
};


class Extension_element
{
private:
    extension_repr repr;

public:
    Extension_element() { };

    Extension_element(const uint64_t lo, const uint64_t hi)
    {
        extension_repr e(hi, lo);
        this->repr = e;
    }

    Extension_element(const extension_repr repr)
    {
        this->repr = repr;
    }

    Extension_element(const Extension_element& e)
    {
        this->repr = e.get_repr();
    }

    Extension_element operator+(const Extension_element &other) const
    {
        return Extension_element(
            global::E.add(this->repr, other.get_repr())
        );
    }

    Extension_element &operator+=(const Extension_element &other)
    {
        this->repr = global::E.add(this->repr, other.get_repr());
        return *this;
    }

    Extension_element operator-(const Extension_element &other) const
    {
        return Extension_element(
            global::E.subtract(this->repr, other.get_repr())
        );
    }

    Extension_element &operator-=(const Extension_element &other)
    {
        this->repr = global::E.subtract(this->repr, other.get_repr());
        return *this;
    }

    Extension_element operator*(const Extension_element &other) const
    {
        extension_repr prod = global::E.mul(this->repr, other.get_repr());
        return Extension_element(global::E.rem(prod));
    }

    Extension_element &operator*=(const Extension_element &other)
    {
        extension_repr prod = global::E.mul(this->repr, other.get_repr());
        this->repr = global::E.rem(prod);
        return *this;
    }

    bool operator==(const Extension_element &other) const
    {
        return this->repr == other.get_repr();
    }

    bool is_even() const
    {
        return this->repr.is_even();
    }

    /* used only on elements that are multiplied by two
     * thus we can just move the hi to low
     * maybe even just return gf element straight away
     * as this (probably?) gets anyways done after div2 */
    /* modify instead of returning new? */
    Extension_element div2() const
    {
        return Extension_element(this->repr.div2());
    }

    GF_element project() const;

    extension_repr get_repr() const { return this->repr; }

    Extension_element &operator=(const Extension_element &other)
    {
        this->repr = other.get_repr();
        return *this;
    }

    bool operator!=(const Extension_element &other) const
    {
        return !(*this == other);
    }

    void print() const
    {
        this->repr.print();
    }
};

namespace util
{
    inline Extension_element tau(Extension_element sigma, Extension_element v)
    {
        return sigma.project().inv().lift() * v;
    }
}

#endif
