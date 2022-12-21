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

/* representation for elements of E(4^n)
 * each bit in lo is the low bit of the mod 4 coefficient.
 * similarly for hi
 */
class extension_repr
{
public:
    uint64_t hi;
    uint64_t lo;

    extension_repr &operator>>=(int n)
    {
        this->hi >>= n;
        this->lo >>= n;

        return *this;
    }

    extension_repr &operator<<=(int n)
    {
        this->hi <<= n;
        this->lo <<= n;

        return *this;
    }

    extension_repr &operator&=(uint64_t m)
    {
        this->hi &= m;
        this->lo &= m;

        return *this;
    }

    extension_repr operator>>(int n)
    {
        return { this->hi >> n, this->lo >> n };
    }

    extension_repr operator<<(int n)
    {
        return { this->hi << n, this->lo << n };
    }

    extension_repr operator&(uint64_t m)
    {
        return { this->hi & m, this->lo & m };
    }

    extension_repr shiftr_and(int n, uint64_t m)
    {
        return {
            (this->hi >> n) & m,
            (this->lo >> n) & m
        };
    }

};

struct kronecker_form
{
    /* MSB */
    uint256_t big;
    /* 32bit LSB */
    uint64_t small;
    /* for E(4^16) */
    uint128_t b16;
};

/* Extension of GF(2^n) to the ring E(4^n).
 * If GF(2^n) = Z2 / <g_2> for irreducible polynomial
 * g_2 of degree n, then if g_4 is g_2 but coefficients
 * projected to Z4 we have E(4^n) = Z4 / <g_4>. */
class Extension
{
private:
    const int n;
    const uint64_t mod;
    uint64_t mask;

    std::vector<uint64_t> mod_ast;
    std::vector<uint64_t> q_plus;

    extension_repr n_prime;
    extension_repr r_squared;

    /* euclidean division, only used once during initialization.
     * b has to be monic for this to work */
    extension_repr quo(extension_repr a, extension_repr b) const
    {
        extension_repr q = { 0, 0 };
        int dega = std::max(util::log2(a.lo), util::log2(a.hi));
        int degb = std::max(util::log2(b.lo), util::log2(b.hi));

        while (dega >= degb)
        {

            extension_repr s = {
                (a.hi & (1ll << dega)) >> degb,
                (a.lo & (1ll << dega)) >> degb
            };

            q = this->add(q, s);
            a = this->subtract(a, this->fast_mul(s, b));

            dega = std::max(util::log2(a.lo), util::log2(a.hi));
        }

        return q;
    }


#define DEGA                                                              \
{                                                                         \
    if (a[1].hi || a[1].lo)                                               \
        dega = 64 + std::max(util::log2(a[1].hi), util::log2(a[1].lo));   \
    else                                                                  \
        dega = std::max(util::log2(a[0].hi), util::log2(a[0].lo));        \
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
            uint64_t lo = a[lohi].lo & (1ll << (dega - lohi*64));
            uint64_t hi = a[lohi].hi & (1ll << (dega - lohi*64));

            extension_repr s;
            if (lohi)
            {
                s = {
                    hi << (64 - degb),
                    lo << (64 - degb)
                };
            }
            else
            {
                s = {
                    hi >> degb,
                    lo >> degb
                };
            }
            q = this->add(q, s);

            /* prod of sb, computed with 2 bit multiplier.
             * same method as with the reference multiplication. */
            int degc = dega - degb;
            extension_repr tmp = { 0, 0 };
            if (hi)
                tmp.hi = 0xFFFFFFFFFll;
            if (lo)
                tmp.lo = 0xFFFFFFFFFll;
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

    void init_varying_size()
    {
        /* "intel rem" distributive law optimization */
        extension_repr q_plus_repr =
            this->quo({0, 1ull << (2*this->n)} , { 0, this->mod });
        for (int i = 0; i < this->n + 1; i++)
        {
            if (((this->mod >> i) & 1) && i < this->n)
                this->mod_ast.push_back(i);

            char qphi = (q_plus_repr.hi >> i) & 1;
            char qplo = (q_plus_repr.lo >> i) & 1;
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

        /* montgomery multiplication */
        this->r_squared = {
            0,
            1ull << (this->n * 2)
        };
        this->r_squared = this->rem(this->r_squared);

        // deg == n
        extension_repr r = { 0x0, 1ull << this->n };
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

        // deg <= 2n - 1
        this->n_prime = this->quo(this->mul(r_prime, r), { 0x0, this->mod });
    }

public:
    Extension() {}

    Extension(const int e, const uin64_t g): n(e), mod(g)
    {
        this->mask = (1ll << this->n) - 1;

        switch (n)
        {
        case 16:
            this->n_prime = { 0x1031, 0xD0BD };
            this->r_squared = { 0x018C, 0x0451 };
            break;
        case 32:
            this->n_prime = { 0x205C7331, 0x7EE4A61D };
            this->r_squared = { 0x000006AC, 0x00004051 };
            break;
        default:
            init_varying_size();
            break;
        }

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
        while (a.lo > this->mask || a.hi > this->mask)
        {
            int shift = std::max(util::log2(a.lo), util::log2(a.hi));
            shift -= this->n;
            /* mod has coefficients modulo 2, thus its negation
             * is just it applied to hi and lo (see negate function)*/
            a = this->add(a, { this->mod << shift, this->mod << shift });
        }
        return a;
    }

    extension_repr intel_rem(extension_repr a) const
    {
        return intel_rem<this->n>(a);
    }

    /* https://dl.acm.org/doi/10.1016/j.ipl.2010.04.011 */
    template <int exp>
    extension_repr intel_rem(extension_repr a) const
    {
        extension_repr hi = { a.hi >> this->n, a.lo >> this->n };
        extension_repr lo = { a.hi & this->mask, a.lo & this->mask };

        switch (exp)
        {
        case 16:
            extension_repr tmp = hi >> 14;
            tmp = this->add(tmp, hi >> 13);
            tmp = this->add(tmp, hi >> 11);
            tmp = this->subtract(hi, tmp);

            extension_repr r = this->add(tmp, tmp << 2);
            r = this->add(r, tmp << 3);
            r = this->add(r, tmp << 5);

            r &= this->mask;
            return this->subtract(lo, r);
        case 32:
            extension_repr tmp = hi >> 30;
            tmp = this->add(tmp, hi >> 29);
            tmp = this->add(tmp, hi >> 25);
            tmp = this->subtract(hi, tmp);

            extension_repr r = this->add(tmp, tmp << 2);
            r = this->add(r, tmp << 3);
            r = this->add(r, tmp << 7);

            r &= this->mask;
            return this->subtract(lo, r);
        default:
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

            r &= this->mask;
            return this->subtract(lo, r);
        }
    }

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
        uint64_t carry = a.lo & b.lo;
        return { carry ^ a.hi ^ b.hi, a.lo ^ b.lo };
    }

    extension_repr negate(extension_repr a) const
    {
        return {
            a.lo ^ a.hi,
            a.lo
        };
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
        return {
            (a.hi & c.lo) ^ (a.lo & c.hi),
            a.lo & c.lo
        };
    }

    extension_repr ref_mul(extension_repr a, extension_repr b) const
    {
        extension_repr c = { 0, 0 };

        for (int i = 0; i <= global::E.get_n(); i++)
        {
            extension_repr tmp = { 0, 0 };
            if ((a.hi >> i) & 1)
                tmp.hi = this->mask;
            if ((a.lo >> i) & 1)
                tmp.lo = this->mask;

            /* 2 bit carryless multiplier */
            extension_repr aib = this->mul_const(b, tmp);

            aib <<= i;
            c = global::E.add(c, aib);
        }

        return c;
    }

    extension_repr fast_mul(extension_repr a, extension_repr b) const
    {
        return fast_mul<this->n>(a, b);
    }

    template <int exp>
    extension_repr fast_mul(extension_repr a, extension_repr b) const
    {
        /* clean this up */
        __m128i aa = _mm_set_epi64x(a.hi, a.lo);
        __m128i bb = _mm_set_epi64x(b.hi, b.lo);

        __m128i alobhi = _mm_clmulepi64_si128(aa, bb, 0x01);
        __m128i ahiblo = _mm_clmulepi64_si128(aa, bb, 0x10);

        uint64_t hi1 = _mm_extract_epi64(ahiblo, 0x0);
        uint64_t hi2 = _mm_extract_epi64(alobhi, 0x0);

        uint64_t hi = 0;
        uint64_t lo = 0;

        /* handle product of lo and lo */
        #pragma GCC unroll 32
        for (int i = 0; i <= exp; i++)
        {
            if ((b.lo >> i)&1)
            {
                hi ^= (a.lo << i) & lo;
                lo ^= (a.lo << i);
            }
        }

        return { hi1 ^ hi2 ^ hi, lo };
    }

    extension_repr kronecker_mul(extension_repr a, extension_repr b) const
    {
        return kronecker_mul<this->n>(a, b);
    }

    /* only works if deg <= 15 for a AND b */
    template <int exp>
    extension_repr kronecker_mul(extension_repr a, extension_repr b) const
    {
        extension_repr ret;
        /* we use different representation of polynomials than before here.
         * each bit string can be split to sets of 2 bits where each set
         * corresponds to a coefficient modulo 4. */
        kronecker_form aa = this->kronecker_substitution(a);
        kronecker_form bb = this->kronecker_substitution(b);

        if (exp == 16)
        {
            uint256_t prod = bit::mul_128bit(aa.b16, bb.b16);

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
            ret.lo = _pext_u64(tmp, loextmask);
            ret.hi = _pext_u64(tmp, hiextmask);
            return ret;
        }

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
        ret.hi = 0; ret.lo = 0;
        for (int i = 0; i < 3; i++)
        {
            ret.hi |= _pext_u64(tmp[i], hiextmask) << (28*i);
            ret.lo |= _pext_u64(tmp[i], loextmask) << (28*i);
        }
        return ret;
    }

    kronecker_form kronecker_substitution(extension_repr x) const
    {
        return kronecker_substitution<this->n>(x);
    }

    template <int exp>
    kronecker_form kronecker_substitution(extension_repr x) const
    {
        /* combine lo and hi to single uint64_t
         * where 2 bits represent single coefficient.
         * the "more traditional" bit representation for polynomials */
        uint64_t extmask = 0x5555555555555555ull;
        uint64_t comb = _pdep_u64(x.lo, extmask);
        comb |= _pdep_u64(x.hi, extmask << 1);

        kronecker_form kron;

        if (exp == 16)
        {
            /* contains the "polynomial" after kronecker substitution.
             * for us it is sufficient that each coefficient has 8 bits,
             * (see details in thesis) thus we need 16*8 = 128 bits
             * for the polynomial after substitution. */
            extmask = 0x0303030303030303ull;
            kron.b16.words[0] = _pdep_u64(comb & 0xFFFF, extmask);
            kron.b16.words[1] = _pdep_u64(comb >> 16, extmask);
            return kron;
        }

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
        this->repr = { hi, lo };
    }

    Extension_element(const extension_repr repr)
    {
        this->repr = repr;
    }

    Extension_element(const Extension_element& e)
    {
        this->repr = { e.get_hi(), e.get_lo() };
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
        return this->repr.lo == other.get_lo() && this->repr.hi == other.get_hi();
    }

    bool is_even() const
    {
        return this->repr.lo == 0x0;
    }

    /* used only on elements that are multiplied by two
     * thus we can just move the hi to low
     * maybe even just return gf element straight away
     * as this (probably?) gets anyways done after div2 */
    /* modify instead of returning new? */
    Extension_element div2() const
    {
        return Extension_element(this->repr.hi, 0x0);
    }

    GF_element project() const;

    uint64_t get_lo() const { return this->repr.lo; }
    uint64_t get_hi() const { return this->repr.hi; }
    extension_repr get_repr() const { return this->repr; }

    Extension_element &operator=(const Extension_element &other)
    {
        this->repr.lo = other.get_lo();
        this->repr.hi = other.get_hi();
        return *this;
    }

    bool operator!=(const Extension_element &other) const
    {
        return !(*this == other);
    }

    void print() const
    {
        for (int i = 8; i >= 0; i--)
        {
            uint64_t v = (this->repr.hi >> i) & 1;
            v <<= 1;
            v |= (this->repr.lo >> i) & 1;
            std::cout << v;
        }
        std::cout << " ";
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
