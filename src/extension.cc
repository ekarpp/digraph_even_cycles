/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>

#include "gf.hh"
#include "global.hh"
#include "extension.hh"

using namespace std;

GR4_n::GR4_n(const int e, const uint64_t g): n(e), mod(g)
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

void GR4_n::init_varying_size()
{
    /* "intel rem" distributive law optimization */
    GR_repr q_plus_repr =
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
    GR_repr r = { 0x0, 1ull << this->n };
    int N = 1 << this->n;
    N -= 1;
    N *= 2;

    // deg <= n-1
    GR_repr r_rem = this->rem(r);
    GR_repr r_prime = r_rem;

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


/* euclidean division, only used once during initialization.
 * b has to be monic for this to work */
GR_repr GR4_n::quo(GR_repr a, GR_repr b) const
{
    GR_repr q = { 0, 0 };
    int dega = std::max(util::log2(a.lo), util::log2(a.hi));
    int degb = std::max(util::log2(b.lo), util::log2(b.hi));

    while (dega >= degb)
    {

        GR_repr s = {
            (a.hi & (1ll << dega)) >> degb,
            (a.lo & (1ll << dega)) >> degb
        };

        q = this->add(q, s);
        a = this->subtract(a, this->fast_mul(s, b));

        dega = std::max(util::log2(a.lo), util::log2(a.hi));
    }

    return q;
}

GR_repr GR4_n::ref_mul(GR_repr a, GR_repr b) const
{
    GR_repr c = { 0, 0 };

    for (int i = 0; i <= global::E->get_n(); i++)
    {
        GR_repr tmp = { 0, 0 };
        if ((a.hi >> i) & 1)
            tmp.hi = this->mask;
        if ((a.lo >> i) & 1)
            tmp.lo = this->mask;

        /* 2 bit carryless multiplier */
        GR_repr aib = this->mul_const(b, tmp);

        aib <<= i;
        c = global::E->add(c, aib);
    }

    return c;
}

GR_repr GR4_n::fast_mul(GR_repr a, GR_repr b) const
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
    for (int i = 0; i <= this->n; i++)
    {
        if ((b.lo >> i)&1)
        {
            hi ^= (a.lo << i) & lo;
            lo ^= (a.lo << i);
        }
    }

    return { hi1 ^ hi2 ^ hi, lo };
}

virtual GR_repr GR4_n::kronecker_mul(GR_repr a, GR_repr b) const
{
    /* we use different representation of polynomials than before here.
     * each bit string can be split to sets of 2 bits where each set
     * corresponds to a coefficient modulo 4. */
    kronecker_form aa = this->kronecker_substitution(a);
    kronecker_form bb = this->kronecker_substitution(b);

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
    GR_repr ret;
    ret.hi = 0; ret.lo = 0;
    for (int i = 0; i < 3; i++)
    {
        ret.hi |= _pext_u64(tmp[i], hiextmask) << (28*i);
        ret.lo |= _pext_u64(tmp[i], loextmask) << (28*i);
    }
    return ret;
}

virtual kronecker_form GR4_n::kronecker_substitution(GR_repr x) const
{
    /* combine lo and hi to single uint64_t
     * where 2 bits represent single coefficient.
     * the "more traditional" bit representation for polynomials */
    uint64_t extmask = 0x5555555555555555ull;
    uint64_t comb = _pdep_u64(x.lo, extmask);
    comb |= _pdep_u64(x.hi, extmask << 1);

    /* contains the "polynomial" after kronecker substitution. */
    kronecker_form kron;

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


kronecker_form GR4_16::kronecker_substitution(GR_repr x) const
{
    /* combine lo and hi to single uint64_t
     * where 2 bits represent single coefficient.
     * the "more traditional" bit representation for polynomials */
    uint64_t extmask = 0x5555555555555555ull;
    uint64_t comb = _pdep_u64(x.lo, extmask);
    comb |= _pdep_u64(x.hi, extmask << 1);

    /* contains the "polynomial" after kronecker substitution.
     * for us it is sufficient that each coefficient has 8 bits,
     * (see details in thesis) thus we need 16*8 = 128 bits
     * for the polynomial after substitution. */
    kronecker_form kron;
    extmask = 0x0303030303030303ull;
    kron.b16.words[0] = _pdep_u64(comb & 0xFFFF, extmask);
    kron.b16.words[1] = _pdep_u64(comb >> 16, extmask);
    return kron;
}

GR_repr GR4_16::kronecker_mul(GR_repr a, GR_repr b) const
{

    /* we use different representation of polynomials than before here.
     * each bit string can be split to sets of 2 bits where each set
     * corresponds to a coefficient modulo 4. */
    kronecker_form aa = this->kronecker_substitution(a);
    kronecker_form bb = this->kronecker_substitution(b);

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
    GR_repr ret;
    ret.lo = _pext_u64(tmp, loextmask);
    ret.hi = _pext_u64(tmp, hiextmask);
    return ret;
}


GR_repr GR4_n::euclid_rem(GR_repr a) const
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

/* https://dl.acm.org/doi/10.1016/j.ipl.2010.04.011 */
virtual GR_repr GR4_n::intel_rem(GR_repr a) const
{
    GR_repr hi = { a.hi >> this->n, a.lo >> this->n };
    GR_repr lo = { a.hi & this->mask, a.lo & this->mask };

    /* deg n-2 * deg n*/
    GR_repr rem = { 0x0, 0x0 };
    for (uint i = 0; i < this->q_plus.size() / 3; i++)
    {
        GR_repr tmp = this->mul_const(
            hi,
            { this->q_plus[3*i + 1], this->q_plus[3*i + 2] }
            );
        tmp <<= this->q_plus[3*i + 0];
        rem = this->add(rem, tmp);
    }
    rem >>= this->n;
    /* deg n-1 * deg n - 2*/
    GR_repr r = { 0x0, 0x0 };
    for (uint i = 0; i < this->mod_ast.size(); i++)
    {
        GR_repr tmp = rem << this->mod_ast[i];
        r = this->add(r, tmp);
    }

    r &= this->mask;
    return this->subtract(lo, r);
}

GR_repr GR4_16::intel_rem(GR_repr a) const
{
    GR_repr hi = a >> 16;
    GR_repr lo = a & 0xFFFF;

    GR_repr tmp = hi >> 14;
    tmp = this->add(tmp, hi >> 13);
    tmp = this->add(tmp, hi >> 11);
    tmp = this->subtract(hi, tmp);

    GR_repr r = this->add(tmp, tmp << 2);
    r = this->add(r, tmp << 3);
    r = this->add(r, tmp << 5);

    r &= 0xFFFF;
    return this->subtract(lo, r);
}

GR_repr GR4_32::intel_rem(GR_repr a) const
{
    GR_repr hi = a >> 32;
    GR_repr lo = a & 0xFFFFFFFF;

    GR_repr tmp = hi >> 30;
    tmp = this->add(tmp, hi >> 29);
    tmp = this->add(tmp, hi >> 25);
    tmp = this->subtract(hi, tmp);

    GR_repr r = this->add(tmp, tmp << 2);
    r = this->add(r, tmp << 3);
    r = this->add(r, tmp << 7);

    r &= 0xFFFFFFFF;
    return this->subtract(lo, r);
}

GR_repr GR4_n::mont_rem(GR_repr a) const
{
    /* n-1 deg + n-1 deg */
    GR_repr u = this->mul(a, this->n_prime);

    u &= this->mask;
    /* n deg + n-1 deg */
    GR_repr c = this->add(a, this->fast_mul(u, {0, this->mod}));

    c >>= this->n;
    return c;
}
