/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <stdint.h>
#include <iostream>
#include <unordered_set>

#include "gf.hh"
#include "extension.hh"
#include "global.hh"

using namespace std;

GF2_n::GF2_n(const int &e, const uint64_t &g): n(e), mod(g)
{
    this->mask = (1ll << this->n) - 1;
    if (this->n != 16 && this->n != 32)
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
uint64_t GF2_n::quo(uint64_t a, const uint64_t b) const
{
    uint64_t q = 0b0;
    const int degb = util::log2(b);
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
uint64_t GF2_n::ext_euclid(const uint64_t a) const
{
    // assert(a != 0)
    uint64_t s0 = 1;
    uint64_t s1 = 0;

    uint64_t r0 = a;
    uint64_t r1 = this->mod;

    /* invariants NOT:
     * x^{shift}*r0 = a*s0 + b*t0
     * x^{shift}*r1 = a*s1 + b*t1
     */

    int shift = __builtin_ctzl(r0);
    r0 >>= shift;

    while (r0 != r1)
    {
        int count = 0;
        if (r0 > r1)
        {
            r0 ^= r1;
            count = __builtin_ctzl(r0);
            r0 >>= count;

            s0 ^= s1;
            s1 <<= count;
        }
        else
        {
            r1 ^= r0;
            count = __builtin_ctzl(r1);
            r1 >>= count;

            s1 ^= s0;
            s0 <<= count;
        }
        shift += count;
    }

    for (int i = 0; i < shift; i++)
    {
        if ((s0 & 1) == 1)
            s0 ^= this->mod;
        s0 >>= 1;
    }

    return s0;
}

/* returns r s.t. for some q,
 * a = q*field.mod + r is the division relation (in Z(2^n))
 */
uint64_t GF2_n::rem(const uint64_t a) const
{
    const uint64_t lo = a & this->mask;
    const uint64_t hi = a >> this->n;

    uint64_t r = this->clmul(hi, this->q_plus);
    r >>= this->n;
    r = this->clmul(r, this->mod_ast);
    r &= this->mask;
    return r ^ lo;
}

uint64_t GF2_16::rem(const uint64_t a) const
{
    const uint64_t lo = a & 0xFFFF;
    const uint64_t hi = a >> 16;

    uint64_t r = hi ^ (hi >> 14) ^ (hi >> 13) ^ (hi >> 11);
    r ^= (r << 2) ^ (r << 3) ^ (r << 5);
    r &= 0xFFFF;
    return r ^ lo;
}

uint64_t GF2_32::rem(const uint64_t a) const
{
    const uint64_t lo = a & 0xFFFFFFFF;
    const uint64_t hi = a >> 32;

    uint64_t r = hi ^ (hi >> 30) ^ (hi >> 29) ^ (hi >> 25);
    r ^= (r << 2) ^ (r << 3) ^ (r << 7);
    r &= 0xFFFFFFFF;
    return r ^ lo;
}

GR_element GF_element::lift() const
{
    return GR_element(0x0, this->repr);
}

namespace util
{
    /* returns n distinct random elements from
     * global::F-> (use LSFR?) */
    std::vector<GF_element> distinct_elements(const int n)
    {
        std::vector<GF_element> vec(n);
        std::unordered_set<uint64_t> have;
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
