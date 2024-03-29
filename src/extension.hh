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
class GR_element;

/* representation for elements of E(4^n)
 * each bit in lo is the low bit of the mod 4 coefficient.
 * similarly for hi
 */
class GR_repr
{
public:
    uint64_t hi;
    uint64_t lo;

    inline GR_repr &operator>>=(const int n)
    {
        this->hi >>= n;
        this->lo >>= n;

        return *this;
    }

    inline GR_repr &operator<<=(const int n)
    {
        this->hi <<= n;
        this->lo <<= n;

        return *this;
    }

    inline GR_repr &operator&=(const uint64_t m)
    {
        this->hi &= m;
        this->lo &= m;

        return *this;
    }

    inline GR_repr operator>>(const int n) const
    {
        return { this->hi >> n, this->lo >> n };
    }

    inline GR_repr operator<<(const int n) const
    {
        return { this->hi << n, this->lo << n };
    }

    inline GR_repr operator&(const uint64_t m) const
    {
        return { this->hi & m, this->lo & m };
    }

    inline GR_repr shiftr_and(const int n, const uint64_t m) const
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
class GR4_n
{
private:
    const int n;
    const uint64_t mod;
    uint64_t mask;

    std::vector<uint64_t> mod_ast;
    std::vector<uint64_t> q_plus;

    GR_repr n_prime;
    GR_repr r_squared;

    /* euclidean division, only used once during initialization.
     * b has to be monic for this to work */
    GR_repr quo(GR_repr a, GR_repr b) const;

    void init_varying_size();

public:
    GR4_n(const int e, const uint64_t g);

    inline GR_repr add(const GR_repr &a, const GR_repr &b) const
    {
        const uint64_t carry = a.lo & b.lo;
        return { carry ^ a.hi ^ b.hi, a.lo ^ b.lo };
    }

    inline GR_repr negate(const GR_repr &a) const
    {
        return {
            a.lo ^ a.hi,
            a.lo
        };
    }

    inline GR_repr subtract(const GR_repr &a, const GR_repr &b) const
    {
        return this->add(a, this->negate(b));
    }

    inline GR_repr mul_const(const GR_repr &a, const GR_repr &c) const
    {
        return {
            (a.hi & c.lo) ^ (a.lo & c.hi),
            a.lo & c.lo
        };
    }

    inline GR_repr mul(const GR_repr &a, const GR_repr &b) const
    {
        return kronecker_mul(a, b);
    }

    GR_repr ref_mul(const GR_repr &a, const GR_repr &b) const;

    GR_repr fast_mul(const GR_repr &a, const GR_repr &b) const;

    virtual GR_repr kronecker_mul(const GR_repr &a, const GR_repr &b) const;
    virtual kronecker_form kronecker_substitution(const GR_repr &x) const;

    inline GR_repr rem(const GR_repr &a) const
    {
        return intel_rem(a);
    }

    GR_repr euclid_rem(const GR_repr &a) const;

    virtual GR_repr intel_rem(const GR_repr &a) const;

    GR_repr mont_rem(const GR_repr &a) const;
    inline GR_repr mont_form(const GR_repr &a) const
    {
        return this->mont_rem(this->mul(a, this->r_squared));
    }
    inline GR_repr mont_reduce(const GR_repr &a) const
    {
        return this->mont_rem(this->mul(a, { 0, 1 }));
    }

    inline int get_n() const { return this->n; }
    inline uint64_t get_mod() const { return this->mod; }
    inline uint64_t get_mask() const { return this->mask; }
};

class GR4_16 : public GR4_n
{
public:
    using GR4_n::GR4_n;

    kronecker_form kronecker_substitution(const GR_repr &x) const override;
    GR_repr kronecker_mul(const GR_repr &a, const GR_repr &b) const override;

    GR_repr intel_rem(const GR_repr &a) const override;
};

class GR4_32 : public GR4_n
{
public:
    using GR4_n::GR4_n;

    GR_repr intel_rem(const GR_repr &a) const override;
};

class GR_element
{
private:
    GR_repr repr;

public:
    GR_element(): repr{ 0, 0 } { };

    GR_element(const uint64_t hi, const uint64_t lo): repr{ hi, lo } { };

    explicit GR_element(const GR_repr repr): repr(repr) { };

    GR_element(const GR_element& e): repr{ e.get_hi(), e.get_lo() } { };

    inline GR_element operator+(const GR_element &other) const
    {
        return GR_element(
            global::E->add(this->repr, other.get_repr())
        );
    }

    inline GR_element &operator+=(const GR_element &other)
    {
        this->repr = global::E->add(this->repr, other.get_repr());
        return *this;
    }

    inline GR_element operator-(const GR_element &other) const
    {
        return GR_element(
            global::E->subtract(this->repr, other.get_repr())
        );
    }

    inline GR_element &operator-=(const GR_element &other)
    {
        this->repr = global::E->subtract(this->repr, other.get_repr());
        return *this;
    }

    inline GR_element operator*(const GR_element &other) const
    {
        GR_repr prod = global::E->mul(this->repr, other.get_repr());
        return GR_element(global::E->rem(prod));
    }

    inline GR_element &operator*=(const GR_element &other)
    {
        GR_repr prod = global::E->mul(this->repr, other.get_repr());
        this->repr = global::E->rem(prod);
        return *this;
    }

    inline bool operator==(const GR_element &other) const
    {
        return this->repr.lo == other.get_lo()
            && this->repr.hi == other.get_hi();
    }

    inline bool is_even() const
    {
        return this->repr.lo == 0x0;
    }

    /* used only on elements that are multiplied by two
     * thus we can just move the hi to low
     * maybe even just return gf element straight away
     * as this (probably?) gets anyways done after div2 */
    /* modify instead of returning new? */
    inline GR_element div2() const
    {
        return GR_element(0x0, this->repr.hi);
    }

    inline uint64_t get_lo() const { return this->repr.lo; }
    inline uint64_t get_hi() const { return this->repr.hi; }
    inline GR_repr get_repr() const { return this->repr; }

    GF_element project() const;

    inline GR_element &operator=(const GR_element &other)
    {
        this->repr.lo = other.get_lo();
        this->repr.hi = other.get_hi();
        return *this;
    }

    inline bool operator!=(const GR_element &other) const
    {
        return !(*this == other);
    }

    void print() const
    {
        for (int i = 15; i >= 0; i--)
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
    inline GR_element tau(const GR_element &sigma, const GR_element &v)
    {
        return sigma.project().inv().lift() * v;
    }

    inline GR_element GR_zero()
    {
        return GR_element(0b0, 0b0);
    }

    inline GR_element GR_one()
    {
        return GR_element(0b0, 0b1);
    }

    inline GR_element GR_random()
    {
        return GR_element(
            global::randgen() & global::E->get_mask(),
            global::randgen() & global::E->get_mask()
        );
    }
}

#endif
