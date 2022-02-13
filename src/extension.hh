#ifndef EXTENSION_H
#define EXTENSION_H

#include <stdint.h>

/* representation for elements of E(4^n)
 * each bit in lo is the low bit of the mod 4 coefficient.
 * similarly for hi
 */
/* vectorize with AVX? */
struct int64_2_t
{
    int64_t hi;
    int64_t lo;
};

class Extension_element;

/* Extension of GF(2^n) to the ring E(4^n).
 * If GF(2^n) = Z2 / <g_2> for irreducible polynomial
 * g_2 of degree n, then if g_4 is g_2 but coefficients
 * projected to Z4 we have E(4^n) = Z4 / <g_4>. */
class Extension
{
private:
    int n;
    int64_t mod;
    int64_t mask;

public:
    Extension() {};
    void init(const int n, const int64_t mod);
    Extension_element zero();
    Extension_element one() const;
    Extension_element random();

    inline int64_2_t add(int64_2_t a, int64_2_t b) const;

    int get_n() { return this->n; }
    int64_t get_mod() { return this->mod; }
    int64_t get_mask() { return this->mask; }

};



class Extension_element
{
private:
    int64_2_t repr;

public:
    Extension_element(const int64_t lo, const int64_t hi);
    Extension_element(const int64_2_t repr);
    Extension_element operator+(const Extension_element &other);
    Extension_element operator-(const Extension_element &other);
    Extension_element operator*(const Extension_element &other);
    bool operator==(const Extension_element &other);

    int64_t get_lo() const { return this->repr.lo; }
    int64_t get_hi() const { return this->repr.hi; }
    int64_2_t get_repr() const { return this->repr; }

    Extension_element &operator=(const Extension_element &other)
    {
        this->repr.lo = other.get_lo();
        this->repr.hi = other.get_hi();
        return *this;
    }

    bool operator!=(const Extension_element &other)
    {
        return !(*this == other);
    }
};

#endif