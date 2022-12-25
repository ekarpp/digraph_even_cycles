/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <iostream>

#include "gf.hh"

class Polynomial
{
private:
    std::vector<GF_element> coeffs;
    const int deg;

public:
    /* be lazy and just store coefficients in vector of length n.
     * dont care if some of the coefficients are zero */
    explicit Polynomial(const int n): coeffs(n+1), deg(n) {};
    explicit Polynomial(const std::vector<GF_element> &P):
        coeffs(P), deg(P.size() - 1) {};

    void div(const GF_element &v);

    Polynomial &operator*=(const GF_element &other);

    Polynomial &operator+=(const Polynomial &other);

    GF_element operator[](int i) const
    {
        return this->coeffs[i];
    }

    /* set coefficient with deg i to val */
    /* const & for val?? */
    void operator()(const int i, GF_element val)
    {
        this->coeffs[i] = val;
    }

    /* eval at point x. only used for testing */
    GF_element eval(const GF_element &x)
    {
        GF_element val = this->coeffs[0];
        GF_element prod = x;
        for (int i = 1; i <= this->deg; i++)
        {
            val += prod * this->coeffs[i];
            prod *= x;
        }
        return val;
    }

    void print() const
    {
        for (int i = 0; i <= this->deg; i++)
        {
            std::cout << i << ": ";
            this->coeffs[i].print();
        }
    }
};

namespace util
{
    Polynomial poly_interpolation(
        const std::vector<GF_element> &gamma,
        const std::vector<GF_element> &delta
    );
}

#endif
