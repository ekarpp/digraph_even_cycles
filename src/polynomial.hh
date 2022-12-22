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
    int deg;

public:
    Polynomial(int n);
    Polynomial(std::vector<GF_element> P);

    void div(GF_element v);

    Polynomial &operator*=(GF_element val);

    Polynomial &operator+=(const Polynomial &other);

    GF_element operator[](int i) const
    {
        return this->coeffs[i];
    }

    /* set coefficient with deg i to val */
    void operator()(int i, GF_element val)
    {
        this->coeffs[i] = val;
    }

    /* eval at point x. only used for testing */
    GF_element eval(GF_element x)
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
