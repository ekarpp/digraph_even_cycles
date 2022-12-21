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
    /* la grange interpolation with gamma and delta
     * note that we are in characteristic 2 and thus
     * - = +. done with the formula (3.3) here:
     * https://doi.org/10.1137/S0036144502417715 */
    inline Polynomial poly_interpolation(
        const std::vector<GF_element> &gamma,
        const std::vector<GF_element> &delta
    )
    {
        // assert(gamma.size() == delta.size())
        // assert(n > 2)
        int n = gamma.size();
        Polynomial interp(n - 1);

        /* weights*/
        std::vector<GF_element> w(n, global::F->one());
        for (int j = 1; j < n; j++)
        {
            for (int k = 0; k < j; k++)
            {
                w[k] *= gamma[k] + gamma[j];
                w[j] *= gamma[k] + gamma[j];
            }
        }

        for (int j = 0; j < n; j++)
            w[j].inv_in_place();

        /* main polynomial [ prod_{i} (x + gamma_i) ]*/
        std::vector<GF_element> P(n+1, global::F->zero());
        P[n] = global::F->one();
        P[n-1] = gamma[0];
        for (int i = 1; i < n; i++)
        {
            for (int j = n - i - 1; j < n - 1; j++)
                P[j] += gamma[i] * P[j+1];
            P[n - 1] += gamma[i];
        }

        for (int i = 0; i < n; i++)
        {
            Polynomial tmp(P);
            tmp.div(gamma[i]);
            tmp *= w[i] * delta[i];
            interp += tmp;
        }

        return interp;
    }
}

#endif
