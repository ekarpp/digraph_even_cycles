/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <vector>

#include "polynomial.hh"
#include "global.hh"

using namespace std;
/* be lazy and just store coefficients in vector of length n.
 * dont care if some of the coefficients are zero */
Polynomial::Polynomial(int n): coeffs(n+1, global::F->zero())
{
    this->deg = n;
}

Polynomial::Polynomial(vector<GF_element> P): coeffs(P)
{
    this->deg = P.size() - 1;
}

/* divides this by monomial (x + v) using synthetic division */
void Polynomial::div(GF_element v)
{
    GF_element prev = this->coeffs[this->deg];
    this->coeffs[this->deg] = global::F->zero();

    for (int i = this->deg - 1; i >= 0; i--)
    {
        GF_element tmp = this->coeffs[i];
        this->coeffs[i] = prev;
        prev *= v;
        prev += tmp;
    }
}

Polynomial &Polynomial::operator*=(GF_element val)
{
    for (int i = 0; i <= this->deg; i++)
        this->coeffs[i] *= val;

    return *this;
}

Polynomial &Polynomial::operator+=(const Polynomial &other)
{
    // assert(this->deg == other.deg)
    for (int i = 0; i <= this->deg; i++)
        this->coeffs[i] += other[i];

    return *this;
}
