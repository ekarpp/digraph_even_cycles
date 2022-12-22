/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <vector>

#include "polynomial.hh"
#include "global.hh"

using namespace std;
/* be lazy and just store coefficients in vector of length n.
 * dont care if some of the coefficients are zero */
Polynomial::Polynomial(int n): coeffs(n+1, util::GF_zero())
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
    this->coeffs[this->deg] = util::GF_zero();

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

namespace util
{
    /* la grange interpolation with gamma and delta
     * note that we are in characteristic 2 and thus
     * - = +. done with the formula (3.3) here:
     * https://doi.org/10.1137/S0036144502417715 */
    Polynomial poly_interpolation(
        const std::vector<GF_element> &gamma,
        const std::vector<GF_element> &delta
        )
    {
        // assert(gamma.size() == delta.size())
        // assert(n > 2)
        int n = gamma.size();
        Polynomial interp(n - 1);

        /* weights*/
        std::vector<GF_element> w(n, util::GF_one());
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
        std::vector<GF_element> P(n+1, util::GF_zero());
        P[n] = util::GF_one();
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
