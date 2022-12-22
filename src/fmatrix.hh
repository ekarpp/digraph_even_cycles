/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef FMATRIX_H
#define FMATRIX_H

#include <iostream>
#include <bitset>
#include <vector>
#include <valarray>

#include "gf.hh"
#include "ematrix.hh"
#include "polynomial.hh"

/* forward declare */
class EMatrix;

class FMatrix : public Matrix<GF_element>
{
public:
    using Matrix::Matrix;

    /* return pcc_{n-1} of the matrix we get when we
     * multiply the diagonal of this matrix by e */
    GF_element pcc(const GF_element &e) const;

    EMatrix lift() const;

    /* multiply diagonal by e. merge this with lift, so
     * that only one new copy is created? lift gets always
     * called after this */
    FMatrix mul_diag(const GF_element &e) const;

    /* uses gaussian elimination with pivoting.
     * modifies the object it is called on. */
    GF_element det();

    /* det of the matrix we get when r1 is multiplied by monomials
     * (1,r,..,r^(n-1)) and r2 by monomials (r^(n-1),..,r,1) */
    Polynomial pdet(int r1, int r2) const;

    inline void mul_gamma(int r1, int r2, const GF_element &gamma)
    {
        GF_element prod = gamma;
        for (int col = 1; col < this->get_n(); col++)
        {
            this->mul(r1, col, prod);
            this->mul(r2, this->get_n() - 1 - col, prod);
            prod *= gamma;
        }
    }

    /* swap rows r1 and r2 starting from column idx */
    inline void swap_rows(int r1, int r2, int idx = 0)
    {
        GF_element tmp;
        for (int col = idx; col < this->get_n(); col++)
        {
            tmp = this->operator()(r1, col);
            this->set(r1, col, this->operator()(r2,col));
            this->set(r2, col, tmp);
        }
    }
};

#endif
