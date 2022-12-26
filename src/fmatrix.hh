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

    EMatrix lift() const;

    /* multiply diagonal by e. merge this with lift, so
     * that only one new copy is created? lift gets always
     * called after this */
    EMatrix mul_diag_lift(const GF_element &e) const;

    void mul_gamma(const int r1, const int r2, const GF_element &gamma);

    void swap_rows(const int r1, const int r2, const int idx = 0);

    /* uses gaussian elimination with pivoting.
     * modifies the object it is called on. */
    GF_element det();

    /* det of the matrix we get when r1 is multiplied by monomials
     * (1,r,..,r^(n-1)) and r2 by monomials (r^(n-1),..,r,1) */
    Polynomial pdet(int r1, int r2) const;

    /* return pcc_{n-1} of the matrix we get when we
     * multiply the diagonal of this matrix by e */
    GF_element pcc(const GF_element &e) const;
};

#endif
