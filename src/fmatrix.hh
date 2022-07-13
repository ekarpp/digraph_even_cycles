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

class FMatrix
{
private:
    std::valarray<GF_element> m;
    int n;

public:
    /* for graph.cc */
    FMatrix() {};
    FMatrix(int n, std::valarray<GF_element> matrix);

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

    /* returns a copy of this */
    FMatrix copy() const;

    const std::valarray<GF_element> &get_m() const { return this->m; }

    int get_n() const { return this->n; }

    void mul(int row, int col, const GF_element &v)
    {
        this->m[row*this->n + col] *= v;
    }

    void mul_row(int row, const GF_element &v)
    {
        for (int col = 0; col < this->n; col++)
            this->m[row*this->n + col] *= v;
    }

    void mul_gamma(int r1, int r2, const GF_element &gamma)
    {
        GF_element prod = gamma;
        for (int col = 1; col < this->n; col++)
        {
            this->mul(r1, col, prod);
            this->mul(r2, this->n - 1 - col, prod);
            prod *= gamma;
        }
    }

    /* subtract v times r1 from r2 */
    void row_op(int r1, int r2, GF_element v)
    {
        for (int col = 0; col < this->n; col++)
            this->m[r2*this->n + col] -= v*this->operator()(r1,col);
    }

    const GF_element &operator()(int row, int col) const
    {
        return this->m[row*this->n + col];
    }

    bool operator==(const FMatrix &other) const
    {
        if (this->n != other.get_n())
            return false;

        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                if (this->operator()(i,j) != other(i,j))
                    return false;

        return true;
    }

    bool operator!=(const FMatrix &other) const
    {
        return !(*this == other);
    }

    void set(int row, int col, GF_element val)
    {
        this->m[row*this->n + col] = val;
    }

    /* swap rows r1 and r2 starting from column idx */
    void swap_rows(int r1, int r2, int idx = 0)
    {
        GF_element tmp;
        for (int col = idx; col < this->n; col++)
        {
            tmp = this->operator()(r1, col);
            this->set(r1, col, this->operator()(r2,col));
            this->set(r2, col, tmp);
        }
    }

    void print() const
    {
        for (int row = 0; row < this->n; row++)
        {
            for (int col = 0; col < this->n; col++)
                std::cout << std::bitset<8>(this->operator()(row, col).get_repr()) << " ";
            std::cout << std::endl;
        }
    }
};

#endif
