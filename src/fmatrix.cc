/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <valarray>
#include <vector>

#include "global.hh"
#include "fmatrix.hh"
#include "ematrix.hh"
#include "extension.hh"

#if GF2_bits == 16
#include "packed_fmatrix.hh"
#endif

using namespace std;

FMatrix::FMatrix(int n, valarray<GF_element> matrix): m(n*n)
{
    this->n = n;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            m[i*n + j] = matrix[i*n + j];
}

EMatrix FMatrix::lift() const
{
    valarray<Extension_element> lifted(this->n * this->n);

    for (int x = 0; x < this->n; x++)
    {
        for (int y = 0; y < this->n; y++)
            lifted[x*this->n + y] = this->operator()(x,y).lift();
    }

    return EMatrix(this->n, lifted);
}

FMatrix FMatrix::mul_diag(const GF_element &e) const
{
    valarray<GF_element> m(this->n * this->n);

    for (int row = 0; row < this->n; row++)
    {
        for (int col = 0; col < this->n; col++)
        {
            if (row == col)
                m[row*this->n + col] = this->operator()(row,col) * e;
            else
                m[row*this->n + col] = this->operator()(row,col);
        }
    }

    return FMatrix(this->n, m);
}

/* simple gaussian elimination with pivoting.
 * we are in characteristic two so pivoting does
 * not affect the determinant. */
GF_element FMatrix::det()
{
    GF_element det = global::F.one();
    for (int col = 0; col < this->n; col++)
    {
        /* pivot */
        GF_element mx = global::F.zero();
        int mxi = -1;
        for (int row = col; row < this->n; row++)
        {
            if (this->operator()(row,col) != global::F.zero())
            {
                mx = this->operator()(row,col);
                mxi = row;
                break;
            }
        }

        if (mx == global::F.zero())
            return global::F.zero();

        if (mxi != col)
            this->swap_rows(mxi, col);

        det *= mx;
        mx.inv_in_place();
        this->mul_row(col, mx);
        for (int row = col+1; row < this->n; row++)
            this->row_op(col, row, this->operator()(row,col));
    }
    return det;
}

/* uses random sampling and la grange interpolation
 * to get the polynomial determinant. */
Polynomial FMatrix::pdet(int r1, int r2) const
{
    /* determinant has deg <= 2*n - 2 */
    vector<GF_element> gamma = util::distinct_elements(2*this->n - 1);
    vector<GF_element> delta(2*this->n - 1);

    for (int i = 0; i < 2*this->n - 1; i++)
    {
#if GF2_bits == 16
        Packed_FMatrix PA(*this);
        PA.mul_gamma(r1, r2, gamma[i]);
        delta[i] = PA.det();
#else
        FMatrix A = this->copy();
        A.mul_gamma(r1, r2, gamma[i]);
        delta[i] = A.det();
#endif
    }

    /* la grange */
    return util::poly_interpolation(gamma, delta);
}

FMatrix FMatrix::copy() const
{
    valarray<GF_element> m(this->n * this->n);

    for (int row = 0; row < this->n; row++)
        for (int col = 0; col < this->n; col++)
            m[row*n + col] = this->operator()(row,col);

    return FMatrix(this->n, m);
}

GF_element FMatrix::pcc(const GF_element &e) const
{
    EMatrix E = this->mul_diag(e).lift();
    Extension_element elem = E.per_m_det();
    return elem.div2().project();
}
