/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <valarray>
#include <vector>

#include "global.hh"
#include "fmatrix.hh"
#include "ematrix.hh"
#include "extension.hh"
#include "packed_fmatrix.hh"

using namespace std;

FMatrix::FMatrix(const int d, const valarray<GF_element> &matrix): n(d), m(d*d)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            this->m[i*n + j] = matrix[i*n + j];
}

EMatrix FMatrix::lift() const
{
    valarray<GR_element> lifted(this->get_n() * this->get_n());

    for (int x = 0; x < this->get_n(); x++)
    {
        for (int y = 0; y < this->get_n(); y++)
            lifted[x*this->get_n() + y] = this->operator()(x,y).lift();
    }

    return EMatrix(this->get_n(), lifted);
}

FMatrix FMatrix::mul_diag(const GF_element &e) const
{
    valarray<GF_element> m(this->get_n() * this->get_n());

    for (int row = 0; row < this->get_n(); row++)
    {
        for (int col = 0; col < this->get_n(); col++)
        {
            if (row == col)
                m[row*this->get_n() + col] = this->operator()(row,col) * e;
            else
                m[row*this->get_n() + col] = this->operator()(row,col);
        }
    }

    return FMatrix(this->get_n(), m);
}

/* simple gaussian elimination with pivoting.
 * we are in characteristic two so pivoting does
 * not affect the determinant. */
GF_element FMatrix::det()
{
    GF_element det = util::GF_one();
    for (int col = 0; col < this->get_n(); col++)
    {
        /* pivot */
        GF_element mx = util::GF_zero();
        int mxi = -1;
        for (int row = col; row < this->get_n(); row++)
        {
            if (this->operator()(row,col) != util::GF_zero())
            {
                mx = this->operator()(row,col);
                mxi = row;
                break;
            }
        }

        if (mx == util::GF_zero())
            return util::GF_zero();

        if (mxi != col)
            this->swap_rows(mxi, col);

        det *= mx;
        mx.inv_in_place();
        this->mul_row(col, mx);
        for (int row = col+1; row < this->get_n(); row++)
            this->row_op(col, row, this->operator()(row,col));
    }
    return det;
}

/* uses random sampling and la grange interpolation
 * to get the polynomial determinant. */
Polynomial FMatrix::pdet(int r1, int r2) const
{
    /* determinant has deg <= 2*n - 2 */
    vector<GF_element> gamma = util::distinct_elements(2*this->get_n() - 1);
    vector<GF_element> delta(2*this->get_n() - 1);

<<<<<<< HEAD
    if (global::F->get_n() != 16)
=======
    for (int i = 0; i < 2*this->get_n() - 1; i++)
>>>>>>> 0f7b10e (add matrix.hh. matrix tests fail)
    {
        FMatrix A(this->n);

        for (int i = 0; i < 2*this->n - 1; i++)
        {
<<<<<<< HEAD
            A.copy(*this);
=======
            FMatrix A = this->copy<FMatrix>();
>>>>>>> 0f7b10e (add matrix.hh. matrix tests fail)
            A.mul_gamma(r1, r2, gamma[i]);
            delta[i] = A.det();
        }
        /* la grange */
        return util::poly_interpolation(gamma, delta);
    }

    Packed_FMatrix PA(this->n);

    for (int i = 0; i < 2*this->n - 1; i++)
    {
        PA.init(*this);
        PA.mul_gamma(r1, r2, gamma[i]);
        delta[i] = PA.det();
    }

    return util::poly_interpolation(gamma, delta);

}

<<<<<<< HEAD
/* copy valarray instead?? */
void FMatrix::copy(const FMatrix &A)
{
    for (int row = 0; row < this->n; row++)
        for (int col = 0; col < this->n; col++)
            this->m[row*n + col] = A(row,col);
}

=======
>>>>>>> 0f7b10e (add matrix.hh. matrix tests fail)
GF_element FMatrix::pcc(const GF_element &e) const
{
    EMatrix E = this->mul_diag(e).lift();
    GR_element elem = E.per_m_det();
    return elem.div2().project();
}
