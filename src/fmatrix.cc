/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <valarray>
#include <vector>

#include "global.hh"
#include "fmatrix.hh"
#include "ematrix.hh"
#include "extension.hh"
#include "packed_fmatrix.hh"

using namespace std;

/* multiply diagonal by e and lift the resulting matrix */
EMatrix FMatrix::mul_diag_lift(const GF_element &e) const
{
    EMatrix m(this->get_n());

    for (int row = 0; row < this->get_n(); row++)
    {
        for (int col = 0; col < this->get_n(); col++)
        {
            if (row == col)
                m.set(row, col, (this->operator()(row,col) * e).lift());
            else
                m.set(row, col, this->operator()(row,col).lift());
        }
    }

    return m;
}

void FMatrix::mul_gamma(const int r1, const int r2, const GF_element &gamma)
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
void FMatrix::swap_rows(const int r1, const int r2, const int idx)
{
    for (int col = idx; col < this->get_n(); col++)
    {
        const GF_element tmp = this->operator()(r1, col);
        this->set(r1, col, this->operator()(r2,col));
        this->set(r2, col, tmp);
    }
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
        int pivot_idx = -1;
        for (int row = col; row < this->get_n(); row++)
        {
            if (this->operator()(row,col) != util::GF_zero())
            {
                pivot_idx = row;
                break;
            }
        }

        if (pivot_idx == -1)
            return util::GF_zero();

        if (pivot_idx != col)
            this->swap_rows(pivot_idx, col, col);

        GF_element pivot = this->operator()(col, col);
        det *= pivot;
        pivot.inv_in_place();
        this->mul_row(col, pivot, col);

        for (int row = col+1; row < this->get_n(); row++)
            /* create new element or do row,col last? */
            this->row_op(col, row, GF_element(this->operator()(row,col)), col);
    }
    return det;
}

/* uses random sampling and la grange interpolation
 * to get the polynomial determinant. rows r1 and r2 are similar. */
Polynomial FMatrix::pdet(const int r1, const int r2) const
{
    /* determinant has deg <= 2*n - 2 */
    const vector<GF_element> gamma = util::distinct_elements(2*this->get_n() - 1);
    vector<GF_element> delta(2*this->get_n() - 1);

    if (global::F->get_n() != 16)
    {
        FMatrix A(this->get_n());

        for (int i = 0; i < 2*this->get_n() - 1; i++)
        {
            A.copy(*this);
            A.mul_gamma(r1, r2, gamma[i]);
            delta[i] = A.det();
        }
    }
    else
    {
        Packed_FMatrix PA(this->get_n(), *this);

        for (int i = 0; i < 2*this->get_n() - 1; i++)
        {
            PA.init();
            PA.mul_gamma(r1, r2, gamma[i]);
            delta[i] = PA.det();
        }
    }

    /* la grange */
    return util::poly_interpolation(gamma, delta);
}

GF_element FMatrix::pcc(const GF_element &e) const
{
    EMatrix E = this->mul_diag_lift(e);
    const GR_element elem = E.per_m_det();
    return elem.div2().project();
}
