/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef EMATRIX_H
#define EMATRIX_H

#include <iostream>
#include <valarray>

#include "extension.hh"
#include "matrix.hh"
#include "fmatrix.hh"

/* forward declare */
class FMatrix;

class EMatrix : public Matrix<GR_element>
{
public:
    using Matrix::Matrix;

    FMatrix project() const;

    /* returns Per(this) - Det(this) as described in chapter 3
     * of the paper*/
    GR_element per_m_det();

    GR_element row_op_per(const int i1, const int j);

    GR_element per_similar(const int i1, const int i2);
};

#endif
