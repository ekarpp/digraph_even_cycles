/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef EMATRIX_H
#define EMATRIX_H

#include <iostream>
#include <valarray>

#include "extension.hh"
#include "fmatrix.hh"

/* forward declare */
class FMatrix;

class EMatrix
{
private:
    std::valarray<GR_element> m;
    int n;

public:
    EMatrix(int n, std::valarray<GR_element> matrix);

    FMatrix project() const;

    EMatrix copy() const;

    /* returns Per(this) - Det(this) as described in chapter 3
     * of the paper*/
    GR_element per_m_det();

    GR_element row_op(int i1, int j);

    int get_n() const { return this->n; }

    const GR_element &operator()(int row, int col) const
    {
        return this->m[row*this->n + col];
    }

    bool operator==(const EMatrix &other) const
    {
        if (this->n != other.get_n())
            return false;

        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                if (this->operator()(i,j) != other(i,j))
                    return false;

        return true;
    }

    bool operator!=(const EMatrix &other) const
    {
        return !(*this == other);
    }

    void set(int row, int col, GR_element val)
    {
        this->m[row*this->n + col] = val;
    }

    void print() const
    {
        for (int row = 0; row < this->n; row++)
        {
            for (int col = 0; col < this->n; col++)
                this->operator()(row, col).print();
            std::cout << std::endl;
        }
    }
};

#endif
