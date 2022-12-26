/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef MATRIX_H
#define MATRIX_H

#include <valarray>

template <typename T>
class Matrix
{
private:
    const int n;
    std::valarray<T> m;

public:
    /* for graph.cc */
    Matrix() { }
    explicit Matrix(const int d): n(d), m(d*d) { }
    Matrix(const int d, const std::valarray<T> &matrix): n(d), m(d*d)
    {
        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                this->m[i*this->n + j] = matrix[i*this->n + j];
    }

    inline int get_n() const { return this->n; }

    inline const T &operator()(const int row, const int col) const
    {
        return this->m[row*this->n + col];
    }

    inline bool operator==(const Matrix<T> &other) const
    {
        if (this->n != other.get_n())
            return false;

        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                if (this->operator()(i,j) != other(i,j))
                    return false;

        return true;
    }

    inline bool operator!=(const Matrix<T> &other) const
    {
        return !(*this == other);
    }

    inline void set(const int row, const int col, const T &val)
    {
        this->m[row*this->n + col] = val;
    }

    /* mul M_{r,c} with v*/
    inline void mul(const int r, const int c, const T &v)
    {
        this->m[r*this->n + c] *= v;
    }

    /* multiply row row with v, starting from column idx */
    inline void mul_row(const int row, const T &v, const int idx = 0)
    {
        for (int col = idx; col < this->get_n(); col++)
            this->m[row*this->n + col] *= v;
    }

    /* subtract v times row r1 from row r2, starting from column idx */
    inline void row_op(const int r1, const int r2, const T &v, const int idx = 0)
    {
        for (int col = idx; col < this->n; col++)
            this->m[r2*this->n + col] -= v * this->operator()(r1, col);
    }

    /* copy values from other to this */
    inline void copy(const Matrix<T> &other)
    {
        // assert this.n == other.n
        for (int row = 0; row < this->n; row++)
            for (int col = 0; col < this->n; col++)
                this->m[row*n + col] = other(row, col);
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
