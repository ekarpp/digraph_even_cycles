#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

/* parent class for square matrices with different elements
   contains the usual arithmetic operations */

template <typename E>
class Matrix
{
private:
    std::vector<std::vector<E>> m;
    int n;
    E one;
    E zero;

public:
    Matrix() {}

    Matrix(int n, E one, std::vector<std::vector<E>> &matrix)
    {
        this->n = n;
        this->one = one;
        this->zero = one - one;
        std::vector<std::vector<E>> v(n, std::vector<E>(n));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                v[i][j] = matrix[i][j];
        this->m = v;
    }

    Matrix &operator+=(const Matrix &other)
    {
        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                this->set(i, j, this->operator()(i,j) + other(i,j));

        return *this;
    }

    Matrix operator+(const Matrix &other) const
    {
        std::vector<std::vector<E>> sum(n, std::vector<E>(n));

        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                sum[i][j] = this->operator()(i,j) + other(i,j);

        return Matrix(this->n, this->one, sum);
    }

    Matrix operator-(const Matrix &other) const
    {
        std::vector<std::vector<E>> sum(n, std::vector<E>(n));

        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                sum[i][j] = this->operator()(i,j) - other(i,j);

        return Matrix(this->n, this->one, sum);
    }

    Matrix operator*(const Matrix &other) const
    {
        std::vector<std::vector<E>> prod(n, std::vector<E>(n));

        for (int i = 0; i < this->n; i++)
        {
            for (int j = 0; j < this->n; j++)
            {
                prod[i][j] = this->zero;
                for (int k = 0; k < this->n; k++)
                    prod[i][j] += this->operator()(i,k) * other(k,j);
            }
        }

        return Matrix(this->n, this->one, prod);
    }

    std::vector<E> operator[](int i) const
    {
        /* do error checking in child classes */
        return this->m[i];
    }

    E operator()(int row, int col) const
    {
        /* append identity matrix in case out of bounds.
         * potentially dangerous */
        if (row >= this->n || col >= this->n)
            return (row == col) ? this->one : this->zero;
        return this->m[row][col];
    }

    int get_n() const { return this->n; }

    bool operator==(const Matrix &other) const
    {
        if (this->n != other.get_n())
            return false;

        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                if (this->operator()(i,j) != other(i,j))
                    return false;

        return true;
    }

    bool operator!=(const Matrix &other) const
    {
        return !(*this == other);
    }

    /* set m[row][col] = val. use this to avoid const issues with [] op */
    void set(int row, int col, E val)
    {
        this->m[row][col] = val;
    }
};

#endif
