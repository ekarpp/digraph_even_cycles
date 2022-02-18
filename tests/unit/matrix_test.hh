/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef MATRIX_TEST_H
#define MATRIX_TEST_H

#include "../../src/gf.hh"
#include "../../src/fmatrix.hh"
#include "test.hh"

class Matrix_test : Test
{
private:
    int dim;

    void test_addition();
    void test_subtraction();
    void test_multiplication();

    FMatrix random_matrix();

    void run()
    {
        test_addition();
        test_subtraction();
        test_multiplication();
    }

public:
    Matrix_test(int dim, int tests);
};

#endif