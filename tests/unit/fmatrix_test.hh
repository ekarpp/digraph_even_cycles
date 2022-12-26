/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef FMATRIX_TEST_H
#define FMATRIX_TEST_H

#include <valarray>

#include "test.hh"
#include "../../src/gf.hh"
#include "../../src/fmatrix.hh"

class FMatrix_test : public Test
{
private:
    int dim = 10;

    GF_element det_heap(const FMatrix &m);
    GF_element term(std::valarray<int> &perm, const FMatrix &m);
    void swap(int i1, int i2, std::valarray<int> &perm);

    bool test_determinant_vandermonde();
    bool test_det_singular();
    bool test_determinant_random();
    bool test_pdet();
    bool test_packed_determinant();
    bool test_packed_determinant_singular();
    bool test_packed_gamma_mul();
    bool test_packed_init();

    FMatrix vandermonde();
    FMatrix random(int n);

public:
    using Test::Test;

    bool run(const int d = 0)
    {
        this->start_tests("fmatrix");
        if (d)
            this->dim = d;

        bool failure = test_pdet() | test_determinant_vandermonde()
            | test_determinant_random() | test_det_singular();

        if (global::F->get_n() == 16)
        {
            failure |= test_packed_init() | test_packed_determinant()
                | test_packed_determinant_singular() | test_packed_gamma_mul();
        }

        return failure;
    }
};

#endif
