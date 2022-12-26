/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef EMATRIX_TEST_H
#define EMATRIX_TEST_H

#include <valarray>

#include "test.hh"
#include "../../src/ematrix.hh"

class EMatrix_test : Test
{
private:
    int dim = 5;

    bool test_per_det();
    bool test_per_det_singular();

    EMatrix random();
    GR_element term(std::valarray<int> &perm, const EMatrix &m);
    void swap(int i1, int i2, std::valarray<int> &perm);
    GR_element per_m_det_heap(const EMatrix &m);

public:
    EMatrix_test(int tests = 0)
    {
        if (tests)
            this->tests = tests;
    }

    bool run()
    {
        this->start_tests("ematrix");

        return test_per_det() | test_per_det_singular();
    }
};

#endif
