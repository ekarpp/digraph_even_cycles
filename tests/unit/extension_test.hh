/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef EXTENSION_TEST_H
#define EXTENSION_TEST_H

#include "test.hh"

class GR_test : Test
{
private:
    int n;

    bool test_add_inverse();
    bool test_associativity();
    bool test_mul();
    bool test_fast_mul();
    bool test_intel_rem();
    bool test_mont_rem();
    bool test_even_tau();
    bool test_is_even();
    bool test_kronecker_mul();

public:
    GR_test(int tests)
    {
        this->tests = tests;
    }

    bool run()
    {
        this->start_tests("extension");

        return test_add_inverse() | test_associativity()
            | test_mul() | test_even_tau() | test_is_even()
            | test_fast_mul() | test_intel_rem() | test_mont_rem()
            | test_kronecker_mul();
    }
};

#endif
