/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef GF_TEST_H
#define GF_TEST_H

#include "test.hh"
#include "../../src/gf.hh"

class GF_test : Test
{
private:
    int n;

    bool test_add_inverse();
    bool test_associativity();
    bool test_mul_id();
    bool test_mul_inverse();
    bool test_lift_project();
    bool test_wide_mul();

public:
    GF_test() { };

    bool run()
    {
        this->start_tests("gf");

        bool failure = test_add_inverse() | test_associativity()
            | test_mul_id() | test_mul_inverse()
            | test_lift_project();

        if (global::F->get_n() == 16)
            failure |= test_wide_mul();

        return failure;
    }
};

#endif
