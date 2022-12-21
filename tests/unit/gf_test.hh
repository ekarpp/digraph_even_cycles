/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef GF_TEST_H
#define GF_TEST_H

#include "test.hh"

class GF_test : Test
{
private:
    int n;

    void test_add_inverse();
    void test_associativity();
    void test_mul_id();
    void test_mul_inverse();
    void test_lift_project();
#if GF2_bits == 16
    void test_wide_mul();
#endif

    void run()
    {
        test_add_inverse();
        test_associativity();
        test_mul_id();
        test_mul_inverse();
        test_lift_project();
#if GF2_bits == 16
        test_wide_mul();
#endif
    }

public:
    GF_test();

};

#endif
