/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef UTIL_TEST_H
#define UTIL_TEST_H

#include "test.hh"

class Util_test : public Test
{
private:
    int n;

    bool test_interpolation();
    bool test_log2();

public:
    using Test::Test;

    bool run(const int deg = 0)
    {
        if (deg)
            this->n = deg;
        this->start_tests("util");

        return test_interpolation() | test_log2();
    }
};

#endif
