/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef UTIL_TEST_H
#define UTIL_TEST_H

#include "test.hh"

class Util_test : Test
{
private:
    const int n;

    bool test_interpolation();
    bool test_log2();

public:
    Util_test(const int n, const int tests = 0): n(n)
    {
        if (tests)
            this->tests = tests;
    }

    bool run()
    {
        this->start_tests("util");

        return test_interpolation() | test_log2();
    }
};

#endif
