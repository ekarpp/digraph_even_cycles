/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef SOLVER_TEST_H
#define SOLVER_TEST_H

#include "test.hh"

class Solver_test : Test
{
private:
    const int n;

    bool test_solver();

public:
    bool run()
    {
        this->start_tests("solver");
        return test_solver();
    }

    Solver_test(const int n, const int tests = 0): n(n)
    {
        if (tests)
            this->tests = tests;
    }
};

#endif
