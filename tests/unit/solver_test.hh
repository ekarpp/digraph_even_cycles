/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef SOLVER_TEST_H
#define SOLVER_TEST_H

#include "test.hh"

class Solver_test : public Test
{
private:
    int n = 5;

    bool test_solver();

public:
    using Test::Test;

    bool run(const int deg = 0)
    {
        if (deg)
            this->n = deg;
        this->start_tests("solver");
        return test_solver();
    }
};

#endif
