/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef GENG_TEST_H
#define GENG_TEST_H

#include <vector>

#include "test.hh"

class Geng_test : Test
{
private:
    int n;

    bool test_geng();

    void store_graph(const std::vector<std::vector<int>> &g, int id);

public:
    Geng_test() { };

    bool run()
    {
        this->start_tests("geng");
        return test_geng();
    }
};

#endif
