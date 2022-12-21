/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>

#include "util_test.hh"
#include "../../src/gf.hh"
#include "../../src/polynomial.hh"
#include "../../src/util.hh"
#include "../../src/global.hh"

using namespace std;

Util_test::Util_test(int n, int tests = 0)
{
    if (tests)
        this->tests = tests;
    this->n = n;
    cout << "------------" << endl;
    cout << "TESTING UTIL" << endl;
    cout << "------------" << endl;

    this->run();
}

void Util_test::test_interpolation()
{
    cout << "polynomial interpolation: ";
    int err = 0;
    for (int t = 0; t < this->tests; t++)
    {
        uint64_t g = global::F->rem(global::randgen());
        uint64_t d = global::F->rem(global::randgen());

        vector<GF_element> gamma(n);
        vector<GF_element> delta(n);

        for (int i = 0; i < n; i++)
        {
            gamma[i] = GF_element(g);
            g = global::F->rem(g + 1);

            delta[i] = GF_element(d);
            d = global::F->rem(d + 1);
        }

        Polynomial p = util::poly_interpolation(gamma, delta);

        for (int i = 0; i < n; i++)
        {
            if (p.eval(gamma[i]) != delta[i])
            {
                err++;
                break;
            }
        }
    }
    end_test(err);
}

void Util_test::test_log2()
{
    cout << "log2: ";
    int err = 0;
    for (int t = 0; t < this->tests; t++)
    {
        int shift = global::randgen() % 64;
        uint64_t bit = 1ll << shift;

        if (shift != util::log2(bit))
            err++;
    }
    end_test(err);
}
