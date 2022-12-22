/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>

#include "extension_test.hh"
#include "../../src/extension.hh"
#include "../../src/global.hh"

using namespace std;


GR_test::GR_test(int tests)
{
    cout << "-----------------" << endl;
    cout << "TESTING EXTENSION" << endl;
    cout << "-----------------" << endl;
    this->tests = tests;
    this->run();
}

void GR_test::test_add_inverse()
{
    cout << "add inverse: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        GR_element e = util::GR_random();
        if (e - e != util::GR_zero())
            err++;
    }
    end_test(err);
}

void GR_test::test_associativity()
{
    cout << "test associativity: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        GR_element a = util::GR_random();
        GR_element b = util::GR_random();
        GR_element c = util::GR_random();
        if (a*(b+c) != c*a + b*a)
            err++;
    }
    end_test(err);
}

void GR_test::test_mul()
{
    cout << "test mul: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        GR_element a = util::GR_random();
        GR_element b = util::GR_random();

        if (a*b != b*a || a*util::GR_one() != a
            || b*util::GR_one() != b
            || a*util::GR_zero() != util::GR_zero()
            || b*util::GR_zero() != util::GR_zero())
            err++;
    }
    end_test(err);
}

void GR_test::test_fast_mul()
{
    cout << "test fast mul: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        GR_element a = util::GR_random();
        GR_element b = util::GR_random();

        GR_repr ref = global::E->ref_mul(a.get_repr(), b.get_repr());
        GR_repr fast = global::E->fast_mul(a.get_repr(), b.get_repr());

        if (fast.hi != ref.hi || fast.lo != ref.lo)
            err++;
    }
    end_test(err);
}

void GR_test::test_kronecker_mul()
{
    cout << "test kronecker mul: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        GR_element a = util::GR_random();
        GR_element b = util::GR_random();

        GR_repr ref = global::E->ref_mul(a.get_repr(), b.get_repr());
        GR_repr kron = global::E->kronecker_mul(a.get_repr(), b.get_repr());

        if (kron.hi != ref.hi || kron.lo != ref.lo)
            err++;
    }
    end_test(err);
}

void GR_test::test_intel_rem()
{
    cout << "test intel rem: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        GR_element a = util::GR_random();
        GR_element b = util::GR_random();
        GR_repr v = global::E->fast_mul(a.get_repr(), b.get_repr());

        GR_repr euclid = global::E->euclid_rem(v);
        GR_repr intel = global::E->intel_rem(v);

        if (euclid.hi != intel.hi || euclid.lo != intel.lo)
            err++;
    }
    end_test(err);
}

void GR_test::test_mont_rem()
{
    cout << "test montgomery multiplication: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        GR_repr a = {
            global::randgen() & global::E->get_mask(),
            global::randgen() & global::E->get_mask()
        };

        GR_repr b = {
            global::randgen() & global::E->get_mask(),
            global::randgen() & global::E->get_mask()
        };

        GR_repr mont = global::E->mont_reduce(
            global::E->mont_rem(
                global::E->mul(
                    global::E->mont_form(a),
                    global::E->mont_form(b)
                )
            )
        );

        GR_repr ref = global::E->euclid_rem(
            global::E->mul(a, b)
        );

        if (ref.hi != mont.hi || ref.lo != mont.lo)
            err++;
    }
    end_test(err);
}

void GR_test::test_even_tau()
{
    cout << "test even tau: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        GR_element sigma = util::GR_random();
        GR_element v = util::GR_random();
        if (sigma.is_even() || v.is_even())
            /* we get here with probability (0.5)^(d-1) */
            continue;
        GR_element e = v - sigma * util::tau(sigma, v);
        if (!e.is_even())
            err++;
    }
    this->end_test(err);
}

void GR_test::test_is_even()
{
    cout << "test is even: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        GR_element e(0x0, global::randgen() & global::E->get_mask());
        if (!e.is_even())
            err++;
    }
    this->end_test(err);
}
