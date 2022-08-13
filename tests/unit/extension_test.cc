/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>

#include "extension_test.hh"
#include "../../src/extension.hh"
#include "../../src/global.hh"

using namespace std;


Extension_test::Extension_test(int tests)
{
    cout << "-----------------" << endl;
    cout << "TESTING EXTENSION" << endl;
    cout << "-----------------" << endl;
    this->tests = tests;
    this->run();
}

void Extension_test::test_add_inverse()
{
    cout << "add inverse: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element e = global::E.random();
        if (e - e != global::E.zero())
            err++;
    }
    end_test(err);
}

void Extension_test::test_associativity()
{
    cout << "test associativity: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element a = global::E.random();
        Extension_element b = global::E.random();
        Extension_element c = global::E.random();
        if (a*(b+c) != c*a + b*a)
            err++;
    }
    end_test(err);
}

void Extension_test::test_mul()
{
    cout << "test mul: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element a = global::E.random();
        Extension_element b = global::E.random();

        if (a*b != b*a || a*global::E.one() != a
            || b*global::E.one() != b
            || a*global::E.zero() != global::E.zero()
            || b*global::E.zero() != global::E.zero())
            err++;
    }
    end_test(err);
}

void Extension_test::test_fast_mul()
{
    cout << "test fast mul: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element a = global::E.random();
        Extension_element b = global::E.random();

        extension_repr ref = global::E.ref_mul(a.get_repr(), b.get_repr());
        extension_repr fast = global::E.fast_mul(a.get_repr(), b.get_repr());

        if (fast.get_hi() != ref.get_hi() || fast.get_lo() != ref.get_lo())
            err++;
    }
    end_test(err);
}

void Extension_test::test_kronecker_mul()
{
    cout << "test kronecker mul: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element a = global::E.random();
        Extension_element b = global::E.random();

        extension_repr ref = global::E.ref_mul(a.get_repr(), b.get_repr());
        extension_repr kron = global::E.kronecker_mul(a.get_repr(), b.get_repr());

        if (kron.get_hi() != ref.get_hi() || kron.get_lo() != ref.get_lo())
            err++;
    }
    end_test(err);
}

void Extension_test::test_intel_rem()
{
    cout << "test intel rem: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element a = global::E.random();
        Extension_element b = global::E.random();
        extension_repr v = global::E.fast_mul(a.get_repr(), b.get_repr());

        extension_repr euclid = global::E.euclid_rem(v);
        extension_repr intel = global::E.intel_rem(v);

        if (euclid.get_hi() != intel.get_hi() || euclid.get_lo() != intel.get_lo())
            err++;
    }
    end_test(err);
}

void Extension_test::test_mont_rem()
{
    cout << "test montgomery multiplication: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        extension_repr a = {
            global::randgen() & global::E.get_mask(),
            global::randgen() & global::E.get_mask()
        };

        extension_repr b = {
            global::randgen() & global::E.get_mask(),
            global::randgen() & global::E.get_mask()
        };

        extension_repr mont = global::E.mont_reduce(
            global::E.mont_rem(
                global::E.mul(
                    global::E.mont_form(a),
                    global::E.mont_form(b)
                )
            )
        );

        extension_repr ref = global::E.euclid_rem(
            global::E.mul(a, b)
        );

        if (ref.get_hi() != mont.get_hi() || ref.get_lo() != mont.get_lo())
            err++;
    }
    end_test(err);
}

void Extension_test::test_even_tau()
{
    cout << "test even tau: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element sigma = global::E.random();
        Extension_element v = global::E.random();
        if (sigma.is_even() || v.is_even())
            /* we get here with probability (0.5)^(d-1) */
            continue;
        Extension_element e = v - sigma * util::tau(sigma, v);
        if (!e.is_even())
            err++;
    }
    this->end_test(err);
}

void Extension_test::test_is_even()
{
    cout << "test is even: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element e(0x0, global::randgen() & global::E.get_mask());
        if (!e.is_even())
            err++;
    }
    this->end_test(err);
}

#if GF2_bits == 16
void Extension_test::test_packed_intel_rem()
{
    cout << "test packed intel rem: ";
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        Extension_element a;
        Extension_element b;
        uint64_t v_lo = 0; uint64_t v_hi = 0;
        uint64_t euclid_lo = 0; uint64_t euclid_hi = 0;

        for (int j = 0; j < 2; j++)
        {
            a = global::E.random();
            b = global::E.random();
            extension_repr tmp =
                global::E.fast_mul(a.get_repr(), b.get_repr());
            v_hi |= tmp.get_hi() << (32*j);
            v_lo |= tmp.get_lo() << (32*j);
            tmp = global::E.euclid_rem(tmp);
            euclid_hi |= tmp.get_hi() << (32*j);
            euclid_lo |= tmp.get_lo() << (32*j);
        }

        extension_repr v(v_hi, v_lo);
        extension_repr intel = global::E.packed_intel_rem(v);

        if (euclid_hi != intel.get_hi() || euclid_lo != intel.get_lo())
            err++;
    }
    end_test(err);
}

void Extension_test::test_packed_fast_mul()
{
    cout << "test packed fast mul: ";
    int err = 0;
    uint64_t mask = 0xFFFF;
    for (int i = 0; i < this->tests; i++)
    {
        extension_repr a;
        extension_repr b;
        uint64_t pa_lo = 0; uint64_t pa_hi = 0;
        uint64_t pb_lo = 0; uint64_t pb_hi = 0;
        uint64_t ref_lo = 0; uint64_t ref_hi = 0;

        for (int j = 0; j < 2; j++)
        {
            uint64_t r = global::randgen();

            a = {
                r & mask,
                (r >> 16) & mask
            };

            b = {
                (r >> 32) & mask,
                (r >> 48) & mask
            };

            pa_hi |= a.get_hi() << (32*j);
            pa_lo |= a.get_lo() << (32*j);

            pb_hi |= b.get_hi() << (32*j);
            pb_lo |= b.get_lo() << (32*j);

            extension_repr tmp = global::E.fast_mul(a, b);
            ref_hi |= tmp.get_hi() << (32*j);
            ref_lo |= tmp.get_lo() << (32*j);
        }

        extension_repr pa(pa_hi, pa_lo);
        extension_repr pb(pb_hi, pb_lo);
        extension_repr fast = global::E.packed_fast_mul(pa, pb);

        if (ref_hi != fast.get_hi() || ref_lo != fast.get_lo())
            err++;
    }
    end_test(err);
}
#endif
