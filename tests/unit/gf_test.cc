/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <immintrin.h>

#include "gf_test.hh"
#include "../../src/gf.hh"
#include "../../src/global.hh"
#include "../../src/extension.hh"

using namespace std;

GF_test::GF_test()
{
    cout << "----------" << endl;
    cout << "TESTING GF" << endl;
    cout << "----------" << endl;

    this->run();
}

void GF_test::test_add_inverse()
{
    cout << "add inverse: ";
    int err = 0;
    uint64_t i = 0;
    while (i <= global::F.get_mask())
    {
        GF_element e(i);
        if (e + e != global::F.zero()
            || e - e != global::F.zero())
            err++;
        i++;
    }
    end_test(err);
}

void GF_test::test_associativity()
{
    cout << "test associativity: ";
    int err = 0;
    for (int i = 0; i < 10000; i++)
    {
        GF_element a = global::F.random();
        GF_element b = global::F.random();
        GF_element c = global::F.random();
        if (a*(b+c) != c*a + b*a)
            err++;
    }
    this->end_test(err);
}

void GF_test::test_mul_id()
{
    cout << "mul with id: ";
    int err = 0;
    uint64_t i = 0;
    while (i <= global::F.get_mask())
    {
        GF_element e(i);
        if (e * global::F.one() != e)
            err++;
        i++;
    }
    this->end_test(err);
}

void GF_test::test_mul_inverse()
{
    cout << "mul with inverse: ";
    int err = 0;
    uint64_t i = 1;
    while (i <= global::F.get_mask())
    {
        GF_element e(i);
        if (e / e != global::F.one())
            err++;
        i++;
    }
    this->end_test(err);
}

void GF_test::test_lift_project()
{
    cout << "lift project: ";
    int err = 0;
    uint64_t i = 0;
    while (i <= global::F.get_mask())
    {
        GF_element e(i);
        Extension_element b(i, global::randgen() & global::E.get_mask());
        Extension_element c(i, 0b0);
        if (e.lift().project() != e
            || b.project() != e
            || e.lift() != c)
            err++;
        i++;
    }
    this->end_test(err);
}

#if GF2_bits == 16
void GF_test::test_packed_rem()
{
    cout << "packed rem: ";
    int err = 0;
    uint64_t mask = global::F.get_mask() |
        (global::F.get_mask() << 32);
    for (int i = 0; i < this->tests; i++)
    {
        uint64_t a = global::randgen() & mask;
        uint64_t b = global::randgen() & mask;

        uint64_t ref = global::F.rem(
            global::F.clmul(
                (a >> 32) & global::F.get_mask(),
                (b >> 32) & global::F.get_mask()
            )
        );
        ref <<= 32;

        ref |= global::F.rem(
            global::F.clmul(
                a & global::F.get_mask(),
                b & global::F.get_mask()
            )
        );

        uint64_t r = global::F.packed_rem(
            global::F.packed_clmul(a, b)
        );

        if (ref != r)
            err++;
    }
    this->end_test(err);
}

void GF_test::test_wide_mul()
{
#define WIDTH 4
    cout << "wide mul: ";
    int err = 0;
    for (int i = 0; i < this->tests / WIDTH; i++)
    {
        uint64_t a[WIDTH];
        uint64_t b[WIDTH];
        int64_t prod[WIDTH];

        for (int j = 0; j < WIDTH; j++)
        {
            a[j] = global::randgen() & global::F.get_mask();
            b[j] = global::randgen() & global::F.get_mask();
            prod[j] = global::F.rem(
                global::F.clmul(a[j], b[j])
            );

            uint64_t ta = global::randgen() & global::F.get_mask();
            uint64_t tb = global::randgen() & global::F.get_mask();
            a[j] |= ta << 32;
            b[j] |= tb << 32;

            prod[j] |= global::F.rem(
                global::F.clmul(ta, tb)
            ) << 32;
        }

        __m256i aa = _mm256_set_epi64x(a[3], a[2], a[1], a[0]);
        __m256i bb = _mm256_set_epi64x(b[3], b[2], b[1], b[0]);
        __m256i pp = global::F.wide_mul(aa, bb);

        if (prod[0] != _mm256_extract_epi64(pp, 0))
            err++;

        if (prod[1] != _mm256_extract_epi64(pp, 1))
            err++;

        if (prod[2] != _mm256_extract_epi64(pp, 2))
            err++;

        if (prod[3] != _mm256_extract_epi64(pp, 3))
            err++;
    }
    this->end_test(err);
}
#endif
