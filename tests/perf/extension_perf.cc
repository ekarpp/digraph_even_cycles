/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <omp.h>

#include "../../src/global.hh"
#include "../../src/extension.hh"
#include "../../src/gf.hh"

#ifndef PARALLEL
#define PARALLEL 0
#endif

#define WARMUP (1 << 20)

#define BENCH_MUL(mul_func, last)                           \
{                                                           \
    extension_repr w = {0, 0};                              \
    if (PARALLEL) {                                         \
        _Pragma("omp parallel for")                         \
            for (uint64_t i = 0; i < WARMUP; i++)           \
            w = global::E.add(w, mul_func(a[i], b[i]));     \
    } else {                                                \
        for (uint64_t i = 0; i < WARMUP; i++)               \
            w = global::E.add(w, mul_func(a[i], b[i]));     \
    }                                                       \
    start = omp_get_wtime();                                \
    if (PARALLEL) {                                         \
        _Pragma("omp parallel for")                         \
        for (uint64_t i = 0; i < t; i++)                    \
            a[i] = mul_func(a[i], b[i]);                    \
    } else {                                                \
        for (uint64_t i = 0; i < t; i++)                    \
            a[i] = mul_func(a[i], b[i]);                    \
    }                                                       \
    end = omp_get_wtime();                                  \
    if (exp == 0x6fabc73829101) {                           \
        cout << a[exp].hi << a[exp].lo << endl;             \
        cout << w.hi << w.lo << endl;                       \
    }                                                       \
    for (uint64_t i = 0; i < t; i++) {                      \
        if (last)                                           \
            a[i] = aa[i];                                   \
        else                                                \
            aa[i] = a[i];                                   \
    }                                                       \
    delta = (end - start);                                  \
    mhz = t / delta;                                        \
    mhz /= 1e6;                                             \
}

#define BENCH_REM(rem_func)                                 \
{                                                           \
    extension_repr w = {0, 0};                              \
    if (PARALLEL) {                                         \
        _Pragma("omp parallel for")                         \
        for (uint64_t i = 0; i < WARMUP; i++)               \
            w = global::E.add(w, rem_func(a[i]));           \
    } else {                                                \
        for (uint64_t i = 0; i < WARMUP; i++)               \
            w = global::E.add(w, rem_func(a[i]));           \
    }                                                       \
    start = omp_get_wtime();                                \
    if (PARALLEL) {                                         \
        _Pragma("omp parallel for")                         \
            for (uint64_t i = 0; i < t; i++)                \
                a[i] = rem_func(a[i]);                      \
    } else {                                                \
        for (uint64_t i = 0; i < t; i++)                    \
            a[i] = rem_func(a[i]);                          \
    }                                                       \
    end = omp_get_wtime();                                  \
    if (exp == 0x6fabc73829101) {                           \
        cout << a[exp].hi << a[exp].lo << endl;             \
        cout << w.hi << w.lo << endl;                       \
    }                                                       \
    for (uint64_t i = 0; i < t; i++)                        \
        a[i] = aa[i];                                       \
    delta = (end - start);                                  \
    mhz = t / delta;                                        \
    mhz /= 1e6;                                             \
}

using namespace std;

util::rand64bit global::randgen;
Extension global::E;
GF2n global::F;
bool global::output = false;

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        cout << "-s $int for seed" << endl;
        cout << "-t $int for amount of tests, exponent for two" << endl;
        return 0;
    }

    uint64_t seed = time(nullptr);

    uint64_t t = 1;
    int exp = 1;
    int opt;
    while ((opt = getopt(argc, argv, "s:t:")) != -1)
    {
        switch (opt)
        {
        case 's':
            seed = stoi(optarg);
            break;
        case 't':
            exp = stoi(optarg);
            t <<= exp;
            break;
        }
    }

    cout << "seed: " << seed << endl;
    global::randgen.init(seed);
#if GF2_bits == 0
    int n = 24;
    cout << "E^" << n << endl;
    uint64_t mod = util::irred_poly(n);
    global::E.init(n, mod);
#else
    global::E.init();
#endif

    vector<extension_repr> a(t);
    vector<extension_repr> b(t);
    vector<extension_repr> aa(t);

    double start;
    double end;

    start = omp_get_wtime();
#if GF2_bits == 16
    #define REPEATS 2
#else
    #define REPEATS 1
#endif

    for (uint64_t i = 0; i < t; i += REPEATS)
    {
        uint64_t ar = global::randgen();
        uint64_t br = global::randgen();
        for (uint64_t j = 0; j < REPEATS; j++)
        {
            uint64_t alo = ar >> (2*GF2_bits*j);
            uint64_t ahi = ar >> (GF2_bits*(2*j+1));
            a[i+j] = {
                ahi & global::E.get_mask(),
                alo & global::E.get_mask()
            };
            aa[i+j] = a[i+j];

            uint64_t blo = br >> (2*GF2_bits*j);
            uint64_t bhi = br >> (GF2_bits*(2*j+1));
            b[i+j] = {
                bhi & global::E.get_mask(),
                blo & global::E.get_mask()
            };
        }
    }
    end = omp_get_wtime();
    cout << "initialized in " << end - start << " s" << endl;

    double delta;
    double mhz;

    BENCH_MUL(global::E.ref_mul, 1);

    cout << t << " ref multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    BENCH_MUL(global::E.fast_mul, 1);

    cout << t << " fast multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    BENCH_MUL(global::E.kronecker_mul, 0);

    cout << t << " kronecker multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    cout << endl << endl;

    BENCH_REM(global::E.euclid_rem);

    cout << t << " euclid remainders in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    BENCH_REM(global::E.mont_rem);

    cout << t << " mont remainders in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    BENCH_REM(global::E.intel_rem);

    cout << t << " intel remainders in time " <<
        delta << " s or " << mhz << " Mhz" << endl;

    cout << endl << endl;

#if GF2_bits == 16
    uint64_t pack_mask = 0xFFFF | (0xFFFFull << 32);
    t /= 2;
    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
    {
        uint64_t ar = global::randgen();
        uint64_t br = global::randgen();

        a[i] = {
            ar & pack_mask,
            (ar >> 16) & pack_mask
        };
        aa[i] = a[i];
        b[i] = {
            br & pack_mask,
            (br >> 16) & pack_mask
        };
    }
    end = omp_get_wtime();
    cout << "initialized in " << end - start << " s" << endl;

    BENCH_MUL(global::E.packed_fast_mul, 0);

    cout << 2*t << " pack fast multiplications in time: " <<
        delta << " s or " << 2*mhz << " Mhz" << endl;

    BENCH_REM(global::E.packed_intel_rem);

    cout << 2*t << " pack intel remainders in time " <<
        delta << " s or " << 2*mhz << " Mhz" << endl;
#endif

    return 0;
}
