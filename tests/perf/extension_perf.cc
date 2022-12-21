/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <omp.h>

#include "../../src/global.hh"
#include "../../src/extension.hh"
#include "../../src/gf.hh"

constexpr int WARMUP = 1 << 20;

#define BENCH_MUL(mul_func, last)                           \
{                                                           \
    GR_repr w = {0, 0};                              \
        _Pragma("omp parallel for")                         \
            for (uint64_t i = 0; i < WARMUP; i++)           \
                w = global::E->add(w, mul_func(a[i], b[i])); \
    start = omp_get_wtime();                                \
        _Pragma("omp parallel for")                         \
        for (uint64_t i = 0; i < t; i++)                    \
            a[i] = mul_func(a[i], b[i]);                    \
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
    GR_repr w = {0, 0};                              \
        _Pragma("omp parallel for")                         \
        for (uint64_t i = 0; i < WARMUP; i++)               \
            w = global::E->add(w, rem_func(a[i]));           \
    start = omp_get_wtime();                                \
        _Pragma("omp parallel for")                         \
            for (uint64_t i = 0; i < t; i++)                \
                a[i] = rem_func(a[i]);                      \
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
GR4_n *global::E;
GF2_n *global::F;
bool global::output = false;

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        cout << "-s $int for seed" << endl;
        cout << "-t $int for amount of tests, exponent for two" << endl;
        cout << "-n $int for size of finite field" << endl;
        cout << "-p $int for number of threads (default = all threads)" << endl;
        return 0;
    }

    uint64_t seed = time(nullptr);

    uint64_t t = 1;
    int exp = 1;
    int n = 16;
    int p = omp_get_num_threads();
    int opt;
    while ((opt = getopt(argc, argv, "s:t:n:p:")) != -1)
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
        case 'n':
            n = stoi(optarg);
            break;
        }
    }

    cout << "seed: " << seed << endl;
    global::randgen.init(seed);
    omp_set_num_threads(p);

    uint64_t mod;
    switch (n)
    {
    case 16:
        /* x^16 + x^5 + x^3 + x^2 +  1 */
        mod = 0x1002D;
        global::F = new GF2_16(16, mod);
        global::E = new GR4_16(16, mod);
        break;
    case 32:
        /* x^32 + x^7 + x^3 + x^2 + 1 */
        mod = 0x10000008D;
        global::F = new GF2_32(32, mod);
        global::E = new GR4_32(32, mod);
        break;
    default:
        mod = util::irred_poly(n);
        global::F = new GF2_n(n, mod);
        global::E = new GR4_n(n, mod);
        break;
    }

    vector<GR_repr> a(t);
    vector<GR_repr> b(t);
    vector<GR_repr> aa(t);

    double start;
    double end;

    start = omp_get_wtime();

    /* use as many random bits as possible to make init faster */
    const uint64_t REPEATS = (2*global::E->get_n()) / 64;
    for (uint64_t i = 0; i < t; i += REPEATS)
    {
        uint64_t ar = global::randgen();
        uint64_t br = global::randgen();
        for (uint64_t j = 0; j < REPEATS; j++)
        {
            uint64_t alo = ar >> (2 * global::E->get_n() * j);
            uint64_t ahi = ar >> (global::E->get_n() * (2 * j + 1));
            a[i+j] = {
                ahi & global::E->get_mask(),
                alo & global::E->get_mask()
            };
            aa[i+j] = a[i+j];

            uint64_t blo = br >> (2 * global::E->get_n() * j);
            uint64_t bhi = br >> (global::E->get_n() * (2 * j + 1));
            b[i+j] = {
                bhi & global::E->get_mask(),
                blo & global::E->get_mask()
            };
        }
    }
    end = omp_get_wtime();
    cout << "initialized in " << end - start << " s" << endl;

    double delta;
    double mhz;

    BENCH_MUL(global::E->ref_mul, 1);

    cout << t << " ref multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    BENCH_MUL(global::E->fast_mul, 1);

    cout << t << " fast multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    BENCH_MUL(global::E->kronecker_mul, 0);

    cout << t << " kronecker multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    cout << endl << endl;

    BENCH_REM(global::E->euclid_rem);

    cout << t << " euclid remainders in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    BENCH_REM(global::E->mont_rem);

    cout << t << " mont remainders in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    BENCH_REM(global::E->intel_rem);

    cout << t << " intel remainders in time " <<
        delta << " s or " << mhz << " Mhz" << endl;

    cout << endl << endl;

    return 0;
}
