/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <omp.h>

#include "../../src/global.hh"
#include "../../src/extension.hh"
#include "../../src/gf.hh"

constexpr uint64_t WARMUP = 1 << 15;

enum Mul_enum { REF_MUL, FAST_MUL, KRONECKER_MUL };
enum Rem_enum { EUCLID_REM, INTEL_REM, MONT_REM };

using namespace std;

util::rand64bit global::randgen;
GR4_n *global::E;
GF2_n *global::F;
bool global::output = false;

template <Mul_enum M>
double bench_mul(
    vector<GR_repr> a,
    vector<GR_repr> b,
    vector<GR_repr> aa,
    bool last,
    uint64_t t)
{
    GR_repr w = {0, 0};
    uint64_t wup = (WARMUP > t) ? t : WARMUP;
    #pragma omp parallel for
    for (uint64_t i = 0; i < wup; i++)
        w = global::E->add(w, global::E->mul(a[i], b[i]));

    double start = omp_get_wtime();
    #pragma omp parallel for
    for (uint64_t i = 0; i < t; i++)
    {
        switch (M)
        {
        case REF_MUL:
            a[i] = global::E->ref_mul(a[i], b[i]);
            break;
        case FAST_MUL:
            a[i] = global::E->fast_mul(a[i], b[i]);
            break;
        case KRONECKER_MUL:
            a[i] = global::E->kronecker_mul(a[i], b[i]);
            break;
        }
    }
    double end = omp_get_wtime();

    if (start > end)
    {
        cout << a[time(nullptr) % t].hi << a[time(nullptr) % t].lo << endl;
        cout << w.hi << w.lo << endl;
    }

    for (uint64_t i = 0; i < t; i++) {
        if (last)
            aa[i] = a[i];
        else
            a[i] = aa[i];
    }

    return end - start;
}

template <Rem_enum R>
double bench_rem(
    vector<GR_repr> a,
    vector<GR_repr> aa,
    uint64_t t)
{
    GR_repr w = {0, 0};

    uint64_t wup = (WARMUP > t) ? t : WARMUP;
    #pragma omp parallel for
    for (uint64_t i = 0; i < wup; i++)
        w = global::E->add(w, global::E->rem(a[i]));

    double start = omp_get_wtime();
    #pragma omp parallel for
    for (uint64_t i = 0; i < t; i++)
    {
        switch (R)
        {
        case EUCLID_REM:
            a[i] = global::E->euclid_rem(a[i]);
            break;
        case INTEL_REM:
            a[i] = global::E->intel_rem(a[i]);
            break;
        case MONT_REM:
            a[i] = global::E->mont_rem(a[i]);
            break;
        }
    }
    double end = omp_get_wtime();

    if (start > end)
    {
        cout << a[time(nullptr) % t].hi << a[time(nullptr) % t].lo << endl;
        cout << w.hi << w.lo << endl;
    }

    for (uint64_t i = 0; i < t; i++)
        a[i] = aa[i];

    return end - start;
}

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

    /* use more random bits to make init faster */
    const uint64_t REPEATS = (global::E->get_n() > 16) ? 1 : 2;
    const uint64_t BITS = (global::E->get_n() > 16) ? 32 : 16;
    for (uint64_t i = 0; i < t; i += REPEATS)
    {
        uint64_t ar = global::randgen();
        uint64_t br = global::randgen();
        for (uint64_t j = 0; j < REPEATS; j++)
        {
            uint64_t alo = ar >> (2 * BITS * j);
            uint64_t ahi = ar >> (BITS * (2 * j + 1));
            a[i+j] = {
                ahi & global::E->get_mask(),
                alo & global::E->get_mask()
            };
            aa[i+j] = a[i+j];

            uint64_t blo = br >> (2 * BITS * j);
            uint64_t bhi = br >> (BITS * (2 * j + 1));
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

    delta = bench_mul<REF_MUL>(a, b, aa, 0, t);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " ref multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    delta = bench_mul<FAST_MUL>(a, b, aa, 0, t);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " fast multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    delta = bench_mul<KRONECKER_MUL>(a, b, aa, 1, t);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " kronecker multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    cout << endl;

    delta = bench_rem<EUCLID_REM>(a, aa, t);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " euclid remainders in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    delta = bench_rem<MONT_REM>(a, aa, t);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " mont remainders in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    delta = bench_rem<INTEL_REM>(a, aa, t);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " intel remainders in time " <<
        delta << " s or " << mhz << " Mhz" << endl;

    cout << endl;

    return 0;
}
