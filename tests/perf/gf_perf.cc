/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <omp.h>

#include "../../src/global.hh"
#include "../../src/gf.hh"
#include "../../src/extension.hh"

typedef long long int long4_t __attribute__ ((vector_size (32)));

constexpr int VECTOR_N = 8;

using namespace std;

util::rand64bit global::randgen;
GF2_n *global::F;
GR4_n *global::E;
bool global::output = false;

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        cout << "-s $int for seed" << endl;
        cout << "-t $int for amount of tests" << endl;
        cout << "-n $int for size of finite field" << endl;
        return 0;
    }
    uint64_t seed = time(nullptr);

    uint64_t t = 1;
    int n = 16;
    int opt;
    while ((opt = getopt(argc, argv, "s:t:n:")) != -1)
    {
        switch (opt)
        {
        case 'n':
            n = stoi(optarg);
            break;
        case 's':
            seed = stoi(optarg);
            break;
        case 't':
            t = stoi(optarg);
            break;
        }
    }

    cout << "seed: " << seed << endl;
    global::randgen.init(seed);

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

    vector<uint64_t> a(t);
    vector<uint64_t> b(t);
    vector<uint64_t> p(t);
    vector<uint64_t> r(t);

    vector<GF_element> aa(t);
    vector<GF_element> bb(t);

    for (uint64_t i = 0; i < t; i++)
    {
        a[i] = global::randgen() & global::F->get_mask();
        b[i] = global::randgen() & global::F->get_mask();
        aa[i] = util::GF_random();
        bb[i] = util::GF_random();
    }

    double start;
    double end;
    double delta;
    double mhz;

    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
        aa[i] *= bb[i];
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " multiplications (whole) in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
        p[i] = global::F->clmul(a[i], b[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
        r[i] = global::F->rem(p[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " remainder in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
        r[i] = global::F->ext_euclid(r[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " inversion in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;


    if (global::F->get_n() != 16)
        return 0;

    vector<long4_t> av(t);
    vector<long4_t> bv(t);
    vector<long4_t> pv(t);
    for (uint64_t i = 0; i < t; i++)
    {
        av[i] = _mm256_set_epi64x(
            global::randgen(),
            global::randgen(),
            global::randgen(),
            global::randgen()
        );
        bv[i] = _mm256_set_epi64x(
            global::randgen(),
            global::randgen(),
            global::randgen(),
            global::randgen()
        );
    }

    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
        pv[i] = global::F->wide_mul(av[i], bv[i]);
    end = omp_get_wtime();
    delta = end - start;
    mhz = VECTOR_N*t / delta;
    mhz /= 1e6;

    cout << VECTOR_N*t << " muls with wide mul in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    return 0;
}
