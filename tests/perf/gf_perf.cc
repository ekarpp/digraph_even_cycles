/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <omp.h>

#include "../../src/global.hh"
#include "../../src/gf.hh"
#include "../../src/extension.hh"

using namespace std;

util::rand64bit global::randgen;
GF2n global::F;
Extension global::E;
bool global::output = false;

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        cout << "-s $int for seed" << endl;
        cout << "-t $int for amount of tests" << endl;
        return 0;
    }
    uint64_t seed = time(nullptr);

    uint64_t t = 1;
    int opt;
    while ((opt = getopt(argc, argv, "s:t:")) != -1)
    {
        switch (opt)
        {
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
#if GF2_bits == 0
    int n = 10;
    uint64_t mod = util::irred_poly(n);
    global::F.init(n, mod);
#else
    global::F.init();
#endif

    vector<uint64_t> a(t);
    vector<uint64_t> b(t);
    vector<uint64_t> p(t);
    vector<uint64_t> r(t);

    vector<GF_element> aa(t);
    vector<GF_element> bb(t);

    for (uint64_t i = 0; i < t; i++)
    {
        a[i] = global::randgen() & global::F.get_mask();
        b[i] = global::randgen() & global::F.get_mask();
        aa[i] = global::F.random();
        bb[i] = global::F.random();
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
        p[i] = global::F.clmul(a[i], b[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " multiplications in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
        r[i] = global::F.rem(p[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " remainder in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    start = omp_get_wtime();
    for (uint64_t i = 0; i < t; i++)
        r[i] = global::F.ext_euclid(r[i]);
    end = omp_get_wtime();
    delta = (end - start);
    mhz = t / delta;
    mhz /= 1e6;

    cout << t << " inversion in time: " <<
        delta << " s or " << mhz << " Mhz" << endl;

    return 0;
}
