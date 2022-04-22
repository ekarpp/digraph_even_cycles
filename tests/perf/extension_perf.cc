/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <getopt.h>

#include "../../src/global.hh"
#include "../../src/extension.hh"
#include "../../src/util.hh"

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
    global::E.init();

    Extension_element a = global::E.random();
    Extension_element b = global::E.one();
    chrono::steady_clock::time_point start =
        chrono::steady_clock::now();
    for (uint64_t i = 0; i < t; i++)
        b *= a;
    chrono::steady_clock::time_point end =
        chrono::steady_clock::now();
    b.print();
    cout << t << " multiplications in time: " <<
        (chrono::duration_cast<chrono::microseconds>(end - start).count()) /1000000.0 << " s" << endl;
    return 0;
}
