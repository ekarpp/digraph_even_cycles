/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <getopt.h>

#include "../../src/global.hh"
#include "../../src/util.hh"
#include "../../src/gf.hh"
#include "../../src/extension.hh"

#include "gf_test.hh"
#include "extension_test.hh"
#include "fmatrix_test.hh"
#include "util_test.hh"
#include "solver_test.hh"
#include "ematrix_test.hh"
#include "geng_test.hh"

util::rand64bit global::randgen;
GR4_n *global::E;
GF2_n *global::F;
bool global::output = false;

using namespace std;

int main(int argc, char** argv)
{
    if (argc == 1)
    {
        cout << "-e for extension tests" << endl;
        cout << "-g for GF tests" << endl;
        cout << "-f for FMatrix tests" << endl;
        cout << "-x for EMatrix tests" << endl;
        cout << "-u for util tests" << endl;
        cout << "-s for solver tests" << endl;
        cout << "-d $int dimension of matrix" << endl;
        cout << "-n $int degree of modulo polynomial. 16 and 32 optimized" << endl;
        cout << "-t $int how many times random tests are done" << endl;
        cout << "-c run geng test. \"geng -q $n | directg -q | listg -aq \" has to be piped to this." << endl;
        cout << "-r $int for seed to pseudo random generator" << endl;
        return 0;
    }

    bool et = false;
    bool gft = false;
    bool fmt = false;
    bool emt = false;
    bool ut = false;
    bool st = false;
    bool geng = false;
    int dim = 10;
    int tests = 10000;
    int opt;
    uint64_t seed = time(nullptr);

    int n = 16;
    while ((opt = getopt(argc, argv, "cxsuegfmn:d:t:r:")) != -1)
    {
        switch (opt)
        {
        case 'n':
            n = stoi(optarg);
            break;
        case 'r':
            seed = stoll(optarg);
            break;
        case 'c':
            geng = true;
            break;
        case 'x':
            emt = true;
            break;
        case 'u':
            ut = true;
            break;
        case 'd':
            dim = stoi(optarg);
            break;
        case 's':
            st = true;
            break;
        case 'e':
            et = true;
            break;
        case 'g':
            gft = true;
            break;
        case 'f':
            fmt = true;
            break;
        case 't':
            tests = stoi(optarg);
            break;
        }
    }

    cout << "seed: " << seed << endl;
    global::randgen.init(seed);

    if (n > 32 || n < 3)
    {
        cout << "please 3 <= n <= 32" << endl;
        return -1;
    }


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

    if (et)
        GR_test e(tests);
    if (gft)
        GF_test f;
    if (fmt)
        FMatrix_test fm(dim, tests);
    if (ut)
        Util_test u(dim, tests);
    if (st)
        Solver_test s(dim, tests);
    if (emt)
        EMatrix_test em(tests);
    if (geng)
        Geng_test g;

    return 0;
}
