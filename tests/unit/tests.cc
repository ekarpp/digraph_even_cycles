/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <getopt.h>
#include <cstring>

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
    if (argc == 1 || (argc == 2 && strcmp(argv[1], "--help") == 0))
    {
        cout << "Usage: digraph-tests [-e] [-g] [-f] [-x] [-u] [-s] [-c] [-d <dimension>] [-n <degree>] [-t <repeats>] [-r <seed>]" << endl;
        cout << endl;
        cout << "Options:" << endl;
        cout << " -e\t execute extension ring tests" << endl;
        cout << " -g\t execute finite field tests" << endl;
        cout << " -f\t execute matrix tests with finite field elements" << endl;
        cout << " -x\t execute matrix tests with extension ring elements" << endl;
        cout << " -u\t execute utility functionality tests" << endl;
        cout << " -s\t execute solver tests" << endl;
        cout << " -c\t run geng tests. example pipe command: \"geng -q $n | directg -q | listg -aq | ./digraph-tests -c -n 16\"" << endl;
        cout << " -d\t dimension of square matrices for matrix tests" << endl;
        cout << " -n\t exponent for the underlying finite field with 3 <= n <= 32. optimized for n=16 or n=32." << endl;
        cout << " -t\t number of repeats for tests" << endl;
        cout << " -r\t seed fed to the random number generator" << endl;
        cout << " --help\t display usage information" << endl;
        cout << endl;
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

    bool failure = false;

    if (et)
    {
        GR_test e(tests);
        failure |= e.run();
    }
    if (gft)
    {
        GF_test f;
        failure |= f.run();
    }
    if (fmt)
    {
        FMatrix_test fm(tests);
        failure |= fm.run(dim);
    }
    if (ut)
    {
        Util_test u(tests);
        failure |= u.run(dim);
    }
    if (st)
    {
        Solver_test s(tests);
        failure |= s.run(dim);
    }
    if (emt)
    {
        EMatrix_test em(tests);
        failure |= em.run();
    }
    if (geng)
    {
        Geng_test g;
        failure |= g.run();
    }

    return failure;
}
