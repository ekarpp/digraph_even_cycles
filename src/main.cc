/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <cstring>

#include "gf.hh"
#include "global.hh"
#include "graph.hh"
#include "util.hh"
#include "fmatrix.hh"
#include "solver.hh"

using namespace std;

util::rand64bit global::randgen;
GF2_n *global::F;
GR4_n *global::E;
bool global::output = true;

bool parse_file(const string &fname, vector<vector<int>> &graph)
{
    ifstream file(fname);

    if (!file.is_open())
    {
        cout << "unable to open file: " << fname << endl;
        return false;
    }

    string line;
    int u = 0;
    while (getline(file, line))
    {
        /* comment */
        if (line[0] == '#')
            continue;
        istringstream iss(line);
        int v;
        vector<int> vec;
        while (iss >> v)
            vec.push_back(v);
        graph.push_back(vec);
        u++;
    }
    file.close();

    return true;
}

int main(const int argc, char **argv)
{
    if (argc == 1 || (argc == 2 && strcmp(argv[1], "--help") == 0))
    {
        cout << "Usage: digraph -f <file> [-b] [-q] [-t] [-u] [-n <field exponent>] [-s <seed>] [-p <threads>]" << endl;
        cout << endl;
        cout << "Options:" << endl;
        cout << " -f\t path to a graph file (custom syntax explained in readme.md)" << endl;
        cout << " -b\t use brute force solver (exponential complexity)" << endl;
        cout << " -q\t do not output progress of computation" << endl;
        cout << " -t\t output computation time" << endl;
        cout << " -u\t direct the input graph (random process)" << endl;
        cout << " -n\t exponent for the underlying finite field with 3 <= n <= 32. optimized for n=16 or n=32." << endl;
        cout << " -p\t number of threads (defaults to 1)" << endl;
        cout << " -s\t seed fed to the random number generator" << endl;
        cout << " --help\t display usage information" << endl;
        cout << endl;
        return 0;
    }

    int opt;
    vector<vector<int>> graph;

    bool brute = false;
    bool duration = false;
    bool direct = false;
    bool file_given = false;
    uint64_t seed = time(nullptr);
    int n = 16;
    int p = 1;

    while ((opt = getopt(argc, argv, "utqbf:s:n:p:")) != -1)
    {
        switch (opt)
        {
        case 'f':
            file_given = true;
            /* error during parsing */
            if (!parse_file(optarg, graph))
                return -1;
            break;
        case 's':
            seed = stoi(optarg);
            break;
        case 'u':
            direct = true;
            break;
        case 'b':
            brute = true;
            break;
        case 't':
            duration = true;
            break;
        case 'q':
            global::output = false;
            break;
        case 'n':
            n = stoi(optarg);
            break;
        case 'p':
            p = stoi(optarg);
            break;
        case '?':
            cout << "call with no arguments for help" << endl;
            return -1;
        }
    }

    if (!file_given)
    {
        cout << "-f is mandatory!" << endl;
        return -1;
    }

    if (n > 32 || n < 3)
    {
        cout << "please 3 <= n <= 32" << endl;
        return -1;
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

    if (direct)
        util::direct_undirected(graph);

    Graph G(graph);
    Solver s;

    const double start = omp_get_wtime();
    const int k = (brute)
        ? s.shortest_even_cycle_brute(G)
        : s.shortest_even_cycle(G);
    const double end = omp_get_wtime();

    cout << k << endl;

    if (duration) {
        const double delta = end - start;
        cout << "computed graph of " << G.get_n() << " vertices in ";
        cout << delta << " seconds." << endl;
    }

    return 0;
}
