/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <omp.h>

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

bool parse_file(string fname, vector<vector<int>> &graph)
{
    string line;
    ifstream file(fname);

    if (!file.is_open())
    {
        cout << "unable to open file: " << fname << endl;
        return false;
    }


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

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        cout << "provide path to graph file with -f" << endl;
        cout << "line i (starting at zero) in graph file" << endl;
        cout << "tells to which nodes there is an edge to" << endl;
        cout << "-b to use the brute force solver instead (SLOW)" << endl;
        cout << "-q for no progress output from solver" << endl;
        cout << "-t to output time spent computing" << endl;
        cout << "-u if the input graph is undirected, directs it randomly" << endl;
        cout << "-n finite field size. 16 and 32 optimized. has to be <= 32" << endl;
        cout << "-p for number of threads (defaults to all)" << endl;
        cout << "-s for seed" << endl;
        return 0;
    }

    int opt;
    vector<vector<int>> graph;

    bool brute = false;
    bool duration = false;
    bool direct = false;
    uint64_t seed = time(nullptr);
    int n = 16;
    int p = omp_get_num_threads();

    while ((opt = getopt(argc, argv, "utqbf:s:n:p:")) != -1)
    {
        switch (opt)
        {
        case 'f':
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

    if (n > 32)
    {
        cout << "please n <= 32" << endl;
        return -1;
    }

    cout << "seed: " << seed << endl;
    global::randgen.init(seed);

    uint64_t mod;

    switch (n)
    {
    case 16:
        /* x^16 + x^5 + x^3 + x^2 +  1 */
        mod = 0x1002D;
        global::F = new GF2_16(mod, 16);
        global::E = new GR4_16(mod, 16);
        break;
    case 32:
        /* x^32 + x^7 + x^3 + x^2 + 1 */
        mod = 0x10000008D;
        global::F = new GF2_32(mod, 32);
        global::E = new GR4_32(mod, 32);
        break;
    default:
        mod = util::irred_poly(n);
        global::F = new GF2_n(mod, n);
        global::E = new GR4_n(mod, n);
        break;
    }

    if (direct)
        util::direct_undirected(graph);

    Graph G(graph);
    Solver s;

    double start;
    double end;

    int k;

    omp_set_num_threads(p);

    if (brute)
    {
        start = omp_get_wtime();
        k = s.shortest_even_cycle_brute(G);
        end = omp_get_wtime();
    }
    else
    {
        start = omp_get_wtime();
        k = s.shortest_even_cycle(G);
        end = omp_get_wtime();
    }

    cout << k << endl;

    if (duration) {
        double delta = end - start;
        cout << "computed graph of " << G.get_n() << " vertices in ";
        cout << delta << " seconds." << endl;
    }

    delete global::F;
    delete global::E;

    return 0;
}
