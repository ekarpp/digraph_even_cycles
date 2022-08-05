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
GF2n global::F;
Extension global::E;
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
        cout << "-s for seed" << endl;
        return 0;
    }

    int opt;
    vector<vector<int>> graph;

    bool brute = false;
    bool duration = false;
    bool direct = false;
    uint64_t seed = time(nullptr);

    while ((opt = getopt(argc, argv, "utqbf:s:")) != -1)
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
        case '?':
            cout << "call with no arguments for help" << endl;
            return -1;
        }
    }

    cout << "seed: " << seed << endl;
    global::randgen.init(seed);

#if GF2_bits == 0
    int n = 24;
    uint64_t mod = util::irred_poly(n);
    global::F.init(n, mod);
    global::E.init(n, mod);
#else
    global::F.init();
    global::E.init();
#endif

    if (direct)
        util::direct_undirected(graph);

    Graph G(graph);
    Solver s;

    double start;
    double end;

    if (brute)
    {
        start = omp_get_wtime();
        int k = s.shortest_even_cycle_brute(G);
        end = omp_get_wtime();
        cout << k << endl;
    }
    else
    {
        start = omp_get_wtime();
        int k = s.shortest_even_cycle(G);
        end = omp_get_wtime();
        cout << k << endl;
    }

    if (duration) {
        double delta = end - start;
        cout << "computed graph of " << G.get_n() << " vertices in ";
        cout << delta << " seconds." << endl;
    }

    return 0;
}
