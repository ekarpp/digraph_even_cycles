/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>

#include "solver_test.hh"
#include "../../src/graph.hh"
#include "../../src/global.hh"
#include "../../src/solver.hh"

using namespace std;

bool Solver_test::test_solver()
{
    cout << "solver random graph test: ";
    int err = 0;
    Solver s;
    for (int t = 0; t < this->tests; t++)
    {
        vector<vector<int>> adj(this->n, vector<int>());
        for (int u = 0; u < this->n; u++)
        {
            for (int v = 0; v < this->n; v++)
            {
                if (u == v)
                    continue;
                if ((global::randgen() & 0b11) == 0x0)
                    adj[u].push_back(v);
            }
        }

        Graph G(adj);

        if (s.shortest_even_cycle(G) != s.shortest_even_cycle_brute(G))
            err++;
    }

    /* TODO: add error margin */
    if (err == 0)
        cout << "\033[32m";
    else
        cout << "\033[31m";
    cout << err << " out of " << this->tests << " failed (";
    cout << (float) err / this->tests * 100 << "%)" << "\033[0m" << endl;

    return false;
}
