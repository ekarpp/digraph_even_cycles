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

    /* real limit is 1 - (1 - 2^{-d})^n, we are rather conservative here. */
    const double error_lim = 1.0 / this->n;
    const double errorp = err * (1.0 / this->tests);

    if (errorp < error_lim)
        cout << "\033[32m";
    else
        cout << "\033[31m";
    cout << err << " out of " << this->tests << " failed (";
    cout << errorp * 100 << "%)" << "\033[0m" << endl;

    return (errorp >= error_lim);
}
