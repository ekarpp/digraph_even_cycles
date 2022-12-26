/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "geng_test.hh"
#include "../../src/graph.hh"
#include "../../src/solver.hh"

using namespace std;

bool Geng_test::test_geng()
{
    string line = "";
    cin >> line;
    this->n = stoi(line);
    Solver s;

    int total = 0;
    int fail = 0;
    cout << "failed:" << endl;
    do
    {
        vector<vector<int>> graph;
        for (int i = 0; i < this->n; i++)
        {
            vector<int> vec;
            cin >> line;
            for (int j = 0; j < this->n; j++)
            {
                if (line[j] == '1')
                    vec.push_back(j);
            }
            graph.push_back(vec);

        }
        Graph G(graph);

        int algo = s.shortest_even_cycle(G);
        int ref = s.shortest_even_cycle_brute(G);
        if (algo != ref)
        {
            this->store_graph(graph, total);
            fail++;
        }
        total++;
        if (total % 1000 == 1)
        {
            cout << "\r" << fail << "/" << total;
            cout << " (" << ((float) fail) / total * 100 << "%)";
            cout << flush;
        }
    } while (cin >> line);
    cout << "\r" << fail << "/" << total;
    cout << " (" << ((float) fail) / total * 100 << "%)";
    cout << endl;

    /* TODO: add error margin */
    return false;
}

void Geng_test::store_graph(const vector<vector<int>> &g, int id)
{
    ofstream file("geng-fail/" + to_string(this->n) + "/" + to_string(id));

    /* just return if cant open */
    if (!file.is_open())
        return;

    for (int i = 0; i < this->n; i++)
    {
        for (uint j = 0; j < g[i].size(); j++)
            file << j << " ";
        file << "\n";
    }

    file.close();
}
