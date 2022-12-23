/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <valarray>

#include "global.hh"
#include "graph.hh"
#include "gf.hh"

using namespace std;

Graph::Graph(vector<vector<int>> &adjacency_list): A(adjacency_list.size())
{
    this->n = adjacency_list.size();
    this->adj = adjacency_list;
    this->sample_adjacency();

    if (global::output)
        cout << "created graph of " << this->n << " vertices:" << endl;
}

/* samples the adjacency matrix with random edge weights from F
 * also creates a loop at each vertex
 */
void Graph::sample_adjacency()
{
    valarray<GF_element> m(util::GF_zero(), this->n * this->n);

    for (int u = 0; u < this->n; u++)
    {
        /* loop at each vertex */
        m[u*this->n + u] = util::GF_random();
        for (uint i = 0; i < this->adj[u].size(); i++)
        {
            const int v = this->adj[u][i];
            m[u*this->n + v] = util::GF_random();
        }
    }

    this->A.copy(FMatrix(this->n, m));
    return;
}

/* goes through all cycles that contain vertex start
 * and updates len accordingly.
 * len contains the length of the shortest found so far */
int Graph::dfs_cycle(const int start,
                     const int depth,
                     const int v,
                     vector<bool> &visited,
                     int len) const
{
    visited[v] = true;
    vector<int> nbors = this->adj[v];
    for (uint i = 0; i < nbors.size(); i++)
    {
        int u = nbors[i];
        if (!visited[u])
            len = this->dfs_cycle(start, depth + 1, u, visited, len);
        else if (u == start && depth % 2 == 0 && depth < len)
            len = depth;
    }
    visited[v] = false;
    return len;
}
