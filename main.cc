/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <cmath>

#include "src/gf.hh"
#include "src/global.hh"
#include "src/graph.hh"
#include "src/util.hh"
#include "src/fmatrix.hh"
#include "src/solver.hh"

using namespace std;

util::rand64bit global::randgen;
GF2n global::F;
Extension global::E;

int main(int argc, char **argv)
{
    uint64_t seed = time(nullptr);
    cout << "seed: " << seed << endl;
    global::randgen.init(seed);

    vector<vector<int>> vec = {
        {1},
        {0}
    };

    int d = 5 * ceil(log(vec.size()) / log(2));
    /* modz2 and gfmul only support upto 32 */
    if (d > 32)
        return -1;

    uint64_t poly = util::irred_poly(d);
    global::F.init(d, poly);
    global::E.init(d, poly);

    Graph G(vec);

    Solver s;

    cout << s.shortest_even_cycle(G) << endl;

    return 0;
}
