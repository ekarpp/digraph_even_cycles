/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef SOLVER_H
#define SOLVER_H

#include "graph.hh"

class Solver
{
public:
    Solver() {}

    int shortest_even_cycle(Graph &G) const;

    int shortest_even_cycle_brute(const Graph &G) const;
};

#endif
