/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <stdint.h>
#include <iostream>

#include "gf.hh"
#include "extension.hh"
#include "global.hh"

using namespace std;

/* GF */
GF_element GF2_n::zero() const
{
    return GF_element(0);
}

GF_element GF2_n::one() const
{
    return GF_element(1);
}

/* this can create zero, is it a problem? */
GF_element GF2_n::random() const
{
    return GF_element(global::randgen() & this->mask);
}

/* TODO: solve forward declaration issue and move to .hh */
GR_element GF_element::lift() const
{
    return GR_element(this->repr, 0b0);
}
