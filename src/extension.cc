/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>

#include "gf.hh"
#include "global.hh"
#include "extension.hh"

using namespace std;

/* Extension */
GR_element GR4_n::zero() const
{
    return GR_element(0b0, 0b0);
}

GR_element GR4_n::one() const
{
    return GR_element(0b1, 0b0);
}

GR_element GR4_n::random() const
{
    return GR_element(
        global::randgen() & this->mask,
        global::randgen() & this->mask
    );
}

/* Extension element */
GF_element GR_element::project() const
{
    return GF_element(this->repr.lo);
}
