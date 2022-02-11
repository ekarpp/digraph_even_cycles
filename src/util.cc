#include <bitset>
#include <cmath>
#include <set>

#include "util.hh"

using namespace std;

/* Ben-Or's irreducible polynomial generator.
 * returns irreducible polynomial of degree deg
 * in Z2[x] encoded as a bitstring
 */
unsigned long long ben_or(int deg)
{
    bitset<64> p;
    p[deg] = true;
    srand(1);

    while (true)
    {
        int i;
        for (i = 0; i < deg; i++)
            p[i] = rand() & 1;
        for (i = 1; i <= deg >> 1; i++)
        {
            if (!gcd1(i, p))
                i = deg;
        }

        if (i < deg)
            break;
    }

    return p.to_ullong();
}

/* is gcd of x^(2^i) - x and p one (in Z2[x])*/
bool gcd1(int i, bitset<64> p)
{
    /* use set to represent polynomials for conveninece of
     * the methods and to save space as we have degree 2^i */
    set<long long> r;
    r.insert(1 << i);
    r.insert(1);

    set<long long> rn;
    for (int j = 0; j < 64; j++)
    {
        if (p[j])
            rn.insert(j);
    }

    if (*r.rbegin() < *rn.rbegin())
        r.swap(rn);

    /* standard Euclid's algo with Euclidean division */
    while (rn.size() != 0)
    {
        long long deg_r = *r.rbegin();
        long long deg_rn = *rn.rbegin();

        while (deg_r >= deg_rn)
        {
            set<long long>::reverse_iterator it = rn.rbegin();
            for ( ; it != rn.rend(); it++)
            {
                long long id = *it + deg_r - deg_rn;
                if (r.count(id) == 1)
                    r.erase(id);
                else
                    r.insert(id);
            }
            if (r.empty())
                deg_r = -1;
            else
                deg_r = *r.rbegin();

        }

        rn.swap(r);
    }

    return r.size() == 1 && r.count(0) == 1;
}
