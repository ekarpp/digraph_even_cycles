/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef BITVECTOR_H
#define BITVECTOR_H

#include <stdint.h>
#include <immintrin.h>

struct uint128_t
{
    unsigned long long words[2];
};

struct uint256_t
{
    unsigned long long words[4];
};

struct uint512_t
{
    unsigned long long words[8];
};

struct uint576_t
{
    unsigned long long words[9];
};

#define ADD_WORDS(n)                                                           \
{                                                                              \
    char carry = 0;                                                            \
    for (int i = 0; i < n; i++)                                                \
        carry = _addcarry_u64(carry, a.words[i], b.words[i], &(sum.words[i])); \
}

namespace bit
{

    inline uint256_t add_256bit(uint256_t a, uint256_t b)
    {
        uint256_t sum;
        ADD_WORDS(4);
        return sum;
    }

    inline uint512_t add_512bit(uint512_t a, uint512_t b)
    {
        uint512_t sum;
        ADD_WORDS(8);
        return sum;
    }

    inline uint576_t add_576bit(uint576_t a, uint576_t b)
    {
        uint576_t sum;
        ADD_WORDS(9);
        return sum;
    }

    inline uint576_t pad_words(uint512_t a, int n)
    {
        uint576_t padded;
        for (int i = 0; i < n; i++)
            padded.words[i] = 0;

        for (int i = n; i < 9; i++)
            padded.words[i] = a.words[i - n];

        return padded;
    }

    inline uint576_t widen_512bits(uint512_t a)
    {
        uint576_t wide;
        for (int i = 0; i < 8; i++)
            wide.words[i] = a.words[i];

        wide.words[8] = 0;

        return wide;
    }

    /* pos <= 64 */
    inline uint576_t lshift_512bit(uint512_t a, int pos)
    {
        uint576_t shift;
        shift.words[0] = a.words[0] << pos;
        for (int i = 1; i < 8; i++)
        {
            shift.words[i] = a.words[i] << pos;
            shift.words[i] |= a.words[i-1] >> (64 - pos);
        }
        shift.words[8] = a.words[7] >> (64 - pos);
        return shift;
    }

    inline char add_128bit_carry(uint128_t a, uint128_t b, unsigned long long *sum)
    {
        char carry = 0;
        carry = _addcarry_u64(carry, a.words[0], b.words[0], sum);
        carry = _addcarry_u64(carry, a.words[1], b.words[1], sum + 1);
        return carry;
    }

    inline uint256_t mul_128bit(uint128_t a, uint128_t b)
    {
        uint256_t prod;
        prod.words[0] = _mulx_u64(a.words[0], b.words[0], prod.words + 1);
        prod.words[2] = _mulx_u64(a.words[1], b.words[1], prod.words + 3);

        uint128_t ahbl, albh;
        ahbl.words[0] = _mulx_u64(a.words[1], b.words[0], ahbl.words + 1);
        albh.words[0] = _mulx_u64(a.words[0], b.words[1], albh.words + 1);

        unsigned long long sum[2];
        /* add carry to hi of product */
        prod.words[3] += add_128bit_carry(ahbl, albh, sum);
        prod.words[3] += add_128bit_carry(
            { sum[0], sum[1] },
            { prod.words[1], prod.words[2] },
            prod.words + 1
        );

        return prod;
    }

    inline uint512_t mul_256bit_64bit(uint256_t a, uint64_t b)
    {
        /*
          (a0 +  a1 x^64  + a2 x^128 + a3 x^192)b =
          a0 b + a1 b x^64 + a2 b x^128 + a3 b x^192
         */
        uint512_t prod_even;
        prod_even.words[0] = _mulx_u64(a.words[0], b, prod_even.words + 1);
        prod_even.words[2] = _mulx_u64(a.words[2], b, prod_even.words + 3);

        prod_even.words[4] = 0; prod_even.words[5] = 0;
        prod_even.words[6] = 0; prod_even.words[7] = 0;

        uint512_t prod_odd;
        prod_odd.words[1] = _mulx_u64(a.words[1], b, prod_odd.words + 2);
        prod_odd.words[3] = _mulx_u64(a.words[3], b, prod_odd.words + 4);

        prod_odd.words[0] = 0; prod_odd.words[5] = 0;
        prod_odd.words[6] = 0; prod_odd.words[7] = 0;

        return add_512bit(prod_even, prod_odd);
    }

    inline uint512_t mul_256bit(uint256_t a, uint256_t b)
    {
        uint256_t ahbh = bit::mul_128bit(
            { a.words[2], a.words[3] },
            { b.words[2], b.words[3] }
        );
        uint256_t ahbl = bit::mul_128bit(
            { a.words[2], a.words[3] },
            { b.words[0], b.words[1] }
        );
        uint256_t albh = bit::mul_128bit(
            { a.words[0], a.words[1] },
            { b.words[2], b.words[3] }
        );
        uint256_t albl = bit::mul_128bit(
            { a.words[0], a.words[1] },
            { b.words[0], b.words[1] }
        );

        uint256_t mid;
        unsigned char carry = 0;
        long long unsigned int sum[2];
        carry = add_128bit_carry(
            { albl.words[2], albl.words[3] },
            { ahbl.words[0], ahbl.words[1] },
            sum
        );
        carry += add_128bit_carry(
            { sum[0], sum[1] },
            { albh.words[0], albh.words[1] },
            sum
        );
        mid.words[0] = sum[0];
        mid.words[1] = sum[1];
        mid.words[2] = carry;

        uint512_t ret;
        ret.words[0] = albl.words[0];
        ret.words[1] = albl.words[1];
        ret.words[2] = mid.words[0];
        ret.words[3] = mid.words[1];

        sum[0] = 0;
        sum[1] = 0;
        carry = add_128bit_carry(
            { ahbh.words[0], ahbh.words[1] },
            { ahbl.words[2], ahbl.words[3] },
            sum
        );
        carry += add_128bit_carry(
            { sum[0], sum[1] },
            { albh.words[2], albh.words[3] },
            sum
        );
        carry += add_128bit_carry(
            { sum[0], sum[1] },
            { mid.words[2], 0 },
            sum
        );

        ret.words[4] = sum[0];
        ret.words[5] = sum[1];

        sum[0] = 0;
        sum[1] = 0;
        carry = add_128bit_carry(
            { carry, 0 },
            { ahbh.words[2], ahbh.words[3] },
            sum
        );

        ret.words[6] = sum[0];
        ret.words[7] = sum[1];

        return ret;
    }
}
#endif
