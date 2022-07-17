/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef P_FMATRIX_H
#define P_FMATRIX_H

#include <valarray>
#include <vector>
#include <immintrin.h>

#include "gf.hh"
#include "global.hh"
#include "fmatrix.hh"

typedef long long int long4_t __attribute__ ((vector_size (32)));

#define VECTOR_N 8

/*
  _mm256_permutevar8x32_epi32(coeffs[col], idx),                        \
  _mm256_permutevar8x32_epi32(coeffs[col + 1], idx),                    \
 */

#define COEFF_LOOP(index)                                               \
    {                                                                   \
    for (int col = 0; col < this->cols - 1; col++)                      \
        coeffs[col] = _mm256_blend_epi32(                               \
            coeffs[col],                                                \
            coeffs[col + 1],                                            \
            0xFF >> index                                               \
        );                                                              \
    }

#define DET_LOOP(index)                                         \
    {                                                           \
        int r0 = VECTOR_N*col + index;                          \
        long4_t mx = this->get(r0, col);                        \
        int mxi = r0;                                           \
        int cmpmsk = 0b11 << (4*(VECTOR_N - 1 - index));        \
        for (int row = r0 + 1; row < this->rows; row++)         \
        {                                                       \
            int cmp = _mm256_movemask_epi8(                     \
                _mm256_cmpgt_epi32(                             \
                    this->get(row, col),                        \
                    mx                                          \
                    )                                           \
                );                                              \
            if (cmp & cmpmsk)                                   \
            {                                                   \
                mx = this->get(row,col);                        \
                mxi = row;                                      \
            }                                                   \
        }                                                       \
        uint64_t mx_ext =                                       \
            _mm256_extract_epi32(mx, VECTOR_N - 1 - index);     \
        if (mx_ext == 0)                                        \
            return global::F.zero();                            \
        if (mxi != r0)                                          \
            this->swap_rows(mxi, r0);                           \
        /* vectorize? */                                        \
        det = global::F.rem(                                    \
            global::F.clmul(det, mx_ext)                        \
        );                                                      \
        mx_ext = global::F.ext_euclid(mx_ext);                  \
        this->mul_row(r0, mx_ext);                              \
        /* vectorize end? */                                    \
        char mask = VECTOR_N - 1 - index;                       \
        long4_t idx = _mm256_set_epi32(                         \
            mask,                                               \
            mask,                                               \
            mask,                                               \
            mask,                                               \
            mask,                                               \
            mask,                                               \
            mask,                                               \
            mask                                                \
        );                                                      \
        for (int row = r0 + 1; row < this->rows; row++)         \
        {                                                       \
            long4_t val = _mm256_permutevar8x32_epi32(          \
                this->get(row, col),                            \
                idx                                             \
            );                                                  \
            this->row_op(r0, row, val);                         \
        }                                                       \
    }

class Packed_FMatrix
{
private:
    int rows;
    int cols;
    // original matrix n moduloe VECTOR_N
    int nmod;
    std::vector<long4_t> m;

    long4_t get(int row, int col) const
    {
        return this->m[row*this->cols + col];
    }

    void set(int row, int col, long4_t v)
    {
        this->m[row*this->cols + col] = v;
    }

public:
    Packed_FMatrix(
        const FMatrix &matrix
    )
    {
        this->nmod = matrix.get_n() % VECTOR_N;
        this->rows = matrix.get_n();
        if (this->rows % VECTOR_N)
            this->rows += VECTOR_N - (matrix.get_n() % VECTOR_N);
        this->cols = this->rows / VECTOR_N;

        this->m.resize(this->rows * this->cols);

        for (int r = 0; r < matrix.get_n(); r++)
        {
            for (int c = 0; c < matrix.get_n() / VECTOR_N; c++)
                this->set(r, c,
                          _mm256_set_epi64x(
                              matrix(r, VECTOR_N*c + 0).get_repr() << 32
                                  | matrix(r, VECTOR_N*c + 1).get_repr(),
                              matrix(r, VECTOR_N*c + 2).get_repr() << 32
                                  | matrix(r, VECTOR_N*c + 3).get_repr(),
                              matrix(r, VECTOR_N*c + 4).get_repr() << 32
                                  | matrix(r, VECTOR_N*c + 5).get_repr(),
                              matrix(r, VECTOR_N*c + 6).get_repr() << 32
                                  | matrix(r, VECTOR_N*c + 7).get_repr()
                          )
                    );
            if (this->nmod)
            {
                int c = this->cols - 1;
                uint64_t elems[4];
                elems[0] = 0; elems[1] = 0; elems[2] = 0; elems[3] = 0;
                for (int i = 0; i < this->nmod; i++)
                    elems[i/2] |= matrix(r, VECTOR_N*c + i).get_repr() << (32*(1 - i%2));

                this->set(r, c, _mm256_set_epi64x(
                              elems[0],
                              elems[1],
                              elems[2],
                              elems[3]
                          )
                );
            }
        }
        for (int r = matrix.get_n(); r < this->rows; r++)
        {
            for (int c = 0; c < this->cols - 1; c++)
                this->set(r, c, _mm256_setzero_si256());

            /* lazy.... */
            switch (r % VECTOR_N)
            {
            case 1:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              1, 0, 0, 0
                         )
                );
                break;
            case 2:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              0, 1ull << 32, 0, 0
                         )
                );
                break;
            case 3:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              0, 1, 0, 0
                         )
                );
                break;
            case 4:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              0, 0, 1ull << 32, 0
                         )
                );
                break;
            case 5:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              0, 0, 1, 0
                         )
                );
                break;
            case 6:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              0, 0, 0, 1ull << 32
                         )
                );
                break;
            case 7:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              0, 0, 0, 1
                         )
                );
                break;
            }

        }
    }

    void mul_gamma(int r1, int r2, const GF_element &gamma)
    {
        /* here we do r1 first left to right and save the auxiliary gamma vectors.
         * then we permute the gamma vectors as required (note this->nmod here),
         * and do r2 left to right */
        /* do gamma multiplication during initialization?? */
        long4_t pac_gamma = _mm256_set_epi64x(
            gamma.get_repr() << 32 | gamma.get_repr(),
            gamma.get_repr() << 32 | gamma.get_repr(),
            gamma.get_repr() << 32 | gamma.get_repr(),
            gamma.get_repr() << 32 | gamma.get_repr()
        );
        // pac_gamma = [gamma^VECTOR_N]
        pac_gamma = global::F.wide_mul(pac_gamma, pac_gamma);
        pac_gamma = global::F.wide_mul(pac_gamma, pac_gamma);
        pac_gamma = global::F.wide_mul(pac_gamma, pac_gamma);

        uint64_t elems[4];
        uint64_t g = 1ull;
        for (int i = 0; i < 4; i++)
        {
            elems[i] = g << 32;
            g = global::F.rem(global::F.clmul(g, gamma.get_repr()));
            elems[i] |= g;
            g = global::F.rem(global::F.clmul(g, gamma.get_repr()));
        }
        long4_t prod = _mm256_set_epi64x(
            elems[0],
            elems[1],
            elems[2],
            elems[3]
        );
        std::vector<long4_t> coeffs(this->cols);
        const long4_t idx = _mm256_set_epi32(
            0b000,
            0b001,
            0b010,
            0b011,
            0b100,
            0b101,
            0b110,
            0b111
        );
        /* first do left to right r1 */
        for (int col = 0; col < this->cols; col++)
        {
            /* already save them in reverse order here and permute,
             * values in reverse order, too*/
            coeffs[this->cols - 1 - col] = _mm256_permutevar8x32_epi32(
                prod,
                idx
            );
            long4_t elem = this->get(r1, col);
            elem = global::F.wide_mul(elem, prod);
            this->set(r1, col, elem);

            prod = global::F.wide_mul(prod, pac_gamma);
        }

        /* handle special permutations required in case not divisible by 4 */
        if (this->nmod)
        {
            switch (this->nmod)
            {
            case 1:
                COEFF_LOOP(1);
                break;
            case 2:
                COEFF_LOOP(2);
                break;
            case 3:
                COEFF_LOOP(3);
                break;
            case 4:
                COEFF_LOOP(4);
                break;
            case 5:
                COEFF_LOOP(5);
                break;
            case 6:
                COEFF_LOOP(6);
                break;
            case 7:
                COEFF_LOOP(7);
                break;
            }

            __m256i idx;
            /* lazy and ugly.. */
            switch (this->nmod)
            {
            case 1:
                idx = _mm256_set_epi32(
                    0b0,
                    0b0,
                    0b0,
                    0b0,
                    0b0,
                    0b0,
                    0b0,
                    0b0
                );
                break;
            case 2:
                idx = _mm256_set_epi32(
                    0b1,
                    0b0,
                    0b0,
                    0b0,
                    0b0,
                    0b0,
                    0b0,
                    0b0
                );
                break;
            case 3:
                idx = _mm256_set_epi32(
                    0b10,
                    0b1,
                    0b0,
                    0b0,
                    0b0,
                    0b0,
                    0b0,
                    0b0
                );
                break;
            case 4:
                idx = _mm256_set_epi32(
                    0b11,
                    0b10,
                    0b1,
                    0b0,
                    0b0,
                    0b0,
                    0b0,
                    0b0
                );
                break;
            case 5:
                idx = _mm256_set_epi32(
                    0b100,
                    0b11,
                    0b10,
                    0b1,
                    0b0,
                    0b0,
                    0b0,
                    0b0
                );
                break;
            case 6:
                idx = _mm256_set_epi32(
                    0b101,
                    0b100,
                    0b11,
                    0b10,
                    0b1,
                    0b0,
                    0b0,
                    0b0
                );
                break;
            case 7:
                idx = _mm256_set_epi32(
                    0b110,
                    0b101,
                    0b100,
                    0b11,
                    0b10,
                    0b1,
                    0b0,
                    0b0
                );
                break;
            }
            coeffs[this->cols - 1] = _mm256_permutevar8x32_epi32(
                coeffs[this->cols - 1],
                idx
            );

        }

        /* and do r2 left to right */
        for (int col = 0; col < this->cols; col++)
        {
            this->set(r2, col,
                      global::F.wide_mul(
                          this->get(r2, col),
                          coeffs[col]
                      )
                );
        }
    }

    void swap_rows(int r1, int r2)
    {
        long4_t tmp;
        for (int c = 0; c < this->cols; c++)
        {
            tmp = this->get(r1, c);
            this->set(r1, c, this->get(r2, c));
            this->set(r2, c, tmp);
        }
    }

    void mul_row(int row, uint64_t v)
    {
        long4_t pack = _mm256_set_epi64x(
            v << 32 | v,
            v << 32 | v,
            v << 32 | v,
            v << 32 | v
        );
        for (int col = 0; col < this->cols; col++)
            this->set(row, col,
                      global::F.wide_mul(this->get(row, col), pack)
                );
    }

    /* subtract v times r1 from r2 */
    void row_op(int r1, int r2, long4_t pack)
    {
        for (int col = 0; col < this->cols; col++)
        {
            long4_t tmp = global::F.wide_mul(this->get(r1, col), pack);

            this->set(r2, col,
                      _mm256_xor_si256(this->get(r2, col), tmp)
                );
        }
    }

    GF_element det()
    {
        uint64_t det = 0x1;
        for (int col = 0; col < this->cols; col++)
        {
            DET_LOOP(0);
            DET_LOOP(1);
            DET_LOOP(2);
            DET_LOOP(3);
            DET_LOOP(4);
            DET_LOOP(5);
            DET_LOOP(6);
            DET_LOOP(7);
        }
        return GF_element(det);
    }

    /* only used for testing */
    FMatrix unpack() const
    {
        std::valarray<GF_element> unpacked(this->rows * this->rows);

        for (int row = 0; row < this->rows; row++)
        {
            for (int col = 0; col < this->cols; col++)
            {
                for (int e = 0; e < VECTOR_N; e++)
                {
                    uint64_t rep;
                    switch (e/2)
                    {
                    case 0:
                        rep = _mm256_extract_epi64(this->get(row, col), 3);
                        break;
                    case 1:
                        rep = _mm256_extract_epi64(this->get(row, col), 2);
                        break;
                    case 2:
                        rep = _mm256_extract_epi64(this->get(row, col), 1);
                        break;
                    case 3:
                        rep = _mm256_extract_epi64(this->get(row, col), 0);
                        break;
                    }
                    if (e%2)
                        rep &= 0xFFFF;
                    else
                        rep >>= 32*(1 - e%2);
                    unpacked[row*this->rows + VECTOR_N*col + e] =
                        GF_element(rep);
                }
            }
        }

        uint64_t n = this->rows;
        uint64_t r = this->rows;
        if (this->nmod)
            n -= VECTOR_N - this->nmod;
        return FMatrix(n, unpacked[std::gslice(0, {n,n}, {r,1})]);
    }
};

#endif
