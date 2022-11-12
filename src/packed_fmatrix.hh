/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef P_FMATRIX_H
#define P_FMATRIX_H

#include <valarray>
#include <vector>
#include <immintrin.h>

#include "gf.hh"
#include "global.hh"
#include "fmatrix.hh"

typedef long long long4_t __attribute__ ((vector_size (32)));

#define VECTOR_N 16

#define MOVE_ID_VEC(index)                              \
    {                                                   \
        if (index%2)                                    \
            onepos = _mm256_slli_epi32(                 \
                _mm256_permutevar8x32_epi32(            \
                    onepos,                             \
                    permute_idx                         \
                ),                                      \
                16                                      \
            );                                          \
        else                                            \
            onepos = _mm256_srli_epi32(onepos, 16);     \
    }

#define COEFF_LOOP(index)                                               \
    {                                                                   \
        for (int col = 0; col < this->cols - 1; col++)                  \
            coeffs[col] = _mm256_blend_epi32(                           \
                _mm256_permutevar8x32_epi32(coeffs[col + 0], idx),      \
                _mm256_permutevar8x32_epi32(coeffs[col + 1], idx),      \
                0xFF >> index                                           \
            );                                                          \
    }

#define DET_LOOP(index)                                             \
    {                                                               \
        int r0 = VECTOR_N*col + index;                              \
        long4_t mx;                                                 \
        int mxi = -1;                                               \
        long4_t cmpmsk = _mm256_set_epi64x(                         \
            global::F.get_mask() << 48,                             \
            0,                                                      \
            0,                                                      \
            0                                                       \
        );                                                          \
        if (index >= (VECTOR_N / 2))                                \
            cmpmsk = _mm256_permute4x64_epi64(cmpmsk, 0x0C);        \
        cmpmsk = _mm256_srli_si256(cmpmsk, 2 * (index % 8));        \
        for (int row = r0; row < this->rows; row++)                 \
        {                                                           \
            char ZF = _mm256_testz_si256(                           \
                cmpmsk,                                             \
                this->get(row,col)                                  \
            );                                                      \
            if (ZF == 0)                                            \
            {                                                       \
                mx = this->get(row,col);                            \
                mxi = row;                                          \
                break;                                              \
            }                                                       \
        }                                                           \
        if (mxi == -1)                                              \
            return global::F.zero();                                \
        uint64_t mx_ext =                                           \
            _mm256_extract_epi16(mx, VECTOR_N - 1 - index);         \
        if (mxi != r0)                                              \
            this->swap_rows(mxi, r0);                               \
        /* vectorize? */                                            \
        det = global::F.rem(                                        \
            global::F.clmul(det, mx_ext)                            \
        );                                                          \
        mx_ext = global::F.ext_euclid(mx_ext);                      \
        this->mul_row(r0, mx_ext);                                  \
        /* vectorize end? */                                        \
        char mask = (VECTOR_N - 1 - index) / 2;                     \
        long4_t idx = _mm256_set1_epi32(mask);                      \
        for (int row = r0 + 1; row < this->rows; row++)             \
        {                                                           \
            long4_t val = _mm256_permutevar8x32_epi32(              \
                this->get(row, col),                                \
                idx                                                 \
            );                                                      \
            if (index % 2)                                          \
                val = _mm256_blend_epi16(                           \
                    val,                                            \
                    _mm256_slli_epi32(val, 16),                     \
                    0xAA                                            \
                );                                                  \
            else                                                    \
                val = _mm256_blend_epi16(                           \
                    _mm256_srli_epi32(val, 16),                     \
                    val,                                            \
                    0xAA                                            \
                );                                                  \
            this->row_op(r0, row, val);                             \
        }                                                           \
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
                          _mm256_set_epi16(
                              matrix(r, VECTOR_N*c +  0).get_repr(),
                              matrix(r, VECTOR_N*c +  1).get_repr(),
                              matrix(r, VECTOR_N*c +  2).get_repr(),
                              matrix(r, VECTOR_N*c +  3).get_repr(),
                              matrix(r, VECTOR_N*c +  4).get_repr(),
                              matrix(r, VECTOR_N*c +  5).get_repr(),
                              matrix(r, VECTOR_N*c +  6).get_repr(),
                              matrix(r, VECTOR_N*c +  7).get_repr(),
                              matrix(r, VECTOR_N*c +  8).get_repr(),
                              matrix(r, VECTOR_N*c +  9).get_repr(),
                              matrix(r, VECTOR_N*c + 10).get_repr(),
                              matrix(r, VECTOR_N*c + 11).get_repr(),
                              matrix(r, VECTOR_N*c + 12).get_repr(),
                              matrix(r, VECTOR_N*c + 13).get_repr(),
                              matrix(r, VECTOR_N*c + 14).get_repr(),
                              matrix(r, VECTOR_N*c + 15).get_repr()
                          )
                    );
            if (this->nmod)
            {
                int c = this->cols - 1;
                uint64_t elems[VECTOR_N];
                for (int i = 0; i < VECTOR_N; i++)
                    elems[i] = 0;
                for (int i = 0; i < this->nmod; i++)
                    elems[i] = matrix(r, VECTOR_N*c + i).get_repr();

                this->set(r, c,
                          _mm256_set_epi16(
                              elems[0],  elems[1],  elems[2],  elems[3],
                              elems[4],  elems[5],  elems[6],  elems[7],
                              elems[8],  elems[9],  elems[10], elems[11],
                              elems[12], elems[13], elems[14], elems[15]
                          )
                );
            }
        }
        /* cycle all one position to right */
        const long4_t permute_idx = _mm256_set_epi32(
            0b000,
            0b111,
            0b110,
            0b101,
            0b100,
            0b011,
            0b010,
            0b001
        );

        long4_t onepos = _mm256_set_epi64x(1ull << 48, 0, 0, 0);
        for (int i = 0; i < this->nmod; i++)
            MOVE_ID_VEC(i);

        for (int r = matrix.get_n(); r < this->rows; r++)
        {
            for (int c = 0; c < this->cols - 1; c++)
                this->set(r, c, _mm256_setzero_si256());
            this->set(r, this->cols - 1, onepos);
            MOVE_ID_VEC(r % VECTOR_N);
        }
    }

    void mul_gamma(int r1, int r2, const GF_element &gamma)
    {
        /* here we do r1 first left to right and save the auxiliary gamma vectors.
         * then we permute the gamma vectors as required (note this->nmod here),
         * and do r2 left to right */
        /* do gamma multiplication during initialization?? */
        long4_t pac_gamma = _mm256_set1_epi16(gamma.get_repr());
        for (int i = 0; i < util::log2(VECTOR_N); i++)
            pac_gamma = global::F.wide_mul(pac_gamma, pac_gamma);

        uint64_t elems[VECTOR_N];
        uint64_t g = 1ull;
        #pragma GCC unroll 32
        for (int i = 0; i < VECTOR_N; i++)
        {
            elems[i] = g;
            g = global::F.rem(global::F.clmul(g, gamma.get_repr()));
        }
        long4_t prod = _mm256_set_epi16(
            elems[0],  elems[1],  elems[2],  elems[3],
            elems[4],  elems[5],  elems[6],  elems[7],
            elems[8],  elems[9],  elems[10], elems[11],
            elems[12], elems[13], elems[14], elems[15]
        );
        std::vector<long4_t> coeffs(this->cols);
        long4_t idx = _mm256_set_epi32(
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
            long4_t tmp = _mm256_permutevar8x32_epi32(
                prod,
                idx
            );
            coeffs[this->cols - 1 - col] = _mm256_blend_epi16(
                _mm256_srli_epi32(tmp, 16),
                _mm256_slli_epi32(tmp, 16),
                0xAA
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
                idx = _mm256_set_epi32(
                    0b000,
                    0b111,
                    0b110,
                    0b101,
                    0b100,
                    0b011,
                    0b010,
                    0b001
                );
                COEFF_LOOP(1);
                break;
            case 2:
                idx = _mm256_set_epi32(
                    0b001,
                    0b000,
                    0b111,
                    0b110,
                    0b101,
                    0b100,
                    0b011,
                    0b010
                );
                COEFF_LOOP(2);
                break;
            case 3:
                idx = _mm256_set_epi32(
                    0b010,
                    0b001,
                    0b000,
                    0b111,
                    0b110,
                    0b101,
                    0b100,
                    0b011
                );
                COEFF_LOOP(3);
                break;
            case 4:
                idx = _mm256_set_epi32(
                    0b011,
                    0b010,
                    0b001,
                    0b000,
                    0b111,
                    0b110,
                    0b101,
                    0b100
                );
                COEFF_LOOP(4);
                break;
            case 5:
                idx = _mm256_set_epi32(
                    0b100,
                    0b011,
                    0b010,
                    0b001,
                    0b000,
                    0b111,
                    0b110,
                    0b101
                );
                COEFF_LOOP(5);
                break;
            case 6:
                idx = _mm256_set_epi32(
                    0b101,
                    0b100,
                    0b011,
                    0b010,
                    0b001,
                    0b000,
                    0b111,
                    0b110
                );
                COEFF_LOOP(6);
                break;
            case 7:
                idx = _mm256_set_epi32(
                    0b110,
                    0b101,
                    0b100,
                    0b011,
                    0b010,
                    0b001,
                    0b000,
                    0b111
                );
                COEFF_LOOP(7);
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
        long4_t pack = _mm256_set1_epi16(v);
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
            #pragma GCC unroll 32
            for (int loop = 0; loop < VECTOR_N; loop++)
                DET_LOOP(loop);
        }
        return GF_element(det);
    }
};

#endif
