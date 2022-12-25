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

constexpr int VECTOR_N = 8;

class Packed_FMatrix
{
private:
    int rows;
    int cols;
    // original matrix n moduloe VECTOR_N
    int nmod;
    std::vector<long4_t> m;
    std::vector<long4_t> coeffs;

    const long4_t &get(const int row, const int col) const
    {
        return this->m[row*this->cols + col];
    }

    void set(const int row, const int col, const long4_t &v)
    {
        this->m[row*this->cols + col] = v;
    }

    template <int index>
    void coeff_loop()
    {
        /* cyclic permutation to the right.
         * maybe could use bit ops? (did not work past 128 bit lanes?) */
        const long4_t cycle_idx = _mm256_set_epi32(
            0b000,
            0b111,
            0b110,
            0b101,
            0b100,
            0b011,
            0b010,
            0b001
        );

        long4_t idx = _mm256_set_epi32(
            0b111,
            0b110,
            0b101,
            0b100,
            0b011,
            0b010,
            0b001,
            0b000
        );

        /* compiler should optimize this */
        for (int i = 0; i < index; i++)
            idx = _mm256_permutevar8x32_epi32(idx, cycle_idx);

        for (int col = 0; col < this->cols - 1; col++)
            this->coeffs[col] = _mm256_blend_epi32(
                _mm256_permutevar8x32_epi32(coeffs[col + 0], idx),
                _mm256_permutevar8x32_epi32(coeffs[col + 1], idx),
                0xFF >> index
            );

        this->coeffs[this->cols - 1] = _mm256_permutevar8x32_epi32(
            this->coeffs[this->cols - 1],
            idx
        );
    }

    /* returns true if zero det */
    template <int index>
    bool det_loop(const int col, uint64_t &det)
    {
        const int r0 = VECTOR_N*col + index;
        int piv_idx = -1;
        /* optimize */
        long4_t cmpmsk = _mm256_set_epi64x(
            0xFFFFull << 32,
            0,
            0,
            0
        );
        if (index >= 4)
            cmpmsk = _mm256_permute4x64_epi64(cmpmsk, 0x0C);
        cmpmsk = _mm256_srli_si256(cmpmsk, 4*(index%4));
        for (int row = r0; row < this->rows; row++)
        {
            const char ZF = _mm256_testz_si256(
                cmpmsk,
                this->get(row,col)
            );
            if (ZF == 0)
            {
                piv_idx = row;
                break;
            }
        }
        if (piv_idx == -1)
        {
            det = 0x0;
            return true;
        }
        if (piv_idx != r0)
            this->swap_rows(piv_idx, r0, col);

        constexpr char mask = VECTOR_N - 1 - index;
        /* rows got swapped */
        uint64_t pivot =
            _mm256_extract_epi32(this->get(r0, col), mask);
        /* vectorize? */
        det = global::F->rem(
            global::F->clmul(det, pivot)
        );
        pivot = global::F->ext_euclid(pivot);
        this->mul_row(r0, col, _mm256_set1_epi32(pivot));

        /* vectorize end? */
        const long4_t idx = _mm256_set1_epi32(mask);

        for (int row = r0 + 1; row < this->rows; row++)
        {
            const long4_t val = _mm256_permutevar8x32_epi32(
                this->get(row, col),
                idx
            );
            this->row_op(r0, row, col, val);
        }

        return false;
    }

    /* starting from column idx */
    inline void swap_rows(const int r1, const int r2, const int idx)
    {
        for (int col = idx; col < this->cols; col++)
        {
            const long4_t tmp = this->get(r1, col);
            this->set(r1, col, this->get(r2, col));
            this->set(r2, col, tmp);
        }
    }

    /* starting from column idx */
    inline void mul_row(const int row, const int idx, const long4_t &pack)
    {
        for (int col = idx; col < this->cols; col++)
            this->set(row, col,
                      global::F->wide_mul(this->get(row, col), pack)
                );
    }

    /* subtract v times r1 from r2, starting from column idx */
    inline void row_op(const int r1,
                       const int r2,
                       const int idx,
                       const long4_t &pack
    )
    {
        for (int col = idx; col < this->cols; col++)
        {
            const long4_t tmp = global::F->wide_mul(this->get(r1, col), pack);

            this->set(r2, col,
                      _mm256_xor_si256(this->get(r2, col), tmp)
            );
        }
    }

public:
    explicit Packed_FMatrix(const int n)
    {
        this->nmod = n % VECTOR_N;
        this->rows = n;
        if (this->rows % VECTOR_N)
            this->rows += VECTOR_N - (n % VECTOR_N);
        this->cols = this->rows / VECTOR_N;

        this->m.resize(this->rows * this->cols);
        this->coeffs.resize(this->cols);
    }

    void init(const FMatrix &matrix)
    {
        for (int r = 0; r < matrix.get_n(); r++)
        {
            for (int c = 0; c < matrix.get_n() / VECTOR_N; c++)
                this->set(r, c,
                          _mm256_set_epi32(
                              matrix(r, VECTOR_N*c + 0).get_repr(),
                              matrix(r, VECTOR_N*c + 1).get_repr(),
                              matrix(r, VECTOR_N*c + 2).get_repr(),
                              matrix(r, VECTOR_N*c + 3).get_repr(),
                              matrix(r, VECTOR_N*c + 4).get_repr(),
                              matrix(r, VECTOR_N*c + 5).get_repr(),
                              matrix(r, VECTOR_N*c + 6).get_repr(),
                              matrix(r, VECTOR_N*c + 7).get_repr()
                              )
                    );
            if (this->nmod)
            {
                const int c = this->cols - 1;
                uint64_t elems[VECTOR_N];
                for (int i = 0; i < VECTOR_N; i++)
                    elems[i] = 0;
                for (int i = 0; i < this->nmod; i++)
                    elems[i] = matrix(r, VECTOR_N*c + i).get_repr();

                this->set(r, c, _mm256_set_epi32(
                              elems[0], elems[1],
                              elems[2], elems[3],
                              elems[4], elems[5],
                              elems[6], elems[7]
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

    void mul_gamma(const int r1, const int r2, const GF_element &gamma)
    {
        /* here we do r1 first left to right and save the auxiliary gamma vectors.
         * then we permute the gamma vectors as required (note this->nmod here),
         * and do r2 left to right */
        /* do gamma multiplication during initialization?? */
        long4_t pac_gamma = _mm256_set1_epi32(gamma.get_repr());
        // pac_gamma = [gamma^VECTOR_N]
        pac_gamma = global::F->wide_mul(pac_gamma, pac_gamma);
        pac_gamma = global::F->wide_mul(pac_gamma, pac_gamma);
        pac_gamma = global::F->wide_mul(pac_gamma, pac_gamma);

        uint64_t elems[VECTOR_N];
        uint64_t g = 1ull;
        for (int i = 0; i < VECTOR_N; i++)
        {
            elems[i] = g;
            g = global::F->rem(global::F->clmul(g, gamma.get_repr()));
        }
        long4_t prod = _mm256_set_epi32(
            elems[0], elems[1],
            elems[2], elems[3],
            elems[4], elems[5],
            elems[6], elems[7]
        );

        /* reverse order permutation of 32 bits */
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
            this->coeffs[this->cols - 1 - col] = _mm256_permutevar8x32_epi32(
                prod,
                idx
            );
            const long4_t elem = global::F->wide_mul(this->get(r1, col), prod);
            this->set(r1, col, elem);

            prod = global::F->wide_mul(prod, pac_gamma);
        }

        /* handle special permutations required in case not divisible by 4 */
        if (this->nmod)
        {
            switch (this->nmod)
            {
            case 1:
                coeff_loop<1>();
                break;
            case 2:
                coeff_loop<2>();
                break;
            case 3:
                coeff_loop<3>();
                break;
            case 4:
                coeff_loop<4>();
                break;
            case 5:
                coeff_loop<5>();
                break;
            case 6:
                coeff_loop<6>();
                break;
            case 7:
                coeff_loop<7>();
                break;
            }
        }

        /* and do r2 left to right */
        for (int col = 0; col < this->cols; col++)
        {
            this->set(r2, col,
                      global::F->wide_mul(
                          this->get(r2, col),
                          this->coeffs[col]
                      )
                );
        }
    }

    GF_element det()
    {
        uint64_t det = 0x1;
        for (int col = 0; col < this->cols; col++)
        {
            /* each "column" is a vector that has 8 real columns.
             * each call returns true if the determinant is zero */
            if (det_loop<0>(col, det) || det_loop<1>(col, det)
                || det_loop<2>(col, det) || det_loop<3>(col, det)
                || det_loop<4>(col, det) || det_loop<5>(col, det)
                || det_loop<6>(col, det) || det_loop<7>(col, det))
                break;
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

        const uint64_t r = this->rows;
        uint64_t n = this->rows;
        if (this->nmod)
            n -= VECTOR_N - this->nmod;
        return FMatrix(n, unpacked[std::gslice(0, {n,n}, {r,1})]);
    }
};

#endif
