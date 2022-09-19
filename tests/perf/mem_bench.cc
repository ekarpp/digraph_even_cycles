/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <omp.h>
#include <stdint.h>
#include <algorithm>

using namespace std;

#define CORES 24
// 4 GiBs
#define N (1ull << 29)
#define WARMUP (1ull << 25)
#define GIBS (N*8 / (1ull << 30))
// one line is 8*64 bits
#define LINE 8

int main(void)
{
    double start_t, end_t, delta;
    uint64_t *vec = (uint64_t*) aligned_alloc(64, sizeof(uint64_t)*N);
    uint64_t *wup = (uint64_t*) aligned_alloc(64, sizeof(uint64_t)*WARMUP);

    uint64_t val = time(nullptr);

    for (uint64_t i = 0; i < WARMUP; i++)
        wup[i] = val;

    start_t = omp_get_wtime();
    for (uint64_t i = 0; i < N; i++)
        vec[i] = val;
    end_t = omp_get_wtime();

    delta = end_t - start_t;
    cout << "wrote sequentially " << GIBS
         << " GiB in " << delta << " seconds, "
         << GIBS / delta << " GiB / s." << endl;

    uint64_t sum = 0;
    for (uint64_t i = 0; i < WARMUP; i++)
        sum += wup[i];

    sum = 0;
    start_t = omp_get_wtime();
    for (uint64_t i = 0; i < N; i++)
        sum += vec[i];
    end_t = omp_get_wtime();
    cout << sum << endl;

    delta = end_t - start_t;
    cout << "read sequentially " << GIBS
         << " GiB in " << delta << " seconds, "
         << GIBS / delta << " GiB / s." << endl;

    vector<uint64_t> idx(N);
    for (uint64_t i = 0; i < N; i++)
        idx[i] = i;

    random_shuffle(begin(idx), end(idx));

    start_t = omp_get_wtime();
    for (uint64_t i = 0; i < N; i++)
        sum += vec[idx[i]];
    end_t = omp_get_wtime();
    cout << sum << endl;

    delta = end_t - start_t;
    cout << "read sequentially random order " << GIBS
         << " GiB in " << delta << " seconds, "
         << GIBS / delta << " GiB / s." << endl;

    vector<uint64_t> line_idx(N / LINE);
    for (uint64_t i = 0; i < N / LINE; i++)
        line_idx[i] = i;

    random_shuffle(begin(line_idx), end(line_idx));

    start_t = omp_get_wtime();
    for (uint64_t i = 0; i < N / LINE; i++)
    {
        uint64_t l = line_idx[i];
        for (uint64_t j = 0; j < LINE; j++)
            sum += vec[l + j];
    }
    end_t = omp_get_wtime();
    cout << sum << endl;

    delta = end_t - start_t;
    cout << "read sequentially random cache lines " << GIBS
         << " GiB in " << delta << " seconds, "
         << GIBS / delta << " GiB / s." << endl;

    cout << "############################" << endl;

    val = time(nullptr);

    #pragma omp parallel for
    for (uint64_t i = 0; i < WARMUP; i++)
        wup[i] = val;

    start_t = omp_get_wtime();
    #pragma omp parallel for
    for (uint64_t i = 0; i < N; i++)
        vec[i] = val;
    end_t = omp_get_wtime();

    delta = end_t - start_t;
    cout << "wrote parallelly " << GIBS
         << " GiB in " << delta << " seconds, "
         << GIBS / delta << " GiB / s." << endl;

    uint64_t sums[CORES];
    int cores = CORES;
    uint64_t block_size = N / cores;
    start_t = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < cores; i++)
    {
        uint64_t csum = 0;
        uint64_t first = block_size * i;
        uint64_t last = block_size * (i + 1) - 1;
        for (uint64_t j = first; j < last; j++)
            csum += vec[j];
        sums[i] = csum;
    }
    end_t = omp_get_wtime();
    for (int i = 0; i < cores; i++)
        sum += sums[i];
    cout << sum << endl;

    delta = end_t - start_t;
    cout << "read parallelly " << GIBS
         << " GiB in " << delta << " seconds, "
         << GIBS / delta << " GiB / s." << endl;

    start_t = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < cores; i++)
    {
        uint64_t csum = 0;
        uint64_t first = block_size * i;
        uint64_t last = block_size * (i + 1) - 1;
        for (uint64_t j = first; j < last; j++)
            csum += vec[idx[j]];
        sums[i] = csum;
    }
    end_t = omp_get_wtime();
    for (int i = 0; i < cores; i++)
        sum += sums[i];
    cout << sum << endl;

    delta = end_t - start_t;
    cout << "read parallelly random order " << GIBS
         << " GiB in " << delta << " seconds, "
         << GIBS / delta << " GiB / s." << endl;

    start_t = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < cores; i++)
    {
        uint64_t csum = 0;
        uint64_t first = block_size * i / LINE;
        uint64_t last = block_size * (i + 1) / LINE - 1;
        for (uint64_t j = first; j < last; j++)
        {
            uint64_t l = line_idx[j];
            for (uint64_t k = 0; k < LINE; k++)
                csum += vec[l + k];
        }
        sums[i] = csum;
    }
    end_t = omp_get_wtime();
    for (int i = 0; i < cores; i++)
        sum += sums[i];
    cout << sum << endl;

    delta = end_t - start_t;
    cout << "read parallelly random cache lines " << GIBS
         << " GiB in " << delta << " seconds, "
         << GIBS / delta << " GiB / s." << endl;

    free(vec);
    free(wup);
    return 0;
}
