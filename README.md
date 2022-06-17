# The esa-matchfinder

The esa-matchfinder is a C99 library for efficient Lempel-Ziv factorization using enhanced suffix array (ESA).

Copyright (c) 2022 Ilya Grebnov <ilya.grebnov@gmail.com>

> * The esa-matchfinder is block based algorithm with maximum supported block size of 512 megabytes finding matches in range of 2..64 bytes using 12x bytes of extra memory. ESA_MATCHFINDER_MATCH_BITS definition could be changed to support larger match finding range, but with reduction in maximum supported block size.
> * The esa-matchfinder does not employ any heuristics or search depth limitations and always finds distance optimal matches even on highly repetitive sources. The only exception is matches at beginning at the block; due to implementation details the esa-matchfinder can not find any matches with offset 0.
> * The esa-matchfinder is fast in best, average and worst cases (see [Benchmarks](#benchmarks) below). But the esa-matchfinder is sensitive to fast memory and software prefetching and might not be suitable for some CPU architectures. The esa-matchfinder might also be suboptimal on specific data types (some DNA sequences in particular), so please benchmark yourself.
> * The esa-matchfinder works with compilers from Microsoft and GNU, but I recommend Clang for best performance. Additionally, the esa-matchfinder is designed for 64-bit systems and will work suboptimally on 32-bit system.

## Algorithm 
> The esa-matchfinder uses methodology of bottom-up traversal of the Longest Common Prefix (LCP) interval tree written in 2014-2015 by Eric Biggers <ebiggers3@gmail.com> and dedicated to the public domain worldwide.

The esa-matchfinder finds all distance optimal matches (between min_match_length and max_match_length inclusive) for every position of the input block using following algorithm:

1. Suffix (SA) and longest common prefix (LCP) arrays are constructed for the input block. Next, interval tree is constructed on top of SA and LCP arrays.
> * The data structure consisting of SA and LCP is often referred as enhanced suffix array (ESA). Hence the name of the match finder.

2. Each interval is a maximum range (could not be further extended to the left or right) of suffixes in SA with common prefix of certain length (LCP). This intervals represent internal nodes of suffix tree.
3. Using interval tree we can now traverse up or down to either wider interval with smaller common prefix or narrower intervals with larger common prefixes.
4. For purpose of Lempel-Ziv factorization we only need to support bottom-up traversal, so during interval tree construction we only need to capture link to parent interval and length of interval's common prefix.
5. LCP array is also pruned by min_match_length and max_match_length to reduce size and depth of interval tree.
6. Additionally, for each position of input block we capture link to a leaf interval corresponding to that position, so we can start bottom-up traversal during factorization phase.
7. SA, LCP and interval tree construction is done during input block parsing phase in linear time with optional multi-threaded optimization using OpenMP.
8. LZ factorization phase is done from left to right by bottom-up traversal of interval tree for each position from input block by reading and updating each corresponding interval with latest offset.

## License
The esa-matchfinder released under the [Apache License Version 2.0](LICENSE "Apache license") and is considered suitable for production use. However, no warranty or fitness for a particular purpose is expressed or implied.

## Changes
* June 12, 2022 (1.0.0)
  * Initial public release of the esa-matchfinder.

## Example of usage (See [esa_matchfinder.h](esa_matchfinder.h) for complete APIs list)
```c
#include "esa_matchfinder.h"

long long multi_pass_optimal_parse(const unsigned char * buffer, int size)
{
    long long total_matches = 0;

    void * mf = esa_matchfinder_create(size, /*min_match_length*/ 2, /*max_match_length*/ 64);
    if (mf != NULL && esa_matchfinder_parse(mf, buffer, size) == ESA_MATCHFINDER_NO_ERROR)
    {
        for (int pass = 0; pass < 2; pass += 1)
        {
            ESA_MATCHFINDER_MATCH matches[ESA_MATCHFINDER_MAX_MATCH_LENGTH];

            esa_matchfinder_rewind(mf, /*position*/ 0);

            for (int position = 0; position < size; position += 1)
            {
                total_matches += esa_matchfinder_find_all_matches(mf, matches) - matches;
            }
        }
    }

    esa_matchfinder_destroy(mf);

    return total_matches;
}
```

---

# Benchmarks #
  
## Methodology ##
  * Input files were capped at 510MB for tests on x86-64 architecture and 128MB for tests on ARMv8 architecture
  * For all match finders maximum match length was set to 64 bytes (other parameters were not changed)
  * Optimality is defined as percentage of optimal matches (longest possible length) found across all possible matches
  * The timings are minimum of five runs measuring single-threaded performance in *optimal* parsing mode

## Specification (x86-64 architecture) ##
  * OS: Microsoft Windows 10 Pro 64-Bit
  * CPU: Intel Core i7-9700K Processor (12M Cache, 5GHz)
  * RAM: 2x8 GB dual-channel DDR4 (4133 MHz, 17-17-17-37)
  * Compiler: Microsoft Visual C++ compiler v14.32 
  * Optimizations: /MD /DNDEBUG /O2 /GL /arch:AVX2

### Silesia Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | optimality | LZMA HC4 v21.07 | optimality | LZMA BT4 v21.07 | optimality |
|:----:|:----:|:----------------------:|:----------:|:---------------:|:----------:|:---------------:|:----------:|
| dickens | 10192446 | **0.715 sec (14.26 MB/s)** | 100.00% | 3.632 sec (2.81 MB/s) | 40.43% | 2.880 sec (3.54 MB/s) | 99.84% |
| mozilla | 51220480 | **3.533 sec (14.50 MB/s)** | 100.00% | 12.364 sec (4.14 MB/s) | 55.98% | 8.546 sec (5.99 MB/s) | 80.08% |
| mr | 9970564 | **0.865 sec (11.53 MB/s)** | 100.00% | 1.680 sec (5.93 MB/s) | 56.13% | 1.688 sec (5.91 MB/s) | 94.35% |
| nci | 33553445 | **3.626 sec (9.25 MB/s)** | 100.00% | 6.155 sec (5.45 MB/s) | 16.25% | 5.957 sec (5.63 MB/s) | 61.62% |
| ooffice | 6152192 | **0.324 sec (18.99 MB/s)** | 100.00% | 1.120 sec (5.49 MB/s) | 70.35% | 0.780 sec (7.89 MB/s) | 84.30% |
| osdb | 10085684 | **0.599 sec (16.84 MB/s)** | 100.00% | 2.191 sec (4.60 MB/s) | 35.58% | 1.716 sec (5.88 MB/s) | 85.83% |
| reymont | 6627202 | **0.478 sec (13.86 MB/s)** | 100.00% | 1.644 sec (4.03 MB/s) | 27.64% | 1.586 sec (4.18 MB/s) | 99.43% |
| samba | 21606400 | **1.484 sec (14.56 MB/s)** | 100.00% | 4.093 sec (5.28 MB/s) | 63.92% | 2.932 sec (7.37 MB/s) | 85.49% |
| sao | 7251944 | **0.405 sec (17.91 MB/s)** | 100.00% | 1.345 sec (5.39 MB/s) | 34.99% | 0.997 sec (7.27 MB/s) | 51.66% |
| webster | 41458703 | **4.388 sec (9.45 MB/s)** | 100.00% | 15.712 sec (2.64 MB/s) | 38.14% | 12.630 sec (3.28 MB/s) | 95.90% |
| x-ray | 8474240 | **0.412 sec (20.57 MB/s)** | 100.00% | 1.751 sec (4.84 MB/s) | 90.99% | 1.081 sec (7.84 MB/s) | 95.25% |
| xml | 5345280 | **0.278 sec (19.23 MB/s)** | 100.00% | 0.777 sec (6.88 MB/s) | 57.68% | 0.667 sec (8.01 MB/s) | 98.50% |

### Large Canterbury Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | optimality | LZMA HC4 v21.07 | optimality | LZMA BT4 v21.07 | optimality |
|:----:|:----:|:----------------------:|:----------:|:---------------:|:----------:|:---------------:|:----------:|
| bible.txt | 4047392 | **0.256 sec (15.81 MB/s)** | 100.00% | 1.021 sec (3.96 MB/s) | 50.12% | 0.806 sec (5.02 MB/s) | 99.83% |
| E.coli | 4638690 | **0.327 sec (14.19 MB/s)** | 100.00% | 0.817 sec (5.68 MB/s) | 2.70% | 1.603 sec (2.89 MB/s) | 99.95% |
| world192.txt | 2473400 | **0.138 sec (17.92 MB/s)** | 100.00% | 0.525 sec (4.71 MB/s) | 63.90% | 0.363 sec (6.81 MB/s) | 99.10% |

### Manzini Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | optimality | LZMA HC4 v21.07 | optimality | LZMA BT4 v21.07 | optimality |
|:----:|:----:|:----------------------:|:----------:|:---------------:|:----------:|:---------------:|:----------:|
| chr22.dna | 34553758 | **4.651 sec (7.43 MB/s)** | 100.00% | 6.019 sec (5.74 MB/s) | 4.67% | 18.142 sec (1.90 MB/s) | 96.22% |
| etext99 | 105277340 | **12.065 sec (8.73 MB/s)** | 100.00% | 46.457 sec (2.27 MB/s) | 24.64% | 46.607 sec (2.26 MB/s) | 98.62% |
| gcc-3.0.tar | 86630400 | **8.230 sec (10.53 MB/s)** | 100.00% | 23.390 sec (3.70 MB/s) | 59.21% | 17.811 sec (4.86 MB/s) | 94.33% |
| howto | 39422105 | **3.259 sec (12.10 MB/s)** | 100.00% | 14.254 sec (2.77 MB/s) | 45.33% | 10.613 sec (3.71 MB/s) | 95.05% |
| jdk13c | 69728899 | **6.522 sec (10.69 MB/s)** | 100.00% | 11.797 sec (5.91 MB/s) | 61.10% | 10.246 sec (6.81 MB/s) | 92.01% |
| linux-2.4.5.tar | 116254720 | **11.075 sec (10.50 MB/s)** | 100.00% | 36.868 sec (3.15 MB/s) | 57.33% | 27.014 sec (4.30 MB/s) | 93.13% |
| rctail96 | 114711151 | **10.375 sec (11.06 MB/s)** | 100.00% | 31.860 sec (3.60 MB/s) | 54.54% | 26.365 sec (4.35 MB/s) | 98.46% |
| rfc | 116421901 | **14.677 sec (7.93 MB/s)** | 100.00% | 37.246 sec (3.13 MB/s) | 35.41% | 33.363 sec (3.49 MB/s) | 87.53% |
| sprot34.dat | 109617186 | **11.628 sec (9.43 MB/s)** | 100.00% | 53.546 sec (2.05 MB/s) | 59.60% | 28.559 sec (3.84 MB/s) | 93.73% |
| w3c2 | 104201579 | **9.480 sec (10.99 MB/s)** | 100.00% | 19.914 sec (5.23 MB/s) | 64.99% | 15.080 sec (6.91 MB/s) | 93.09% |

### Large Text Compression Benchmark Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | optimality | LZMA HC4 v21.07 | optimality | LZMA BT4 v21.07 | optimality |
|:----:|:----:|:----------------------:|:----------:|:---------------:|:----------:|:---------------:|:----------:|
| enwik8 | 100000000 | **10.241 sec (9.76 MB/s)** | 100.00% | 48.950 sec (2.04 MB/s) | 33.41% | 40.285 sec (2.48 MB/s) | 96.56% |
| enwik9 | 534773760 | **81.713 sec (6.54 MB/s)** | 100.00% | 277.775 sec (1.93 MB/s) | 24.76% | 257.795 sec (2.07 MB/s) | 91.65% |

### The Gauntlet Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | optimality | LZMA HC4 v21.07 | optimality | LZMA BT4 v21.07 | optimality |
|:----:|:----:|:----------------------:|:----------:|:---------------:|:----------:|:---------------:|:----------:|
| abac | 200000 | 0.013 sec (15.38 MB/s) | 100.00% | **0.006 sec (33.33 MB/s)** | 100.00% | 0.008 sec (25.00 MB/s) | 100.00% |
| abba | 10500596 | **1.045 sec (10.05 MB/s)** | 100.00% | 2.926 sec (3.59 MB/s) | 0.03% | 2.789 sec (3.77 MB/s) | 97.64% |
| book1x20 | 15375420 | **0.952 sec (16.15 MB/s)** | 100.00% | 4.008 sec (3.84 MB/s) | 30.82% | 2.912 sec (5.28 MB/s) | 99.87% |
| fib_s14930352 | 14930352 | **0.861 sec (17.34 MB/s)** | 100.00% | 1.844 sec (8.10 MB/s) | 96.99% | 1.057 sec (14.13 MB/s) | 100.00% |
| fss10 | 12078908 | **0.660 sec (18.30 MB/s)** | 100.00% | 1.682 sec (7.18 MB/s) | 95.65% | 0.989 sec (12.21 MB/s) | 100.00% |
| fss9 | 2851443 | **0.128 sec (22.28 MB/s)** | 100.00% | 0.397 sec (7.18 MB/s) | 95.65% | 0.233 sec (12.24 MB/s) | 100.00% |
| houston | 3839141 | 0.257 sec (14.94 MB/s) | 100.00% | **0.141 sec (27.23 MB/s)** | 97.64% | 0.173 sec (22.19 MB/s) | 100.00% |
| paper5x80 | 956322 | **0.036 sec (26.56 MB/s)** | 100.00% | 0.073 sec (13.10 MB/s) | 94.59% | 0.071 sec (13.47 MB/s) | 99.94% |
| test1 | 2097152 | 0.080 sec (26.21 MB/s) | 100.00% | **0.070 sec (29.96 MB/s)** | 100.00% | 0.101 sec (20.76 MB/s) | 100.00% |
| test2 | 2097152 | 0.080 sec (26.21 MB/s) | 100.00% | **0.070 sec (29.96 MB/s)** | 100.00% | 0.101 sec (20.76 MB/s) | 100.00% |
| test3 | 2097088 | 0.077 sec (27.23 MB/s) | 100.00% | **0.071 sec (29.54 MB/s)** | 100.00% | 0.086 sec (24.38 MB/s) | 100.00% |

### Pizza & Chilli Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | optimality | LZMA HC4 v21.07 | optimality | LZMA BT4 v21.07 | optimality |
|:----:|:----:|:----------------------:|:----------:|:---------------:|:----------:|:---------------:|:----------:|
| dblp.xml | 296135874 | **39.566 sec (7.48 MB/s)** | 100.00% | 100.764 sec (2.94 MB/s) | 46.52% | 74.943 sec (3.95 MB/s) | 93.38% |
| dna | 403927746 | 98.964 sec (4.08 MB/s) | 100.00% | **72.058 sec (5.61 MB/s)** | 0.51% | 314.605 sec (1.28 MB/s) | 63.97% |
| english.1024MB | 534773760 | **94.228 sec (5.68 MB/s)** | 100.00% | 266.114 sec (2.01 MB/s) | 15.29% | 309.063 sec (1.73 MB/s) | 95.59% |
| pitches | 55832855 | **3.698 sec (15.10 MB/s)** | 100.00% | 16.901 sec (3.30 MB/s) | 79.95% | 10.341 sec (5.40 MB/s) | 95.11% |
| proteins | 534773760 | **69.993 sec (7.64 MB/s)** | 100.00% | 463.257 sec (1.15 MB/s) | 43.66% | 320.663 sec (1.67 MB/s) | 99.85% |
| sources | 210866607 | **22.116 sec (9.53 MB/s)** | 100.00% | 73.429 sec (2.87 MB/s) | 54.69% | 56.835 sec (3.71 MB/s) | 93.25% |

### Pizza & Chilli Repetitive Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | optimality | LZMA HC4 v21.07 | optimality | LZMA BT4 v21.07 | optimality |
|:----:|:----:|:----------------------:|:----------:|:---------------:|:----------:|:---------------:|:----------:|
| cere | 461286644 | **68.413 sec (6.74 MB/s)** | 100.00% | 77.037 sec (5.99 MB/s) | 7.30% | 266.363 sec (1.73 MB/s) | 94.04% |
| coreutils | 205281778 | **18.967 sec (10.82 MB/s)** | 100.00% | 54.280 sec (3.78 MB/s) | 34.14% | 40.171 sec (5.11 MB/s) | 87.59% |
| einstein.de.txt | 92758441 | **6.620 sec (14.01 MB/s)** | 100.00% | 10.236 sec (9.06 MB/s) | 80.28% | 9.088 sec (10.21 MB/s) | 96.82% |
| einstein.en.txt | 467626544 | **39.117 sec (11.95 MB/s)** | 100.00% | 56.046 sec (8.34 MB/s) | 77.54% | 48.619 sec (9.62 MB/s) | 97.64% |
| Escherichia_Coli | 112689515 | **16.191 sec (6.96 MB/s)** | 100.00% | 19.889 sec (5.67 MB/s) | 0.26% | 63.204 sec (1.78 MB/s) | 96.93% |
| influenza | 154808555 | 44.042 sec (3.52 MB/s) | 100.00% | **25.125 sec (6.16 MB/s)** | 28.19% | 27.703 sec (5.59 MB/s) | 98.49% |
| kernel | 257961616 | **21.931 sec (11.76 MB/s)** | 100.00% | 81.326 sec (3.17 MB/s) | 24.72% | 59.968 sec (4.30 MB/s) | 93.18% |
| para | 429265758 | **69.888 sec (6.14 MB/s)** | 100.00% | 73.130 sec (5.87 MB/s) | 3.94% | 251.836 sec (1.70 MB/s) | 93.09% |
| world_leaders | 46968181 | **5.758 sec (8.16 MB/s)** | 100.00% | 6.124 sec (7.67 MB/s) | 40.96% | 6.126 sec (7.67 MB/s) | 62.59% |
| dblp.xml.00001.1 | 104857600 | **12.065 sec (8.69 MB/s)** | 100.00% | 18.851 sec (5.56 MB/s) | 32.29% | 16.655 sec (6.30 MB/s) | 96.79% |
| dblp.xml.00001.2 | 104857600 | **12.765 sec (8.21 MB/s)** | 100.00% | 18.972 sec (5.53 MB/s) | 32.01% | 16.767 sec (6.25 MB/s) | 96.99% |
| dblp.xml.0001.1 | 104857600 | **13.064 sec (8.03 MB/s)** | 100.00% | 18.952 sec (5.53 MB/s) | 32.36% | 16.712 sec (6.27 MB/s) | 96.76% |
| dblp.xml.0001.2 | 104857600 | **15.651 sec (6.70 MB/s)** | 100.00% | 19.258 sec (5.44 MB/s) | 30.31% | 17.409 sec (6.02 MB/s) | 96.99% |
| dna.001.1 | 104857600 | **16.519 sec (6.35 MB/s)** | 100.00% | 18.601 sec (5.64 MB/s) | 0.23% | 33.080 sec (3.17 MB/s) | 99.49% |
| english.001.2 | 104857600 | **11.087 sec (9.46 MB/s)** | 100.00% | 29.621 sec (3.54 MB/s) | 36.64% | 21.313 sec (4.92 MB/s) | 99.32% |
| proteins.001.1 | 104857600 | 12.562 sec (8.35 MB/s) | 100.00% | 24.271 sec (4.32 MB/s) | 89.70% | **12.411 sec (8.45 MB/s)** | 99.99% |
| sources.001.2 | 104857600 | **10.484 sec (10.00 MB/s)** | 100.00% | 20.817 sec (5.04 MB/s) | 43.87% | 16.775 sec (6.25 MB/s) | 97.24% |
| fib41 | 267914296 | 19.852 sec (13.50 MB/s) | 100.00% | 33.175 sec (8.08 MB/s) | 96.99% | **19.070 sec (14.05 MB/s)** | 100.00% |
| rs.13 | 216747218 | **16.024 sec (13.53 MB/s)** | 100.00% | 30.276 sec (7.16 MB/s) | 95.65% | 17.824 sec (12.16 MB/s) | 100.00% |
| tm29 | 268435456 | **19.926 sec (13.47 MB/s)** | 100.00% | 51.541 sec (5.21 MB/s) | 89.58% | 26.480 sec (10.14 MB/s) | 100.00% |

## Specification (ARMv8 architecture) ##
  * OS: Ubuntu 20.04 LTS 64-Bit
  * CPU: ODROID-N2+ Amlogic S922X (Cortex-A73 2.4Ghz)
  * RAM: 4GB LPDDR4 (2666 MHz) 
  * Compiler: Clang v10.0.0
  * Optimizations: -DNDEBUG -O3 -flto=thin -mcpu=native

### Silesia Corpus (ARMv8 architecture) ###

| file | size | esa-matchfinder v1.0.0 | optimality | LZMA HC4 v21.07 | optimality | LZMA BT4 v21.07 | optimality |
|:----:|:----:|:----------------------:|:----------:|:---------------:|:----------:|:---------------:|:----------:|
| dickens | 10192446 | **4.595 sec (2.22 MB/s)** | 100.00% | 18.955 sec (0.54 MB/s) | 40.43% | 9.678 sec (1.05 MB/s) | 99.84% |
| mozilla | 51220480 | **18.435 sec (2.78 MB/s)** | 100.00% | 43.410 sec (1.18 MB/s) | 55.98% | 24.890 sec (2.06 MB/s) | 80.08% |
| mr | 9970564 | **4.070 sec (2.45 MB/s)** | 100.00% | 7.299 sec (1.37 MB/s) | 56.13% | 5.935 sec (1.68 MB/s) | 94.35% |
| nci | 33553445 | 18.982 sec (1.77 MB/s) | 100.00% | 18.183 sec (1.85 MB/s) | 16.25% | **16.916 sec (1.98 MB/s)** | 61.62% |
| ooffice | 6152192 | **1.884 sec (3.27 MB/s)** | 100.00% | 5.239 sec (1.17 MB/s) | 70.35% | 2.779 sec (2.21 MB/s) | 84.30% |
| osdb | 10085684 | **4.090 sec (2.47 MB/s)** | 100.00% | 9.079 sec (1.11 MB/s) | 35.58% | 5.522 sec (1.83 MB/s) | 85.83% |
| reymont | 6627202 | **3.094 sec (2.14 MB/s)** | 100.00% | 6.582 sec (1.01 MB/s) | 27.64% | 5.023 sec (1.32 MB/s) | 99.43% |
| samba | 21606400 | **8.063 sec (2.68 MB/s)** | 100.00% | 16.016 sec (1.35 MB/s) | 63.92% | 8.947 sec (2.42 MB/s) | 85.49% |
| sao | 7251944 | **2.553 sec (2.84 MB/s)** | 100.00% | 6.461 sec (1.12 MB/s) | 34.99% | 3.753 sec (1.93 MB/s) | 51.66% |
| webster | 41458703 | **23.900 sec (1.73 MB/s)** | 100.00% | 67.306 sec (0.62 MB/s) | 38.14% | 36.756 sec (1.13 MB/s) | 95.90% |
| x-ray | 8474240 | **2.221 sec (3.82 MB/s)** | 100.00% | 7.558 sec (1.12 MB/s) | 90.99% | 3.921 sec (2.16 MB/s) | 95.25% |
| xml | 5345280 | **1.590 sec (3.36 MB/s)** | 100.00% | 2.380 sec (2.25 MB/s) | 57.68% | 1.910 sec (2.80 MB/s) | 98.50% |

### Large Canterbury Corpus (ARMv8 architecture) ###

| file | size | esa-matchfinder v1.0.0 | optimality | LZMA HC4 v21.07 | optimality | LZMA BT4 v21.07 | optimality |
|:----:|:----:|:----------------------:|:----------:|:---------------:|:----------:|:---------------:|:----------:|
| bible.txt | 4047392 | **1.619 sec (2.50 MB/s)** | 100.00% | 4.741 sec (0.85 MB/s) | 50.12% | 2.688 sec (1.51 MB/s) | 99.83% |
| E.coli | 4638690 | 2.068 sec (2.24 MB/s) | 100.00% | **1.988 sec (2.33 MB/s)** | 2.70% | 5.537 sec (0.84 MB/s) | 99.95% |
| world192.txt | 2473400 | **0.914 sec (2.71 MB/s)** | 100.00% | 2.583 sec (0.96 MB/s) | 63.90% | 1.243 sec (1.99 MB/s) | 99.10% |

### Manzini Corpus (ARMv8 architecture) ###

| file | size | esa-matchfinder v1.0.0 | optimality | LZMA HC4 v21.07 | optimality | LZMA BT4 v21.07 | optimality |
|:----:|:----:|:----------------------:|:----------:|:---------------:|:----------:|:---------------:|:----------:|
| chr22.dna | 34553758 | 24.388 sec (1.42 MB/s) | 100.00% | **14.290 sec (2.42 MB/s)** | 4.67% | 54.272 sec (0.64 MB/s) | 96.22% |
| etext99 | 105277340 | **66.408 sec (1.59 MB/s)** | 100.00% | 216.948 sec (0.49 MB/s) | 24.64% | 133.282 sec (0.79 MB/s) | 98.62% |
| gcc-3.0.tar | 86630400 | **41.198 sec (2.10 MB/s)** | 100.00% | 89.721 sec (0.97 MB/s) | 59.21% | 51.345 sec (1.69 MB/s) | 94.33% |
| howto | 39422105 | **18.222 sec (2.16 MB/s)** | 100.00% | 61.739 sec (0.64 MB/s) | 45.33% | 31.992 sec (1.23 MB/s) | 95.05% |
| jdk13c | 69728899 | 32.194 sec (2.17 MB/s) | 100.00% | 39.770 sec (1.75 MB/s) | 61.10% | **30.121 sec (2.31 MB/s)** | 92.01% |
| linux-2.4.5.tar | 116254720 | **55.748 sec (2.09 MB/s)** | 100.00% | 135.994 sec (0.85 MB/s) | 57.33% | 75.110 sec (1.55 MB/s) | 93.13% |
| rctail96 | 114711151 | **55.181 sec (2.08 MB/s)** | 100.00% | 133.541 sec (0.86 MB/s) | 54.54% | 80.031 sec (1.43 MB/s) | 98.46% |
| rfc | 116421901 | **73.199 sec (1.59 MB/s)** | 100.00% | 149.776 sec (0.78 MB/s) | 35.41% | 96.395 sec (1.21 MB/s) | 87.53% |
| sprot34.dat | 109617186 | **58.204 sec (1.88 MB/s)** | 100.00% | 179.813 sec (0.61 MB/s) | 59.60% | 80.105 sec (1.37 MB/s) | 93.73% |
| w3c2 | 104201579 | 49.571 sec (2.10 MB/s) | 100.00% | 70.244 sec (1.48 MB/s) | 64.99% | **45.104 sec (2.31 MB/s)** | 93.09% |

### Large Text Compression Benchmark Corpus (ARMv8 architecture) ###

| file | size | esa-matchfinder v1.0.0 | optimality | LZMA HC4 v21.07 | optimality | LZMA BT4 v21.07 | optimality |
|:----:|:----:|:----------------------:|:----------:|:---------------:|:----------:|:---------------:|:----------:|
| enwik8 | 100000000 | **53.981 sec (1.85 MB/s)** | 100.00% | 206.851 sec (0.48 MB/s) | 33.41% | 112.780 sec (0.89 MB/s) | 96.56% |

### The Gauntlet Corpus (ARMv8 architecture) ###

| file | size | esa-matchfinder v1.0.0 | optimality | LZMA HC4 v21.07 | optimality | LZMA BT4 v21.07 | optimality |
|:----:|:----:|:----------------------:|:----------:|:---------------:|:----------:|:---------------:|:----------:|
| abac | 200000 | 0.045 sec (4.40 MB/s) | 100.00% | **0.023 sec (8.78 MB/s)** | 100.00% | 0.046 sec (4.31 MB/s) | 100.00% |
| abba | 10500596 | 7.005 sec (1.50 MB/s) | 100.00% | **6.204 sec (1.69 MB/s)** | 0.03% | 7.848 sec (1.34 MB/s) | 97.64% |
| book1x20 | 15375420 | **6.223 sec (2.47 MB/s)** | 100.00% | 25.273 sec (0.61 MB/s) | 30.82% | 12.392 sec (1.24 MB/s) | 99.87% |
| fib_s14930352 | 14930352 | **3.682 sec (4.06 MB/s)** | 100.00% | 6.112 sec (2.44 MB/s) | 96.99% | 4.503 sec (3.32 MB/s) | 100.00% |
| fss10 | 12078908 | **3.123 sec (3.87 MB/s)** | 100.00% | 5.428 sec (2.23 MB/s) | 95.65% | 3.771 sec (3.20 MB/s) | 100.00% |
| fss9 | 2851443 | **0.757 sec (3.77 MB/s)** | 100.00% | 1.277 sec (2.23 MB/s) | 95.65% | 0.891 sec (3.20 MB/s) | 100.00% |
| houston | 3839141 | 0.788 sec (4.87 MB/s) | 100.00% | **0.502 sec (7.65 MB/s)** | 97.64% | 0.815 sec (4.71 MB/s) | 100.00% |
| paper5x80 | 956322 | **0.179 sec (5.35 MB/s)** | 100.00% | 0.192 sec (4.97 MB/s) | 94.59% | 0.217 sec (4.41 MB/s) | 99.94% |
| test1 | 2097152 | 0.481 sec (4.36 MB/s) | 100.00% | **0.245 sec (8.55 MB/s)** | 100.00% | 0.494 sec (4.24 MB/s) | 100.00% |
| test2 | 2097152 | 0.479 sec (4.38 MB/s) | 100.00% | **0.246 sec (8.54 MB/s)** | 100.00% | 0.494 sec (4.24 MB/s) | 100.00% |
| test3 | 2097088 | **0.300 sec (6.98 MB/s)** | 100.00% | 0.461 sec (4.55 MB/s) | 100.00% | 0.506 sec (4.15 MB/s) | 100.00% |
