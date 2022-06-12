# The esa-matchfinder

The esa-matchfinder is a C99 library for efficient Lempel-Ziv factorization using enhanced suffix array (ESA).

Copyright (c) 2022 Ilya Grebnov <ilya.grebnov@gmail.com>

> * The esa-matchfinder is block based algorithm with maximum supported block size of 512 megabytes finding matches in range of 2..64 bytes using 12x bytes of memory. Note, ESA_MATCHFINDER_MATCH_BITS definition could be updated to support larger match finding range, but with corresponding reduction in maximum supported block size.
> * The esa-matchfinder is fast in best, average and worst cases (see [Benchmarks](#benchmarks) below). But the esa-matchfinder is sensitive to fast memory and software prefetching and might not be suitable for some CPU architectures. The esa-matchfinder might also be slower than other algorithms on specific data types (DNA sequences in particular), so please benchmark yourself.
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
  * Input files were capped at 64MB for medium distance tests and 510MB for long distance tests
  * MMC source code was modified to increase dictionary size to 64MB
  * For all match finders maximum match length was set to 64 bytes (other parameters were not changed)
  * The timings are minimum of five runs measuring single-threaded performance in *optimal* parsing mode

## Specification (x86-64 architecture) ##
  * OS: Microsoft Windows 10 Pro 64-Bit
  * CPU: Intel Core i7-9700K Processor (12M Cache, 5GHz)
  * RAM: 2x8 GB dual-channel DDR4 (4133 MHz, 17-17-17-37)
  * Compiler: Microsoft Visual C++ compiler v14.32 
  * Optimizations: /MD /DNDEBUG /O2 /GL /arch:AVX2

### Silesia Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | MMC v0.2.0 | LZMA BT4 v21.07 |
|:----:|:----:|:----------------------:|:----------:|:---------------:|
| dickens | 10192446 | **0.729 sec (13.98 MB/s)** | 1.736 sec (5.87 MB/s) | 2.984 sec (3.42 MB/s) |
| mozilla | 51220480 | **3.610 sec (14.19 MB/s)** | 14.278 sec (3.59 MB/s) | 8.782 sec (5.83 MB/s) |
| mr | 9970564 | **0.872 sec (11.43 MB/s)** | 2.517 sec (3.96 MB/s) | 1.787 sec (5.58 MB/s) |
| nci | 33553445 | **3.666 sec (9.15 MB/s)** | 6.939 sec (4.84 MB/s) | 6.103 sec (5.50 MB/s) |
| ooffice | 6152192 | **0.330 sec (18.64 MB/s)** | 1.086 sec (5.67 MB/s) | 0.812 sec (7.58 MB/s) |
| osdb | 10085684 | **0.609 sec (16.56 MB/s)** | 2.926 sec (3.45 MB/s) | 1.770 sec (5.70 MB/s) |
| reymont | 6627202 | **0.487 sec (13.61 MB/s)** | 1.110 sec (5.97 MB/s) | 1.640 sec (4.04 MB/s) |
| samba | 21606400 | **1.509 sec (14.32 MB/s)** | 5.023 sec (4.30 MB/s) | 3.044 sec (7.10 MB/s) |
| sao | 7251944 | **0.414 sec (17.52 MB/s)** | 0.952 sec (7.62 MB/s) | 1.039 sec (6.98 MB/s) |
| webster | 41458703 | **4.514 sec (9.18 MB/s)** | 7.739 sec (5.36 MB/s) | 13.111 sec (3.16 MB/s) |
| x-ray | 8474240 | **0.418 sec (20.27 MB/s)** | 2.076 sec (4.08 MB/s) | 1.140 sec (7.43 MB/s) |
| xml | 5345280 | **0.283 sec (18.89 MB/s)** | 1.050 sec (5.09 MB/s) | 0.677 sec (7.90 MB/s) |

### Large Canterbury Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | MMC v0.2.0 | LZMA BT4 v21.07 |
|:----:|:----:|:----------------------:|:----------:|:---------------:|
| bible.txt | 4047392 | **0.260 sec (15.57 MB/s)** | 0.735 sec (5.51 MB/s) | 0.840 sec (4.82 MB/s) |
| E.coli | 4638690 | **0.331 sec (14.01 MB/s)** | 0.449 sec (10.33 MB/s) | 1.684 sec (2.75 MB/s) |
| world192.txt | 2473400 | **0.140 sec (17.67 MB/s)** | 0.625 sec (3.96 MB/s) | 0.377 sec (6.56 MB/s) |

### Manzini Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | MMC v0.2.0 | LZMA BT4 v21.07 |
|:----:|:----:|:----------------------:|:----------:|:---------------:|
| chr22.dna | 34553758 | 4.792 sec (7.21 MB/s) | **2.431 sec (14.21 MB/s)** | 18.922 sec (1.83 MB/s) |
| etext99 | 67108864 | **7.025 sec (9.55 MB/s)** | 11.336 sec (5.92 MB/s) | 27.841 sec (2.41 MB/s) |
| gcc-3.0.tar | 67108864 | **6.521 sec (10.29 MB/s)** | 31.405 sec (2.14 MB/s) | 14.909 sec (4.50 MB/s) |
| howto | 39422105 | **3.325 sec (11.86 MB/s)** | 7.638 sec (5.16 MB/s) | 11.008 sec (3.58 MB/s) |
| jdk13c | 67108864 | **6.269 sec (10.70 MB/s)** | 59.649 sec (1.13 MB/s) | 9.990 sec (6.72 MB/s) |
| linux-2.4.5.tar | 67108864 | **5.953 sec (11.27 MB/s)** | 25.448 sec (2.64 MB/s) | 14.641 sec (4.58 MB/s) |
| rctail96 | 67108864 | **5.581 sec (12.02 MB/s)** | 2890.772 sec (0.02 MB/s) | 15.210 sec (4.41 MB/s) |
| rfc | 67108864 | **7.624 sec (8.80 MB/s)** | 15.425 sec (4.35 MB/s) | 18.162 sec (3.70 MB/s) |
| sprot34.dat | 67108864 | **6.654 sec (10.09 MB/s)** | 19.633 sec (3.42 MB/s) | 16.409 sec (4.09 MB/s) |
| w3c2 | 67108864 | **6.057 sec (11.08 MB/s)** | 46.052 sec (1.46 MB/s) | 10.667 sec (6.29 MB/s) |

### Large Text Compression Benchmark Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | MMC v0.2.0 | LZMA BT4 v21.07 |
|:----:|:----:|:----------------------:|:----------:|:---------------:|
| enwik8 | 67108864 | **6.412 sec (10.47 MB/s)** | 13.014 sec (5.16 MB/s) | 25.787 sec (2.60 MB/s) |
| enwik9 | 67108864 | **6.430 sec (10.44 MB/s)** | 13.004 sec (5.16 MB/s) | 25.746 sec (2.61 MB/s) |

### The Gauntlet Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | MMC v0.2.0 | LZMA BT4 v21.07 |
|:----:|:----:|:----------------------:|:----------:|:---------------:|
| abac | 200000 | 0.013 sec (15.38 MB/s) | 29.413 sec (0.01 MB/s) | **0.011 sec (18.18 MB/s)** |
| abba | 10500596 | 1.045 sec (10.05 MB/s) | **1.010 sec (10.40 MB/s)** | 2.853 sec (3.68 MB/s) |
| book1x20 | 15375420 | **0.965 sec (15.93 MB/s)** | 2.871 sec (5.36 MB/s) | 3.135 sec (4.90 MB/s) |
| fib_s14930352 | 14930352 | **0.837 sec (17.84 MB/s)** | 0.849 sec (17.59 MB/s) | 1.038 sec (14.38 MB/s) |
| fss10 | 12078908 | **0.670 sec (18.03 MB/s)** | 1.143 sec (10.57 MB/s) | 0.988 sec (12.23 MB/s) |
| fss9 | 2851443 | **0.129 sec (22.10 MB/s)** | 0.424 sec (6.73 MB/s) | 0.232 sec (12.29 MB/s) |
| houston | 3839141 | 0.259 sec (14.82 MB/s) | 415.787 sec (0.01 MB/s) | **0.201 sec (19.10 MB/s)** |
| paper5x80 | 956322 | **0.036 sec (26.56 MB/s)** | 0.580 sec (1.65 MB/s) | 0.074 sec (12.92 MB/s) |
| test1 | 2097152 | **0.080 sec (26.21 MB/s)** | 134.310 sec (0.02 MB/s) | 0.127 sec (16.51 MB/s) |
| test2 | 2097152 | **0.079 sec (26.55 MB/s)** | 134.844 sec (0.02 MB/s) | 0.126 sec (16.64 MB/s) |
| test3 | 2097088 | **0.075 sec (27.96 MB/s)** | 1.263 sec (1.66 MB/s) | 0.081 sec (25.89 MB/s) |

### Pizza & Chilli Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | MMC v0.2.0 | LZMA BT4 v21.07 |
|:----:|:----:|:----------------------:|:----------:|:---------------:|
| dblp.xml | 67108864 | **6.846 sec (9.80 MB/s)** | 23.530 sec (2.85 MB/s) | 14.313 sec (4.69 MB/s) |
| dna | 67108864 | 10.188 sec (6.59 MB/s) | **4.701 sec (14.28 MB/s)** | 44.323 sec (1.51 MB/s) |
| english.1024MB | 67108864 | **7.543 sec (8.90 MB/s)** | 11.462 sec (5.85 MB/s) | 25.163 sec (2.67 MB/s) |
| pitches | 55832855 | **3.762 sec (14.84 MB/s)** | 10.213 sec (5.47 MB/s) | 10.658 sec (5.24 MB/s) |
| proteins | 67108864 | **5.996 sec (11.19 MB/s)** | 21.488 sec (3.12 MB/s) | 32.817 sec (2.04 MB/s) |
| sources | 67108864 | **5.937 sec (11.30 MB/s)** | 11.738 sec (5.72 MB/s) | 16.119 sec (4.16 MB/s) |

### Pizza & Chilli Repetitive Corpus (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | MMC v0.2.0 | LZMA BT4 v21.07 |
|:----:|:----:|:----------------------:|:----------:|:---------------:|
| cere | 67108864 | 7.822 sec (8.58 MB/s) | **4.505 sec (14.90 MB/s)** | 38.015 sec (1.77 MB/s) |
| coreutils | 67108864 | **5.670 sec (11.84 MB/s)** | 18.372 sec (3.65 MB/s) | 12.600 sec (5.33 MB/s) |
| einstein.de.txt | 67108864 | **4.640 sec (14.46 MB/s)** | 124.802 sec (0.54 MB/s) | 6.549 sec (10.25 MB/s) |
| einstein.en.txt | 67108864 | **4.667 sec (14.38 MB/s)** | 145.440 sec (0.46 MB/s) | 6.515 sec (10.30 MB/s) |
| Escherichia_Coli | 67108864 | 8.791 sec (7.63 MB/s) | **4.028 sec (16.66 MB/s)** | 38.775 sec (1.73 MB/s) |
| influenza | 67108864 | 16.802 sec (3.99 MB/s) | **4.361 sec (15.39 MB/s)** | 12.353 sec (5.43 MB/s) |
| kernel | 67108864 | **5.174 sec (12.97 MB/s)** | 13.090 sec (5.13 MB/s) | 15.134 sec (4.43 MB/s) |
| para | 67108864 | 8.682 sec (7.73 MB/s) | **4.463 sec (15.04 MB/s)** | 40.073 sec (1.67 MB/s) |
| world_leaders | 46968181 | **5.755 sec (8.16 MB/s)** | 73.836 sec (0.64 MB/s) | 6.437 sec (7.30 MB/s) |
| dblp.xml.00001.1 | 67108864 | **6.420 sec (10.45 MB/s)** | 24.752 sec (2.71 MB/s) | 10.972 sec (6.12 MB/s) |
| dblp.xml.00001.2 | 67108864 | **6.695 sec (10.02 MB/s)** | 24.660 sec (2.72 MB/s) | 11.007 sec (6.10 MB/s) |
| dblp.xml.0001.1 | 67108864 | **7.787 sec (8.62 MB/s)** | 24.644 sec (2.72 MB/s) | 10.987 sec (6.11 MB/s) |
| dblp.xml.0001.2 | 67108864 | **8.720 sec (7.70 MB/s)** | 22.591 sec (2.97 MB/s) | 11.315 sec (5.93 MB/s) |
| dna.001.1 | 67108864 | 8.583 sec (7.82 MB/s) | **5.563 sec (12.06 MB/s)** | 22.291 sec (3.01 MB/s) |
| english.001.2 | 67108864 | **6.531 sec (10.28 MB/s)** | 16.256 sec (4.13 MB/s) | 14.501 sec (4.63 MB/s) |
| proteins.001.1 | 67108864 | **6.758 sec (9.93 MB/s)** | 35.937 sec (1.87 MB/s) | 8.630 sec (7.78 MB/s) |
| sources.001.2 | 67108864 | **6.246 sec (10.74 MB/s)** | 16.539 sec (4.06 MB/s) | 11.051 sec (6.07 MB/s) |
| fib41 | 67108864 | 4.459 sec (15.05 MB/s) | **3.127 sec (21.46 MB/s)** | 4.667 sec (14.38 MB/s) |
| rs.13 | 67108864 | **4.422 sec (15.18 MB/s)** | 5.475 sec (12.26 MB/s) | 5.475 sec (12.26 MB/s) |
| tm29 | 67108864 | 4.577 sec (14.66 MB/s) | **2.427 sec (27.65 MB/s)** | 6.700 sec (10.02 MB/s) |

### Long Distance Corpus (large files from above) (x86-64 architecture) ###

| file | size | esa-matchfinder v1.0.0 | LZMA BT4 v21.07 |
|:----:|:----:|:----------------------:|:---------------:|
| etext99 | 105277340 | **12.017 sec (8.76 MB/s)** | 46.514 sec (2.26 MB/s) |
| gcc-3.0.tar | 86630400 | **8.193 sec (10.57 MB/s)** | 17.775 sec (4.87 MB/s) |
| jdk13c | 69728899 | **6.483 sec (10.76 MB/s)** | 10.214 sec (6.83 MB/s) |
| linux-2.4.5.tar | 116254720 | **11.063 sec (10.51 MB/s)** | 26.948 sec (4.31 MB/s) |
| rctail96 | 114711151 | **10.595 sec (10.83 MB/s)** | 27.143 sec (4.23 MB/s) |
| rfc | 116421901 | **14.592 sec (7.98 MB/s)** | 34.368 sec (3.39 MB/s) |
| sprot34.dat | 109617186 | **11.570 sec (9.47 MB/s)** | 28.412 sec (3.86 MB/s) |
| w3c2 | 104201579 | **9.449 sec (11.03 MB/s)** | 15.009 sec (6.94 MB/s) |
| enwik8 | 100000000 | **10.187 sec (9.82 MB/s)** | 40.141 sec (2.49 MB/s) |
| enwik9 | 534773760 | **81.368 sec (6.57 MB/s)** | 257.039 sec (2.08 MB/s) |
| dblp.xml | 296135874 | **38.583 sec (7.68 MB/s)** | 74.522 sec (3.97 MB/s) |
| dna | 403927746 | **98.562 sec (4.10 MB/s)** | 304.514 sec (1.33 MB/s) |
| english.1024MB | 534773760 | **93.898 sec (5.70 MB/s)** | 318.262 sec (1.68 MB/s) |
| proteins | 534773760 | **70.014 sec (7.64 MB/s)** | 319.924 sec (1.67 MB/s) |
| sources | 210866607 | **22.047 sec (9.56 MB/s)** | 56.558 sec (3.73 MB/s) |
| cere | 461286644 | **67.885 sec (6.80 MB/s)** | 266.890 sec (1.73 MB/s) |
| coreutils | 205281778 | **18.859 sec (10.89 MB/s)** | 39.999 sec (5.13 MB/s) |
| einstein.de.txt | 92758441 | **6.553 sec (14.16 MB/s)** | 9.088 sec (10.21 MB/s) |
| einstein.en.txt | 467626544 | **39.140 sec (11.95 MB/s)** | 48.841 sec (9.57 MB/s) |
| Escherichia_Coli | 112689515 | **15.994 sec (7.05 MB/s)** | 63.415 sec (1.78 MB/s) |
| influenza | 154808555 | 44.249 sec (3.50 MB/s) | **27.765 sec (5.58 MB/s)** |
| kernel | 257961616 | **21.829 sec (11.82 MB/s)** | 59.935 sec (4.30 MB/s) |
| para | 429265758 | **69.957 sec (6.14 MB/s)** | 260.698 sec (1.65 MB/s) |
| dblp.xml.00001.1 | 104857600 | **10.742 sec (9.76 MB/s)** | 16.543 sec (6.34 MB/s) |
| dblp.xml.00001.2 | 104857600 | **10.843 sec (9.67 MB/s)** | 16.651 sec (6.30 MB/s) |
| dblp.xml.0001.1 | 104857600 | **13.253 sec (7.91 MB/s)** | 16.623 sec (6.31 MB/s) |
| dblp.xml.0001.2 | 104857600 | **15.713 sec (6.67 MB/s)** | 17.324 sec (6.05 MB/s) |
| dna.001.1 | 104857600 | **16.521 sec (6.35 MB/s)** | 33.202 sec (3.16 MB/s) |
| english.001.2 | 104857600 | **11.010 sec (9.52 MB/s)** | 21.430 sec (4.89 MB/s) |
| proteins.001.1 | 104857600 | **12.536 sec (8.36 MB/s)** | 12.564 sec (8.35 MB/s) |
| sources.001.2 | 104857600 | **10.417 sec (10.07 MB/s)** | 16.700 sec (6.28 MB/s) |
| fib41 | 267914296 | 19.807 sec (13.53 MB/s) | **18.883 sec (14.19 MB/s)** |
| rs.13 | 216747218 | **16.162 sec (13.41 MB/s)** | 17.746 sec (12.21 MB/s) |
| tm29 | 268435456 | **19.737 sec (13.60 MB/s)** | 26.593 sec (10.09 MB/s) |

## Specification (ARMv8 architecture) ##
  * OS: Ubuntu 20.04 LTS 64-Bit
  * CPU: ODROID-N2+ Amlogic S922X (Cortex-A73 2.4Ghz)
  * RAM: 4GB LPDDR4 (2666 MHz) 
  * Compiler: Clang v10.0.0
  * Optimizations: -DNDEBUG -O3 -flto=thin -mcpu=native

### Silesia Corpus (ARMv8 architecture) ###

| file | size | esa-matchfinder v1.0.0 | MMC v0.2.0 | LZMA BT4 v21.07 |
|:----:|:----:|:----------------------:|:----------:|:---------------:|
| dickens | 10192446 | **4.968 sec (2.05 MB/s)** | 5.531 sec (1.84 MB/s) | 10.160 sec (1.00 MB/s) |
| mozilla | 51220480 | **19.923 sec (2.57 MB/s)** | 45.908 sec (1.12 MB/s) | 26.223 sec (1.95 MB/s) |
| mr | 9970564 | **4.253 sec (2.34 MB/s)** | 10.262 sec (0.97 MB/s) | 6.155 sec (1.62 MB/s) |
| nci | 33553445 | 20.503 sec (1.64 MB/s) | 21.555 sec (1.56 MB/s) | **17.308 sec (1.94 MB/s)** |
| ooffice | 6152192 | **2.020 sec (3.05 MB/s)** | 3.232 sec (1.90 MB/s) | 2.926 sec (2.10 MB/s) |
| osdb | 10085684 | **4.032 sec (2.50 MB/s)** | 9.461 sec (1.07 MB/s) | 5.731 sec (1.76 MB/s) |
| reymont | 6627202 | 3.082 sec (2.15 MB/s) | **2.971 sec (2.23 MB/s)** | 4.965 sec (1.33 MB/s) |
| samba | 21606400 | **8.050 sec (2.68 MB/s)** | 17.298 sec (1.25 MB/s) | 8.933 sec (2.42 MB/s) |
| sao | 7251944 | **2.520 sec (2.88 MB/s)** | 2.701 sec (2.68 MB/s) | 3.685 sec (1.97 MB/s) |
| webster | 41458703 | 24.075 sec (1.72 MB/s) | **22.343 sec (1.86 MB/s)** | 36.771 sec (1.13 MB/s) |
| x-ray | 8474240 | **2.226 sec (3.81 MB/s)** | 5.884 sec (1.44 MB/s) | 3.935 sec (2.15 MB/s) |
| xml | 5345280 | **1.594 sec (3.35 MB/s)** | 2.676 sec (2.00 MB/s) | 1.912 sec (2.80 MB/s) |

### Large Canterbury Corpus (ARMv8 architecture) ###

| file | size | esa-matchfinder v1.0.0 | MMC v0.2.0 | LZMA BT4 v21.07 |
|:----:|:----:|:----------------------:|:----------:|:---------------:|
| bible.txt | 4047392 | **1.698 sec (2.38 MB/s)** | 1.831 sec (2.21 MB/s) | 2.800 sec (1.45 MB/s) |
| E.coli | 4638690 | 2.140 sec (2.17 MB/s) | **0.893 sec (5.19 MB/s)** | 5.889 sec (0.79 MB/s) |
| world192.txt | 2473400 | **0.955 sec (2.59 MB/s)** | 1.723 sec (1.44 MB/s) | 1.291 sec (1.92 MB/s) |

### Manzini Corpus (ARMv8 architecture) ###

| file | size | esa-matchfinder v1.0.0 | MMC v0.2.0 | LZMA BT4 v21.07 |
|:----:|:----:|:----------------------:|:----------:|:---------------:|
| chr22.dna | 34553758 | 24.481 sec (1.41 MB/s) | **5.581 sec (6.19 MB/s)** | 55.387 sec (0.62 MB/s) |
| etext99 | 67108864 | 37.418 sec (1.79 MB/s) | **34.982 sec (1.92 MB/s)** | 79.987 sec (0.84 MB/s) |
| gcc-3.0.tar | 67108864 | **32.636 sec (2.06 MB/s)** | 100.431 sec (0.67 MB/s) | 42.323 sec (1.59 MB/s) |
| howto | 39422105 | **18.345 sec (2.15 MB/s)** | 22.167 sec (1.78 MB/s) | 32.921 sec (1.20 MB/s) |
| jdk13c | 67108864 | 30.748 sec (2.18 MB/s) | 186.255 sec (0.36 MB/s) | **28.664 sec (2.34 MB/s)** |
| linux-2.4.5.tar | 67108864 | **30.195 sec (2.22 MB/s)** | 82.927 sec (0.81 MB/s) | 40.778 sec (1.65 MB/s) |
| rctail96 | 67108864 | **30.204 sec (2.22 MB/s)** | N/A | 44.304 sec (1.51 MB/s) |
| rfc | 67108864 | **38.302 sec (1.75 MB/s)** | 45.890 sec (1.46 MB/s) | 52.394 sec (1.28 MB/s) |
| sprot34.dat | 67108864 | **34.425 sec (1.95 MB/s)** | 60.810 sec (1.10 MB/s) | 46.422 sec (1.45 MB/s) |
| w3c2 | 67108864 | 31.657 sec (2.12 MB/s) | 145.700 sec (0.46 MB/s) | **31.388 sec (2.14 MB/s)** |

### Large Text Compression Benchmark Corpus (ARMv8 architecture) ###

| file | size | esa-matchfinder v1.0.0 | MMC v0.2.0 | LZMA BT4 v21.07 |
|:----:|:----:|:----------------------:|:----------:|:---------------:|
| enwik8 | 67108864 | **34.192 sec (1.96 MB/s)** | 39.068 sec (1.72 MB/s) | 71.199 sec (0.94 MB/s) |

### Long Distance Corpus (large files from above) (ARMv8 architecture) ###

| file | size | esa-matchfinder v1.0.0 | LZMA BT4 v21.07 |
|:----:|:----:|:----------------------:|:---------------:|
| etext99 | 105277340 | **62.856 sec (1.67 MB/s)** | 133.241 sec (0.79 MB/s) |
| gcc-3.0.tar | 86630400 | **41.193 sec (2.10 MB/s)** | 51.308 sec (1.69 MB/s) |
| jdk13c | 69728899 | 32.292 sec (2.16 MB/s) | **30.031 sec (2.32 MB/s)** |
| linux-2.4.5.tar | 116254720 | **56.905 sec (2.04 MB/s)** | 75.345 sec (1.54 MB/s) |
| rctail96 | 114711151 | **56.626 sec (2.03 MB/s)** | 80.160 sec (1.43 MB/s) |
| rfc | 116421901 | **75.687 sec (1.54 MB/s)** | 96.740 sec (1.20 MB/s) |
| sprot34.dat | 109617186 | **60.338 sec (1.82 MB/s)** | 80.496 sec (1.36 MB/s) |
| w3c2 | 104201579 | 49.771 sec (2.09 MB/s) | **45.112 sec (2.31 MB/s)** |
| enwik8 | 100000000 | **54.221 sec (1.84 MB/s)** | 113.541 sec (0.88 MB/s) |