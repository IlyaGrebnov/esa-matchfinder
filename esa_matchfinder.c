/*--

This file is a part of esa-matchfinder, a library for efficient
Lempel-Ziv factorization using enhanced suffix array (ESA).

   Copyright (c) 2022 Ilya Grebnov <ilya.grebnov@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Please see the file LICENSE for full copyright and license details.

--*/

// This file uses the libsais library for linear time suffix array (SA)
// and permuted longest common prefix array (PLCP) construction.
//
// See https://github.com/IlyaGrebnov/libsais for more information.
//
// The libsais library is released under Apache License 2.0.
// Copyright (c) 2021-2022 Ilya Grebnov <ilya.grebnov@gmail.com>
//

#include "esa_matchfinder.h"
#include "libsais/libsais.h"

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#if defined(_OPENMP)
    #include <omp.h>

    #define ESA_MF_NUM_THREADS_MAX      (256)
#else
    #define ESA_MF_UNUSED(_x)           (void)(_x)
    #define ESA_MF_NUM_THREADS_MAX      (1)
#endif

#define ESA_MF_TOTAL_BITS               (64)

#define ESA_MF_LCP_BITS                 (ESA_MATCHFINDER_MATCH_BITS)
#define ESA_MF_LCP_MAX                  (((uint64_t)1 << ESA_MF_LCP_BITS) - 1)
#define ESA_MF_LCP_SHIFT                (ESA_MF_TOTAL_BITS - ESA_MF_LCP_BITS)
#define ESA_MF_LCP_MASK                 (ESA_MF_LCP_MAX << ESA_MF_LCP_SHIFT)

#define ESA_MF_OFFSET_BITS              (ESA_MF_LCP_SHIFT / 2)
#define ESA_MF_OFFSET_MAX               (((uint64_t)1 << ESA_MF_OFFSET_BITS) - 1)
#define ESA_MF_OFFSET_SHIFT             (ESA_MF_TOTAL_BITS - ESA_MF_LCP_BITS - ESA_MF_OFFSET_BITS)
#define ESA_MF_OFFSET_MASK              (ESA_MF_OFFSET_MAX << ESA_MF_OFFSET_SHIFT)

#define ESA_MF_PARENT_BITS              (ESA_MF_OFFSET_SHIFT)
#define ESA_MF_PARENT_MAX               (((uint64_t)1 << ESA_MF_PARENT_BITS) - 1)
#define ESA_MF_PARENT_SHIFT             (ESA_MF_TOTAL_BITS - ESA_MF_LCP_BITS - ESA_MF_OFFSET_BITS - ESA_MF_PARENT_BITS)
#define ESA_MF_PARENT_MASK              (ESA_MF_PARENT_MAX << ESA_MF_PARENT_SHIFT)

#define ESA_MF_STORAGE_PADDING          (64)

#if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunreachable-code"
    #pragma clang diagnostic ignored "-Wstrict-aliasing"
    #pragma clang diagnostic ignored "-Wuninitialized"
#elif defined(__GNUC__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunreachable-code"
    #pragma GCC diagnostic ignored "-Wstrict-aliasing"
    #pragma GCC diagnostic ignored "-Wuninitialized"
#elif defined(_MSC_VER)
    #pragma warning(push)
    #pragma warning(disable: 4127)
    #pragma warning(disable: 4820)
#endif

typedef struct ESA_MF_THREAD_STATE
{
    ptrdiff_t               interval_tree_start;
    ptrdiff_t               interval_tree_end;
} ESA_MF_THREAD_STATE;

typedef struct ESA_MF_CONTEXT
{
    uint64_t                prefetch[4][8];
    uint64_t                position;

    uint64_t *              sa_parent_link;
    uint32_t *              plcp_leaf_link;
    uint64_t                min_match_length_minus_1;

    int32_t *               esa_storage;
    void *                  libsais_ctx;

    int32_t                 block_size;
    int32_t                 max_block_size;
    int32_t                 min_match_length;
    int32_t                 max_match_length;
    int32_t                 num_threads;

    ESA_MF_THREAD_STATE     threads[ESA_MF_NUM_THREADS_MAX];
} ESA_MF_CONTEXT;

#if defined(__GNUC__) || defined(__clang__)
    #define ESA_MF_RESTRICT __restrict__
#elif defined(_MSC_VER) || defined(__INTEL_COMPILER)
    #define ESA_MF_RESTRICT __restrict
#else
    #error Your compiler, configuration or platform is not supported.
#endif

#if defined(__has_builtin)
    #if __has_builtin(__builtin_prefetch)
        #define ESA_MF_HAS_BUILTIN_PREFECTCH
    #endif
#elif defined(__GNUC__) && (((__GNUC__ == 3) && (__GNUC_MINOR__ >= 2)) || (__GNUC__ >= 4))
    #define ESA_MF_HAS_BUILTIN_PREFECTCH
#endif

#if defined(ESA_MF_HAS_BUILTIN_PREFECTCH)
    #define esa_matchfinder_prefetchr(address) __builtin_prefetch((const void *)(address), 0, 0)
    #define esa_matchfinder_prefetchw(address) __builtin_prefetch((const void *)(address), 1, 0)
#elif defined (_M_IX86) || defined (_M_AMD64)
    #include <intrin.h>
    #define esa_matchfinder_prefetchr(address) _mm_prefetch((const void *)(address), _MM_HINT_NTA)
    #define esa_matchfinder_prefetchw(address) _m_prefetchw((const void *)(address))
#elif defined (_M_ARM)
    #include <intrin.h>
    #define esa_matchfinder_prefetchr(address) __prefetch((const void *)(address))
    #define esa_matchfinder_prefetchw(address) __prefetchw((const void *)(address))
#elif defined (_M_ARM64)
    #include <intrin.h>
    #define esa_matchfinder_prefetchr(address) __prefetch2((const void *)(address), 1)
    #define esa_matchfinder_prefetchw(address) __prefetch2((const void *)(address), 17)
#else
    #error Your compiler, configuration or platform is not supported.
#endif

#if !defined(__LITTLE_ENDIAN__) && !defined(__BIG_ENDIAN__)
    #if defined(_LITTLE_ENDIAN) \
            || (defined(BYTE_ORDER) && defined(LITTLE_ENDIAN) && BYTE_ORDER == LITTLE_ENDIAN) \
            || (defined(_BYTE_ORDER) && defined(_LITTLE_ENDIAN) && _BYTE_ORDER == _LITTLE_ENDIAN) \
            || (defined(__BYTE_ORDER) && defined(__LITTLE_ENDIAN) && __BYTE_ORDER == __LITTLE_ENDIAN) \
            || (defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
        #define __LITTLE_ENDIAN__
    #elif defined(_BIG_ENDIAN) \
            || (defined(BYTE_ORDER) && defined(BIG_ENDIAN) && BYTE_ORDER == BIG_ENDIAN) \
            || (defined(_BYTE_ORDER) && defined(_BIG_ENDIAN) && _BYTE_ORDER == _BIG_ENDIAN) \
            || (defined(__BYTE_ORDER) && defined(__BIG_ENDIAN) && __BYTE_ORDER == __BIG_ENDIAN) \
            || (defined(__BYTE_ORDER__) && defined(__ORDER_BIG_ENDIAN__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
        #define __BIG_ENDIAN__
    #elif defined(_WIN32)
        #define __LITTLE_ENDIAN__
    #else
        #error Your compiler, configuration or platform is not supported.
    #endif
#endif

static void * esa_matchfinder_align_up(const void * address, size_t alignment)
{
    return (void *)((((ptrdiff_t)address) + ((ptrdiff_t)alignment) - 1) & (-((ptrdiff_t)alignment)));
}

static void * esa_matchfinder_alloc_aligned(size_t size, size_t alignment)
{
    void * address = malloc(size + sizeof(short) + alignment - 1);
    if (address != NULL)
    {
        void * aligned_address = esa_matchfinder_align_up((void *)((ptrdiff_t)address + (ptrdiff_t)(sizeof(short))), alignment);
        ((short *)aligned_address)[-1] = (short)((ptrdiff_t)aligned_address - (ptrdiff_t)address);

        return aligned_address;
    }

    return NULL;
}

static void esa_matchfinder_free_aligned(void * aligned_address)
{
    if (aligned_address != NULL)
    {
        free((void *)((ptrdiff_t)aligned_address - ((short *)aligned_address)[-1]));
    }
}

static void esa_matchfinder_set_position(ESA_MF_CONTEXT * matchfinder_ctx, uint64_t position)
{
    matchfinder_ctx->position = position;
    memset(matchfinder_ctx->prefetch, 0, sizeof(matchfinder_ctx->prefetch));
}

static ESA_MF_CONTEXT * esa_matchfinder_alloc_ctx(int32_t max_block_size, int32_t min_match_length, int32_t max_match_length, int32_t num_threads)
{
    num_threads                             = num_threads < ESA_MF_NUM_THREADS_MAX ? num_threads : ESA_MF_NUM_THREADS_MAX;
    max_block_size                          = (max_block_size + ESA_MF_STORAGE_PADDING - 1) & (-ESA_MF_STORAGE_PADDING);

    ESA_MF_CONTEXT *    matchfinder_ctx     = (ESA_MF_CONTEXT *)esa_matchfinder_alloc_aligned(sizeof(ESA_MF_CONTEXT), ESA_MF_STORAGE_PADDING);
    int32_t *           esa_storage         = (int32_t *)esa_matchfinder_alloc_aligned((2 * ESA_MF_STORAGE_PADDING + 3 * (size_t)max_block_size) * sizeof(int32_t), ESA_MF_STORAGE_PADDING);

#if defined(_OPENMP)
    void *              libsais_ctx         = libsais_create_ctx_omp(num_threads);
#else
    void *              libsais_ctx         = libsais_create_ctx();
#endif

    if (matchfinder_ctx != NULL && esa_storage != NULL && libsais_ctx != NULL)
    {
        matchfinder_ctx->esa_storage                = esa_storage;
        matchfinder_ctx->libsais_ctx                = libsais_ctx;

        matchfinder_ctx->block_size                 = -1;
        matchfinder_ctx->max_block_size             = max_block_size;
        matchfinder_ctx->min_match_length           = min_match_length;
        matchfinder_ctx->max_match_length           = max_match_length;
        matchfinder_ctx->num_threads                = num_threads;

        matchfinder_ctx->sa_parent_link             = (uint64_t *)(void *)(matchfinder_ctx->esa_storage + ESA_MF_STORAGE_PADDING) + 0 * matchfinder_ctx->max_block_size;
        matchfinder_ctx->plcp_leaf_link             = (uint32_t *)(void *)(matchfinder_ctx->esa_storage + ESA_MF_STORAGE_PADDING) + 2 * matchfinder_ctx->max_block_size;
        matchfinder_ctx->min_match_length_minus_1   = (uint64_t)matchfinder_ctx->min_match_length - 1;

        esa_matchfinder_set_position(matchfinder_ctx, (uint64_t)-1);

        return matchfinder_ctx;
    }

    libsais_free_ctx(libsais_ctx);

    esa_matchfinder_free_aligned(esa_storage);
    esa_matchfinder_free_aligned(matchfinder_ctx);

    return NULL;
}

static void esa_matchfinder_free_ctx(ESA_MF_CONTEXT * matchfinder_ctx)
{
    if (matchfinder_ctx != NULL)
    {
        libsais_free_ctx(matchfinder_ctx->libsais_ctx);

        esa_matchfinder_free_aligned(matchfinder_ctx->esa_storage);
        esa_matchfinder_free_aligned(matchfinder_ctx);
    }
}

static void esa_matchfinder_convert_right_to_left_32u_to_64u(uint32_t * S, uint64_t * D, ptrdiff_t omp_block_start, ptrdiff_t omp_block_size)
{
    ptrdiff_t i, j; for (i = omp_block_start + omp_block_size - 1, j = omp_block_start; i >= j; i -= 1) { D[i] = (uint64_t)S[i]; }
}

static void esa_matchfinder_convert_left_to_right_32u_to_64u(uint32_t * ESA_MF_RESTRICT S, uint64_t * ESA_MF_RESTRICT D, ptrdiff_t omp_block_start, ptrdiff_t omp_block_size)
{
    ptrdiff_t i, j; for (i = omp_block_start, j = omp_block_start + omp_block_size; i < j; i += 1) { D[i] = (uint64_t)S[i]; }
}

static void esa_matchfinder_convert_inplace_32u_to_64u_omp(uint32_t * S, uint64_t * D, ptrdiff_t n, ptrdiff_t num_threads)
{
    while (n >= 65536)
    {
        ptrdiff_t block_size = n >> 1; n -= block_size;

#if defined(_OPENMP)
        #pragma omp parallel num_threads(num_threads) if(num_threads > 1)
#endif
        {
#if defined(_OPENMP)
            ptrdiff_t omp_thread_num      = omp_get_thread_num();
            ptrdiff_t omp_num_threads     = omp_get_num_threads();
#else
            ESA_MF_UNUSED(num_threads);

            ptrdiff_t omp_thread_num      = 0;
            ptrdiff_t omp_num_threads     = 1;
#endif
            ptrdiff_t omp_block_stride    = (block_size / omp_num_threads) & (-16);
            ptrdiff_t omp_block_start     = omp_thread_num * omp_block_stride;
            ptrdiff_t omp_block_size      = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

            esa_matchfinder_convert_left_to_right_32u_to_64u(S, D, n + omp_block_start, omp_block_size);
        }
    }

    esa_matchfinder_convert_right_to_left_32u_to_64u(S, D, 0, n);
}

static void esa_matchfinder_reset_interval_tree(uint64_t * ESA_MF_RESTRICT sa_parent_link, ptrdiff_t omp_block_start, ptrdiff_t omp_block_size)
{
    ptrdiff_t i, j; for (i = omp_block_start, j = omp_block_start + omp_block_size; i < j; i += 1) { sa_parent_link[i] &= (~ESA_MF_OFFSET_MASK); }
}

static void esa_matchfinder_reset_interval_tree_omp(uint64_t * ESA_MF_RESTRICT sa_parent_link, ptrdiff_t n, ptrdiff_t num_threads)
{
#if defined(_OPENMP)
    #pragma omp parallel num_threads(num_threads) if(num_threads > 1 && n >= 65536)
#endif
    {
#if defined(_OPENMP)
        ptrdiff_t omp_thread_num      = omp_get_thread_num();
        ptrdiff_t omp_num_threads     = omp_get_num_threads();
#else
        ESA_MF_UNUSED(num_threads);

        ptrdiff_t omp_thread_num      = 0;
        ptrdiff_t omp_num_threads     = 1;
#endif
        ptrdiff_t omp_block_stride    = (n / omp_num_threads) & (-16);
        ptrdiff_t omp_block_start     = omp_thread_num * omp_block_stride;
        ptrdiff_t omp_block_size      = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        esa_matchfinder_reset_interval_tree(sa_parent_link, omp_block_start, omp_block_size);
    }
}

static void esa_matchfinder_fast_forward
(
    uint64_t * ESA_MF_RESTRICT  sa_parent_link,
    uint32_t * ESA_MF_RESTRICT  plcp_leaf_link,
    uint64_t                    target_position
)
{
    const uint64_t  prefetch_distance   = 32;
    uint64_t        position            = target_position - 1;

    for (; position >= prefetch_distance; position -= 1)
    {
        esa_matchfinder_prefetchr(&plcp_leaf_link[position - 2 * prefetch_distance]);
        esa_matchfinder_prefetchw(&sa_parent_link[plcp_leaf_link[position - prefetch_distance]]);

        const uint64_t offset           = (uint64_t)position << ESA_MF_OFFSET_SHIFT;
        uint64_t reference              = plcp_leaf_link[position];
        uint64_t interval               = sa_parent_link[reference];

        while ((interval & ESA_MF_OFFSET_MASK) == 0)
        {
            sa_parent_link[reference]   = interval + offset;
            reference                   = (uint32_t)interval;
            interval                    = sa_parent_link[reference];
        }
    }

    for (; position > 0; position -= 1)
    {
        const uint64_t offset           = (uint64_t)position << ESA_MF_OFFSET_SHIFT;
        uint64_t reference              = plcp_leaf_link[position];
        uint64_t interval               = sa_parent_link[reference];

        while ((interval & ESA_MF_OFFSET_MASK) == 0)
        {
            sa_parent_link[reference]   = interval + offset;
            reference                   = (uint32_t)interval;
            interval                    = sa_parent_link[reference];
        }
    }
}

static ptrdiff_t esa_matchfinder_build_interval_tree
(
    uint64_t * ESA_MF_RESTRICT  sa_parent_link,
    uint32_t * ESA_MF_RESTRICT  plcp_leaf_link,
    uint64_t                    min_match_length,
    uint64_t                    max_match_length,
    ptrdiff_t                   omp_block_start,
    ptrdiff_t                   omp_block_size
)
{
    uint64_t intervals[2 * ESA_MATCHFINDER_MAX_MATCH_LENGTH];

    const ptrdiff_t             prefetch_distance       = 32;
    uint64_t * ESA_MF_RESTRICT  stack                   = intervals;
    uint64_t                    top_interval            = stack[0] = 0;
    uint64_t                    next_interval_index     = (uint64_t)(omp_block_start + omp_block_size - 1);

    min_match_length -= 1;
    max_match_length -= min_match_length;

    for (ptrdiff_t i = omp_block_start + omp_block_size - 1; i >= omp_block_start; i -= 1)
    {
        esa_matchfinder_prefetchr(&sa_parent_link[i - 2 * prefetch_distance]);
        esa_matchfinder_prefetchw(&plcp_leaf_link[sa_parent_link[i - prefetch_distance]]);

        uint64_t next_pos                   =  sa_parent_link[i];
        uint64_t next_lcp                   =  (uint64_t)plcp_leaf_link[next_pos] - min_match_length;

        if ((int64_t)next_lcp < 0)          {  next_lcp = 0; }
        if (next_lcp > max_match_length)    {  next_lcp = max_match_length; }

        uint64_t next_interval              =  (next_lcp << ESA_MF_LCP_SHIFT) + next_interval_index;
        uint64_t top_interval_lcp           =  top_interval >> ESA_MF_LCP_SHIFT;

        stack[1]                            =  next_interval;
        top_interval                        =  next_lcp > top_interval_lcp ? next_interval : top_interval;
        next_interval_index                 -= next_lcp > top_interval_lcp;
        stack                               += next_lcp > top_interval_lcp;

        plcp_leaf_link[next_pos]            =  (uint32_t)top_interval;

        while (next_lcp < top_interval_lcp)
        {
            uint64_t closed_interval        =  top_interval;
            
            stack                           =  stack - 1;
            top_interval                    =  stack[0];
            top_interval_lcp                =  top_interval >> ESA_MF_LCP_SHIFT;

            stack[1]                        =  next_interval;
            top_interval                    =  next_lcp > top_interval_lcp ? next_interval : top_interval;
            next_interval_index             -= next_lcp > top_interval_lcp;
            stack                           += next_lcp > top_interval_lcp;
            
            sa_parent_link[(uint32_t)closed_interval] = (uint32_t)top_interval + (closed_interval & ESA_MF_LCP_MASK);
        }
    }

    return (ptrdiff_t)(next_interval_index + 1);
}

#if defined(_OPENMP)

static ptrdiff_t esa_matchfinder_find_breakpoint
(
    uint64_t * ESA_MF_RESTRICT  sa_parent_link,
    uint32_t * ESA_MF_RESTRICT  plcp_leaf_link,
    uint32_t                    min_match_length,
    ptrdiff_t                   omp_block_start,
    ptrdiff_t                   omp_block_size
)
{
    const ptrdiff_t prefetch_distance = 32;

    for (ptrdiff_t i = omp_block_start + omp_block_size - 1; i >= omp_block_start; i -= 1)
    {
        esa_matchfinder_prefetchr(&sa_parent_link[i - 2 * prefetch_distance]);
        esa_matchfinder_prefetchr(&plcp_leaf_link[sa_parent_link[i - prefetch_distance]]);

        if (plcp_leaf_link[sa_parent_link[i]] < min_match_length)
        {
            return i;
        }
    }

    return -1;
}

#endif

static void esa_matchfinder_build_interval_tree_omp
(
    uint64_t * ESA_MF_RESTRICT  sa_parent_link,
    uint32_t * ESA_MF_RESTRICT  plcp_leaf_link,
    uint64_t                    min_match_length,
    uint64_t                    max_match_length,
    ptrdiff_t                   n,
    ptrdiff_t                   num_threads,
    ESA_MF_THREAD_STATE *       threads
)
{
#if defined(_OPENMP)
    ptrdiff_t breakpoints[ESA_MF_NUM_THREADS_MAX];

    for (ptrdiff_t thread = 0; thread < num_threads; thread += 1)
    {
        threads[thread].interval_tree_start = 0;
        threads[thread].interval_tree_end   = 0;
    }

    #pragma omp parallel num_threads(num_threads) if(num_threads > 1 && n >= 65536)
#endif
    {
#if defined(_OPENMP)
        ptrdiff_t omp_thread_num      = omp_get_thread_num();
        ptrdiff_t omp_num_threads     = omp_get_num_threads();
#else
        ESA_MF_UNUSED(num_threads);

        ptrdiff_t omp_thread_num      = 0;
        ptrdiff_t omp_num_threads     = 1;
#endif
        ptrdiff_t omp_block_stride    = (n / omp_num_threads) & (-16);
        ptrdiff_t omp_block_start     = omp_thread_num * omp_block_stride;
        ptrdiff_t omp_block_size      = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;
        ptrdiff_t omp_block_end       = omp_block_start + omp_block_size;

        if (omp_num_threads == 1)
        {
            threads[omp_thread_num].interval_tree_end   = omp_block_end;
            threads[omp_thread_num].interval_tree_start = esa_matchfinder_build_interval_tree(
                sa_parent_link,
                plcp_leaf_link,
                min_match_length,
                max_match_length,
                omp_block_start,
                omp_block_end - omp_block_start);
        }
#if defined(_OPENMP)
        else
        {
            {
                breakpoints[omp_thread_num] = omp_thread_num < omp_num_threads - 1
                    ? esa_matchfinder_find_breakpoint(sa_parent_link, plcp_leaf_link, (uint32_t)min_match_length, omp_block_start, omp_block_end - omp_block_start)
                    : n;
            }

            #pragma omp barrier

            {
                if (breakpoints[omp_thread_num] != -1)
                {
                    omp_block_end       = breakpoints[omp_thread_num];
                    omp_block_start     = 0;

                    for (ptrdiff_t thread = omp_thread_num - 1; thread >= 0; thread -= 1)
                    { 
                        if (breakpoints[thread] != -1) { omp_block_start = breakpoints[thread]; break; }
                    }

                    if (omp_block_start < omp_block_end)
                    {
                        threads[omp_thread_num].interval_tree_end   = omp_block_end;
                        threads[omp_thread_num].interval_tree_start = esa_matchfinder_build_interval_tree(
                            sa_parent_link,
                            plcp_leaf_link,
                            min_match_length,
                            max_match_length,
                            omp_block_start,
                            omp_block_end - omp_block_start);
                    }
                }
            }
        }
#endif
    }

    {
        sa_parent_link[0] = ESA_MF_OFFSET_MASK;
    }
}

void * esa_matchfinder_create(int32_t max_block_size, int32_t min_match_length, int32_t max_match_length)
{
    if ((max_block_size     < 0) ||
        (max_block_size     > ESA_MATCHFINDER_MAX_BLOCK_SIZE) ||
        (min_match_length   < ESA_MATCHFINDER_MIN_MATCH_LENGTH) ||
        (max_match_length   > (int32_t)ESA_MF_LCP_MAX + min_match_length - 1) ||
        (max_match_length   < min_match_length))
    {
        return NULL;
    }

    return (void *)esa_matchfinder_alloc_ctx(max_block_size, min_match_length, max_match_length, 1);
}

#if defined(_OPENMP)

void * esa_matchfinder_create_omp(int32_t max_block_size, int32_t min_match_length, int32_t max_match_length, int32_t num_threads)
{
    if ((max_block_size     < 0) ||
        (max_block_size     > ESA_MATCHFINDER_MAX_BLOCK_SIZE) ||
        (min_match_length   < ESA_MATCHFINDER_MIN_MATCH_LENGTH) ||
        (max_match_length   > (int32_t)ESA_MF_LCP_MAX + min_match_length - 1) ||
        (max_match_length   < min_match_length) ||
        (num_threads        < 0))
    {
        return NULL;
    }

    num_threads = num_threads > 0 ? num_threads : omp_get_max_threads();
    return (void *)esa_matchfinder_alloc_ctx(max_block_size, min_match_length, max_match_length, num_threads);
}

#endif

void esa_matchfinder_destroy(void * mf)
{
    esa_matchfinder_free_ctx((ESA_MF_CONTEXT *)mf);
}

int32_t esa_matchfinder_parse(void * mf, const uint8_t * block, int32_t block_size)
{
    ESA_MF_CONTEXT * matchfinder_ctx = (ESA_MF_CONTEXT *)mf;

    if ((matchfinder_ctx == NULL) || (block == NULL) || (block_size < 0) || (block_size > matchfinder_ctx->max_block_size))
    {
        return ESA_MATCHFINDER_BAD_PARAMETER;
    }

    matchfinder_ctx->block_size = block_size;
    memset(matchfinder_ctx->esa_storage + 0 * ESA_MF_STORAGE_PADDING + 0 * matchfinder_ctx->max_block_size + 0 * matchfinder_ctx->block_size, 0, ESA_MF_STORAGE_PADDING * sizeof(int32_t));
    memset(matchfinder_ctx->esa_storage + 1 * ESA_MF_STORAGE_PADDING + 2 * matchfinder_ctx->max_block_size + 1 * matchfinder_ctx->block_size, 0, ESA_MF_STORAGE_PADDING * sizeof(int32_t));

    int32_t result = libsais_ctx(
        matchfinder_ctx->libsais_ctx,
        block,
        (int32_t *)(void *)matchfinder_ctx->sa_parent_link,
        matchfinder_ctx->block_size,
        (2 * matchfinder_ctx->max_block_size) - matchfinder_ctx->block_size,
        NULL);

    if (result == ESA_MATCHFINDER_NO_ERROR)
    {
#if defined(_OPENMP)
        result = libsais_plcp_omp(
            block,
            (int32_t *)(void *)matchfinder_ctx->sa_parent_link,
            (int32_t *)(void *)matchfinder_ctx->plcp_leaf_link,
            matchfinder_ctx->block_size,
            matchfinder_ctx->num_threads);
#else
        result = libsais_plcp(
            block,
            (int32_t *)(void *)matchfinder_ctx->sa_parent_link,
            (int32_t *)(void *)matchfinder_ctx->plcp_leaf_link,
            block_size);
#endif

        if (result == ESA_MATCHFINDER_NO_ERROR)
        {
            esa_matchfinder_convert_inplace_32u_to_64u_omp(
                (uint32_t *)(void *)matchfinder_ctx->sa_parent_link,
                (uint64_t *)(void *)matchfinder_ctx->sa_parent_link,
                matchfinder_ctx->block_size,
                matchfinder_ctx->num_threads);

            esa_matchfinder_build_interval_tree_omp(
                matchfinder_ctx->sa_parent_link,
                matchfinder_ctx->plcp_leaf_link,
                (uint64_t)matchfinder_ctx->min_match_length,
                (uint64_t)matchfinder_ctx->max_match_length,
                matchfinder_ctx->block_size,
                matchfinder_ctx->num_threads,
                matchfinder_ctx->threads);

            esa_matchfinder_set_position(matchfinder_ctx, 0);
        }
    }

    return result;
}

int32_t esa_matchfinder_get_position(void * mf)
{
    return (int32_t)((ESA_MF_CONTEXT *)mf)->position;
}

int32_t esa_matchfinder_rewind(void * mf, int32_t position)
{
    ESA_MF_CONTEXT * matchfinder_ctx = (ESA_MF_CONTEXT *)mf;

    if ((matchfinder_ctx == NULL) || (position < 0) || (position >= matchfinder_ctx->block_size))
    {
        return ESA_MATCHFINDER_BAD_PARAMETER;
    }

    if (matchfinder_ctx->position != (uint64_t)position)
    {
        if (matchfinder_ctx->position != 0)
        {
            for (ptrdiff_t thread = 0; thread < matchfinder_ctx->num_threads; thread += 1)
            {
                ptrdiff_t interval_tree_start   = matchfinder_ctx->threads[thread].interval_tree_start;
                ptrdiff_t interval_tree_end     = matchfinder_ctx->threads[thread].interval_tree_end;

                if (interval_tree_start < interval_tree_end)
                {
                    esa_matchfinder_reset_interval_tree_omp(
                        matchfinder_ctx->sa_parent_link + interval_tree_start,
                        interval_tree_end - interval_tree_start,
                        matchfinder_ctx->num_threads);
                }
            }
        }

        if (position > 0)
        {
            esa_matchfinder_fast_forward(matchfinder_ctx->sa_parent_link, matchfinder_ctx->plcp_leaf_link, (uint64_t)position);
        }

        esa_matchfinder_set_position(matchfinder_ctx, (uint64_t)position);
    }

    return ESA_MATCHFINDER_NO_ERROR;
}

ESA_MATCHFINDER_MATCH * esa_matchfinder_find_all_matches(void * mf, ESA_MATCHFINDER_MATCH * matches)
{
    ESA_MF_CONTEXT * ESA_MF_RESTRICT const          matchfinder_ctx     = (ESA_MF_CONTEXT *)mf;

    const ptrdiff_t                                 prefetch_distance   = 4;
    const uint64_t                                  position            = matchfinder_ctx->position++;

    uint64_t * ESA_MF_RESTRICT const                sa_parent_link      = matchfinder_ctx->sa_parent_link;
    uint32_t * ESA_MF_RESTRICT const                plcp_leaf_link      = matchfinder_ctx->plcp_leaf_link;
    uint64_t * ESA_MF_RESTRICT const                prefetch            = &matchfinder_ctx->prefetch[position & (prefetch_distance - 1)][0];
    ESA_MATCHFINDER_MATCH * ESA_MF_RESTRICT         next_match          = matches;

    esa_matchfinder_prefetchw(&sa_parent_link[              (sa_parent_link[prefetch[0]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[0] = (sa_parent_link[prefetch[1]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[1] = (sa_parent_link[prefetch[2]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[2] = (sa_parent_link[prefetch[3]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[3] = (sa_parent_link[prefetch[4]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[4] = (sa_parent_link[prefetch[5]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[5] = (sa_parent_link[prefetch[6]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[6] = (plcp_leaf_link[position + 8 * prefetch_distance])]);
    esa_matchfinder_prefetchr(&plcp_leaf_link[position + 9 * prefetch_distance]);

    const uint64_t min_match_length = (uint64_t)matchfinder_ctx->min_match_length_minus_1;
    const uint64_t new_offset       = (uint64_t)position << ESA_MF_OFFSET_SHIFT;
    uint64_t best_match             = ESA_MATCHFINDER_MAX_MATCH_LENGTH;
    uint64_t reference              = plcp_leaf_link[position];

    while (reference != 0)
    {
        const uint64_t interval     = sa_parent_link[reference];
        const uint64_t match        = min_match_length + (interval >> ESA_MF_LCP_SHIFT) + ((interval & ESA_MF_OFFSET_MASK) << (32 - ESA_MF_OFFSET_SHIFT));

#if defined(__LITTLE_ENDIAN__) && !defined(__BIG_ENDIAN__)
        if (offsetof(ESA_MATCHFINDER_MATCH, length) == 0 && offsetof(ESA_MATCHFINDER_MATCH, offset) == 4)
        {
            *(uint64_t *)(void *)next_match = match;
        }
        else
#endif
        {
            next_match->length      = (int32_t)(match      );
            next_match->offset      = (int32_t)(match >> 32);
        }

        next_match                 += match > best_match;
        best_match                  = match;

        sa_parent_link[reference]   = (interval & (~ESA_MF_OFFSET_MASK)) + new_offset;
        reference                   = interval & ESA_MF_PARENT_MASK;
    }

    return next_match;
}

ESA_MATCHFINDER_MATCH esa_matchfinder_find_best_match(void * mf)
{
    ESA_MF_CONTEXT * ESA_MF_RESTRICT const          matchfinder_ctx     = (ESA_MF_CONTEXT *)mf;

    const ptrdiff_t                                 prefetch_distance   = 4;
    const uint64_t                                  position            = matchfinder_ctx->position++;

    uint64_t * ESA_MF_RESTRICT const                sa_parent_link      = matchfinder_ctx->sa_parent_link;
    uint32_t * ESA_MF_RESTRICT const                plcp_leaf_link      = matchfinder_ctx->plcp_leaf_link;
    uint64_t * ESA_MF_RESTRICT const                prefetch            = &matchfinder_ctx->prefetch[position & (prefetch_distance - 1)][0];

    esa_matchfinder_prefetchw(&sa_parent_link[              (sa_parent_link[prefetch[0]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[0] = (sa_parent_link[prefetch[1]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[1] = (sa_parent_link[prefetch[2]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[2] = (sa_parent_link[prefetch[3]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[3] = (sa_parent_link[prefetch[4]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[4] = (sa_parent_link[prefetch[5]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[5] = (sa_parent_link[prefetch[6]] & ESA_MF_PARENT_MASK)]);
    esa_matchfinder_prefetchw(&sa_parent_link[prefetch[6] = (plcp_leaf_link[position + 8 * prefetch_distance])]);
    esa_matchfinder_prefetchr(&plcp_leaf_link[position + 9 * prefetch_distance]);

    const uint64_t min_match_length = (uint64_t)matchfinder_ctx->min_match_length_minus_1;
    const uint64_t new_offset       = (uint64_t)position << ESA_MF_OFFSET_SHIFT;
    uint64_t best_match             = 0;
    uint64_t reference              = plcp_leaf_link[position];

    while (reference != 0)
    {
        const uint64_t interval     = sa_parent_link[reference];
              uint64_t match        = min_match_length + (interval >> ESA_MF_LCP_SHIFT) + ((interval & ESA_MF_OFFSET_MASK) << (32 - ESA_MF_OFFSET_SHIFT));

        match                       = interval & ESA_MF_OFFSET_MASK ? match : best_match;
        best_match                  = best_match == 0               ? match : best_match;

        sa_parent_link[reference]   = (interval & (~ESA_MF_OFFSET_MASK)) + new_offset;
        reference                   = interval & ESA_MF_PARENT_MASK;
    }

    {
        ESA_MATCHFINDER_MATCH match;

#if defined(__LITTLE_ENDIAN__) && !defined(__BIG_ENDIAN__)
        if (offsetof(ESA_MATCHFINDER_MATCH, length) == 0 && offsetof(ESA_MATCHFINDER_MATCH, offset) == 4)
        {
            *(uint64_t *)(void *)&match = best_match;
        }
        else
#endif
        {
            match.length            = (int32_t)(best_match      );
            match.offset            = (int32_t)(best_match >> 32);
        }

        return match;
    }
}

void esa_matchfinder_advance(void * mf, int32_t count)
{
    ESA_MF_CONTEXT * ESA_MF_RESTRICT const          matchfinder_ctx     = (ESA_MF_CONTEXT *)mf;

    const ptrdiff_t                                 prefetch_distance   = 4;
    const uint64_t                                  current_position    = matchfinder_ctx->position;
    const uint64_t                                  target_position     = matchfinder_ctx->position += (uint64_t)count;

    uint64_t * ESA_MF_RESTRICT const                sa_parent_link      = matchfinder_ctx->sa_parent_link;
    uint32_t * ESA_MF_RESTRICT const                plcp_leaf_link      = matchfinder_ctx->plcp_leaf_link;

    for (uint64_t position = current_position; position < target_position; position += 1)
    {
        uint64_t * ESA_MF_RESTRICT const prefetch = &matchfinder_ctx->prefetch[position & (prefetch_distance - 1)][0];

        esa_matchfinder_prefetchw(&sa_parent_link[              (sa_parent_link[prefetch[0]] & ESA_MF_PARENT_MASK)]);
        esa_matchfinder_prefetchw(&sa_parent_link[prefetch[0] = (sa_parent_link[prefetch[1]] & ESA_MF_PARENT_MASK)]);
        esa_matchfinder_prefetchw(&sa_parent_link[prefetch[1] = (sa_parent_link[prefetch[2]] & ESA_MF_PARENT_MASK)]);
        esa_matchfinder_prefetchw(&sa_parent_link[prefetch[2] = (sa_parent_link[prefetch[3]] & ESA_MF_PARENT_MASK)]);
        esa_matchfinder_prefetchw(&sa_parent_link[prefetch[3] = (sa_parent_link[prefetch[4]] & ESA_MF_PARENT_MASK)]);
        esa_matchfinder_prefetchw(&sa_parent_link[prefetch[4] = (sa_parent_link[prefetch[5]] & ESA_MF_PARENT_MASK)]);
        esa_matchfinder_prefetchw(&sa_parent_link[prefetch[5] = (sa_parent_link[prefetch[6]] & ESA_MF_PARENT_MASK)]);
        esa_matchfinder_prefetchw(&sa_parent_link[prefetch[6] = (plcp_leaf_link[position + 8 * prefetch_distance])]);
        esa_matchfinder_prefetchr(&plcp_leaf_link[position + 9 * prefetch_distance]);

        const uint64_t new_offset       = (uint64_t)position << ESA_MF_OFFSET_SHIFT;
        uint64_t reference              = plcp_leaf_link[position];

        while (reference != 0)
        {
            uint64_t interval           = sa_parent_link[reference];

            sa_parent_link[reference]   = (interval & (~ESA_MF_OFFSET_MASK)) + new_offset;
            reference                   = interval & ESA_MF_PARENT_MASK;
        }
    }
}

#if defined(__clang__)
    #pragma clang diagnostic pop
#elif defined(__GNUC__)
    #pragma GCC diagnostic pop
#elif defined(_MSC_VER)
    #pragma warning(pop)
#endif
