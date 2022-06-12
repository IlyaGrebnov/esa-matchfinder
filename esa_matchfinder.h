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

#ifndef ESA_MATCHFINDER_H
#define ESA_MATCHFINDER_H 1

#define ESA_MATCHFINDER_MATCH_BITS          (6)
#define ESA_MATCHFINDER_MAX_BLOCK_SIZE      (1 << ((64 - ESA_MATCHFINDER_MATCH_BITS) / 2))
#define ESA_MATCHFINDER_MIN_MATCH_LENGTH    (2)
#define ESA_MATCHFINDER_MAX_MATCH_LENGTH    (1 << ESA_MATCHFINDER_MATCH_BITS)

#define ESA_MATCHFINDER_NO_ERROR            (0)
#define ESA_MATCHFINDER_BAD_PARAMETER       (-1)

#ifdef __cplusplus
extern "C" {
#endif

    #include <stdint.h>

    typedef struct ESA_MATCHFINDER_MATCH
    {
        int32_t     length;
        int32_t     offset;
    } ESA_MATCHFINDER_MATCH;

    /**
    * Creates the enhanced suffix array (ESA) based match-finder for Lempel-Ziv factorization.
    * @param max_block_size The maximum block size to support (must be less or equal to ESA_MATCHFINDER_MAX_BLOCK_SIZE).
    * @param min_match_length The minimum match length to find (must be greater or equal to ESA_MATCHFINDER_MIN_MATCH_LENGTH).
    * @param max_match_length The maximum match length to find (must be less or equal to ESA_MATCHFINDER_MAX_MATCH_LENGTH).
    * @return The enhanced suffix array (ESA) based match-finder, NULL otherwise.
    */
    void * esa_matchfinder_create(int32_t max_block_size, int32_t min_match_length, int32_t max_match_length);

#if defined(_OPENMP)
    /**
    * Creates the enhanced suffix array (ESA) based match-finder for Lempel-Ziv factorization with multi-threaded optimization using OpenMP.
    * @param max_block_size The maximum block size to support (must be less or equal to ESA_MATCHFINDER_MAX_BLOCK_SIZE).
    * @param min_match_length The minimum match length to find (must be greater or equal to ESA_MATCHFINDER_MIN_MATCH_LENGTH).
    * @param max_match_length The maximum match length to find (must be less or equal to ESA_MATCHFINDER_MAX_MATCH_LENGTH).
    * @param num_threads The number of OpenMP threads to use (can be 0 for default number of OpenMP threads).
    * @return The enhanced suffix array (ESA) based match-finder, NULL otherwise.
    */
    void * esa_matchfinder_create_omp(int32_t max_block_size, int32_t min_match_length, int32_t max_match_length, int32_t num_threads);
#endif

    /**
    * Destroys the match-finder and frees previously allocated memory.
    * @param mf The enhanced suffix array (ESA) based match-finder.
    */
    void esa_matchfinder_destroy(void * mf);

    /**
    * Parses the input block by building enhanced suffix array (ESA) to speed up subsequent match-finding operations.
    * @param mf The enhanced suffix array (ESA) based match-finder.
    * @param block The input block to parse.
    * @param block_size The size of input block to parse.
    * @return 0 if no error occurred, -1 otherwise.
    */
    int32_t esa_matchfinder_parse(void * mf, const uint8_t * block, int32_t block_size);

    /**
    * Gets the current match-finder position.
    * @param mf The enhanced suffix array (ESA) based match-finder.
    * @return The current match-finder position.
    */
    int32_t esa_matchfinder_get_position(void * mf);

    /**
    * Rewinds the match-finder forward or backward to the specified position.
    * @param mf The enhanced suffix array (ESA) based match-finder.
    * @param position The match-finder position to rewind to.
    * @return 0 if no error occurred, -1 otherwise.
    */
    int32_t esa_matchfinder_rewind(void * mf, int32_t position);

    /**
    * Finds all distance optimal matches at the current match-finder position and advances position by one byte. The recorded
    * matches will be sorted by strictly decreasing length and strictly increasing offset from the beginning of the block.
    * @param mf The enhanced suffix array (ESA) based match-finder.
    * @param matches The output array to record the matches (array must be of ESA_MATCHFINDER_MAX_MATCH_LENGTH size).
    * @return The pointer to the end of recorded matches array (if no matches were found, this will be the same as matches).
    */
    ESA_MATCHFINDER_MATCH * esa_matchfinder_find_all_matches(void * mf, ESA_MATCHFINDER_MATCH * matches);

    /**
    * Finds the best match at the current match-finder position and advances position by one byte.
    * @param mf The enhanced suffix array (ESA) based match-finder.
    * @return The best match found (match of zero length and zero offset is returned if no matches were found).
    */
    ESA_MATCHFINDER_MATCH esa_matchfinder_find_best_match(void * mf);

    /**
    * Advances the match-finder position forward by the specified number of bytes without recording matches.
    * @param mf The enhanced suffix array (ESA) based match-finder.
    * @param count The number of bytes to advance.
    */
    void esa_matchfinder_advance(void * mf, int32_t count);

#ifdef __cplusplus
}
#endif

#endif
