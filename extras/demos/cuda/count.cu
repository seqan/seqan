// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_EXTRAS_DEMOS_CUDA_COUNT_H
#define SEQAN_EXTRAS_DEMOS_CUDA_COUNT_H

#include <seqan/basic_extras.h>
#include <seqan/sequence_extras.h>
#include <seqan/index_extras.h>

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class CUDAFMIndexConfig
// --------------------------------------------------------------------------
// Select the traits for the FM-index.

struct CUDAFMIndexConfig
{
    typedef TwoLevels<void>      TValuesSpec;
    typedef Naive<void>          TSentinelsSpec;

    static const unsigned SAMPLING = 8;
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function count()
// --------------------------------------------------------------------------
// Count the occurrences of a set of needles in a indexed haystack.

template <typename TIndex, typename TNeedles>
typename Size<TIndex>::Type
count(TIndex & index, TNeedles & needles)
{
    // Select the algorithm type.
    typedef Multiple<FinderSTree>                       TAlgorithmSpec;
    typedef Pattern<TNeedles, TAlgorithmSpec>           TPattern;
    typedef Finder2<TIndex, TPattern, TAlgorithmSpec>   TFinder;
    typedef OccurrencesCounter<TIndex>                  TCounter;

    // Instantiate a finder object holding the context of the search algorithm.
    TFinder finder(index);

    // Instantiate a pattern object holding the needles.
    TPattern pattern(needles);

    // Instantiate a functor object counting the number of occurrences.
    TCounter counter(pattern);

    // Find all needles in haystack and call counter() on match.
    find(finder, pattern, counter);

    // Return the number of occurrences.
    return getCount(counter);
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    // ----------------------------------------------------------------------
    // Parse input data.
    // ----------------------------------------------------------------------

    // Select the input types.
    typedef DnaString                                       THaystack;
    typedef StringSet<DnaString, Owner<ConcatDirect<> > >   TNeedles;

    if (argc < 3)
    {
        std::cerr << "USAGE: " << argv[0] << " <TEXT> <PATTERN> [<PATTERN> ...]" << std::endl;
        return 1;
    }
    
    // Create a haystack.
    THaystack haystack = argv[1];

    // Create a set of needles.
    TNeedles needles;
    for (int i = 2; i < argc; i++)
        appendValue(needles, argv[i]);

    // ----------------------------------------------------------------------
    // Build the FM-index on the CPU.
    // ----------------------------------------------------------------------

    // Select the index type.
    typedef Index<THaystack, FMIndex<void, CUDAFMIndexConfig> > TIndex;

    // Build the index over the reversed haystack.
    TIndex index(haystack);
    reverse(haystack);
    indexCreate(index);
    reverse(haystack);

    // ----------------------------------------------------------------------
    // Count on the CPU.
    // ----------------------------------------------------------------------

    omp_set_num_threads(8);
    std::cout << "CPU Occurrences: " << count(index, needles) << std::endl;

    // ----------------------------------------------------------------------
    // Copy data to the GPU.
    // ----------------------------------------------------------------------

    // Select the GPU types.
    typedef Device<TNeedles>::Type     TDeviceNeedles;
    typedef Device<TIndex>::Type       TDeviceIndex;

    // Copy the needles to the GPU.
    TDeviceNeedles deviceNeedles;
    assign(deviceNeedles, needles);

    // Copy the index to the GPU.
    TDeviceIndex deviceIndex;
    assign(deviceIndex, index);

    // ----------------------------------------------------------------------
    // Count on the GPU.
    // ----------------------------------------------------------------------

    std::cout << "GPU Occurrences: " << count(deviceIndex, deviceNeedles) << std::endl;

    return 0;
}

#endif  // SEQAN_EXTRAS_DEMOS_CUDA_COUNT_H