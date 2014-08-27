// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_READ_MAPPING_IMPL_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_READ_MAPPING_IMPL_H_

#include <memory>
#include <vector>

#include <seqan/sequence.h>

// ============================================================================
// Forwards
// ============================================================================

class FilteringBestMapperImpl;
class FilteringAllBestMapperImpl;
class FilteringAllMapperImpl;
namespace seqan { class BamAlignmentRecord; }

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class FilteringAllBestMapper
// ----------------------------------------------------------------------------

// Returns all matches for a read.
//
// Match management is limited to creating at most one match for each contig.

class FilteringAllMapper
{
public:
    struct Options
    {
        double errorRate { 0.03 };
        unsigned maxMatches { 100 };
        bool debug { false };

        Options() = default;
        Options(double errorRate, unsigned maxMatches, bool debug = false) :
                errorRate(errorRate), maxMatches(maxMatches), debug(debug)
        {}
    };

    FilteringAllMapper();
    explicit FilteringAllMapper(Options options);
    ~FilteringAllMapper();  // for pimpl

    // Perform read mapping.
    //
    // Note that _qID is properly set in result.
    //
    // contigs is non-const since we need to reverse-complement it, reads and quals are non-const because of
    // Index/Holder const issues.
    //
    // The records in result are sorted by read ID.
    void run(std::vector<seqan::BamAlignmentRecord> & result,
             seqan::StringSet<seqan::Dna5String> & contigs,  // non-const, will be reverse-complemented
             seqan::StringSet<seqan::Dna5String> /*const*/ & reads,
             seqan::StringSet<seqan::CharString> /*const*/ & quals,
             seqan::StringSet<seqan::CharString> const & readNames);

private:

    std::unique_ptr<FilteringAllMapperImpl> impl;
};

// ----------------------------------------------------------------------------
// Class FilteringAllBestMapper
// ----------------------------------------------------------------------------

// Returns all best matches for a read.
//
// Match management is slightly more complex than for FilteringBestMapper.

class FilteringAllBestMapper
{
public:
    struct Options
    {
        double errorRate { 0.03 };
        unsigned maxMatches { 100 };

        Options() = default;
        Options(double errorRate, unsigned maxMatches) :
                errorRate(errorRate), maxMatches(maxMatches)
        {}
    };

    FilteringAllBestMapper();
    explicit FilteringAllBestMapper(Options options);
    ~FilteringAllBestMapper();  // for pimpl

    // Perform read mapping.
    //
    // Note that _qID is properly set in result.
    //
    // contigs is non-const since we need to reverse-complement it, reads and quals are non-const because of
    // Index/Holder const issues.
    //
    // The records in result are sorted by read ID.
    void run(std::vector<seqan::BamAlignmentRecord> & result,
             seqan::StringSet<seqan::Dna5String> & contigs,  // non-const, will be reverse-complemented
             seqan::StringSet<seqan::Dna5String> /*const*/ & reads,
             seqan::StringSet<seqan::CharString> /*const*/ & quals,
             seqan::StringSet<seqan::CharString> const & readNames);

private:

    std::unique_ptr<FilteringAllBestMapperImpl> impl;
};

// ----------------------------------------------------------------------------
// Class FilteringBestMapper
// ----------------------------------------------------------------------------

// Maps reads against a contig in a best-mapper fashion using a pigeonhole filter and banded Myers verification.
//
// The implementation of match management is greatly simplified since the result contains one BAM record for each input
// read.

class FilteringBestMapper
{
public:
    FilteringBestMapper(double errorRate);
    ~FilteringBestMapper();  // for pimpl

    // Perform read mapping.
    //
    // contigs is non-const since we need to reverse-complement it, reads and quals are non-const because of
    // Index/Holder const issues.
    void run(std::vector<seqan::BamAlignmentRecord> & result,
             seqan::StringSet<seqan::Dna5String> & contigs,  // non-const, will be reverse-complemented
             seqan::StringSet<seqan::Dna5String> /*const*/ & reads,
             seqan::StringSet<seqan::CharString> /*const*/ & quals,
             seqan::StringSet<seqan::CharString> const & readNames);
             

private:
    std::unique_ptr<FilteringBestMapperImpl> impl;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_READ_MAPPING_IMPL_H_
