// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_ANISE_STORE_UTILS_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_ANISE_STORE_UTILS_H_

#include <memory>
#include <vector>

#include <seqan/sequence.h>

#include "asm/frag_store.h"

// ============================================================================
// Forwards
// ============================================================================

class DifferingReadsRemoverImpl;
class ContigSplitterImpl;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Typedef TProfileChar
// --------------------------------------------------------------------------

typedef seqan::ProfileChar<seqan::Dna5> TProfileChar;

// --------------------------------------------------------------------------
// Class DifferingReadsRemover
// --------------------------------------------------------------------------

// Removes reads from fragment store that have a too high error rate compared to the contig.
//
// These reads are moved to their own contig (with contigStore and contigNameStore resized appropriately.  Note that we
// we also remove gaps and shift alignments appropriately.
//
// Also splits at zero coverage regions.

class DifferingReadsRemover
{
public:
    DifferingReadsRemover(double maxErrorRate);
    ~DifferingReadsRemover();  // for pimpl

    void run(TFragmentStore & store, std::vector<unsigned> & readIDs) const;

private:
    std::unique_ptr<DifferingReadsRemoverImpl> impl;
};

// --------------------------------------------------------------------------
// Class ContigSplitter
// --------------------------------------------------------------------------

// Split contigs of fragment store at zero-coverage positions and at positions where a split would improve the pairwise
// alignment score (given a scheme with match=1, mismatch=-2, gaps=-2).

class ContigSplitter
{
public:
    ContigSplitter(TFragmentStore & store, bool debug = false);
    ~ContigSplitter();  // for pimpl

    // Perform the separation
    void run();

private:

    std::unique_ptr<ContigSplitterImpl> impl;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function countErrors()
// ----------------------------------------------------------------------------

// Number of errors from read alignment vs. contig.

int countErrors(TFragmentStore /*const*/ & store, TAlignedReadStoreElement /*const*/ & el);

// ----------------------------------------------------------------------------
// Function printStore()
// ----------------------------------------------------------------------------

// Print fragment store to output.

void printStore(std::ostream & out, TFragmentStore const & store);

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_ANISE_STORE_UTILS_H_
