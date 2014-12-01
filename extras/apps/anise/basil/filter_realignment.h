// ==========================================================================
//                                 BASIL
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

#ifndef SANDBOX_HERBARIUM_APPS_BASIL_FILTER_REALIGNMENT_H_
#define SANDBOX_HERBARIUM_APPS_BASIL_FILTER_REALIGNMENT_H_

#include <memory>
#include <vector>
#include <seqan/sequence.h>

#include "library_info.h"

// ============================================================================
// Forwards
// ============================================================================

class RealignmentFilterImpl;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Typedef TBamIOContext
// ----------------------------------------------------------------------------

typedef seqan::BamIOContext<seqan::StringSet<seqan::CharString> > TBamIOContext;

// ----------------------------------------------------------------------------
// Class RealignmentFilter
// ----------------------------------------------------------------------------

class RealignmentFilter
{
public:
    RealignmentFilter(int maxFragmentSize, BamLibraryInfo::Orientation orientation,
                      seqan::CharString const & refFilename, TBamIOContext & bamIOContext,
                      int numThreads, int chunkSize);
    ~RealignmentFilter();  // for pimpl

    void filter(std::vector<seqan::BamAlignmentRecord *> & out,
                std::vector<seqan::BamAlignmentRecord *> const & in);

    void finish(std::vector<seqan::BamAlignmentRecord *> & out);

private:
    std::unique_ptr<RealignmentFilterImpl> impl;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_BASIL_FILTER_REALIGNMENT_H_
