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

#include "realignment.h"

#include <seqan/realign.h>

void realign(TFragmentStore & store, TWindowBorders const & windowBorders, int bandwidth, int environment, bool debug)
{
    seqan::RealignmentOptions_ roptions;
    roptions.method = seqan::RealignmentOptions_::ANSON_MYERS_GOTOH;
    roptions.bandwidth = bandwidth;
    roptions.environment = environment;
    roptions.includeReference = false;
    roptions.debug = debug;
    roptions.printTiming = debug;

    // TODO(holtgrew): use make_adjacent_iterators() here.
    for (unsigned i = 0; i < length(store.contigStore); ++i)
        for (unsigned j = 0; j + 1 < length(windowBorders[i]); ++j)
        {
            unsigned windowBegin = windowBorders[i][j], windowEnd = windowBorders[i][j + 1];
            seqan::AnsonMyersRealigner_<TFragmentStore> realigner(store, roptions);
            realigner.run(i, windowBegin, windowEnd);
    }
}

// Realign the store, this will also recreate the contigs.
void realign(TFragmentStore & store, std::vector<unsigned> const & contigIDs,
             int bandwidth, int environment, bool debug)
{
    seqan::RealignmentOptions_ roptions;
    roptions.method = seqan::RealignmentOptions_::ANSON_MYERS_GOTOH;
    roptions.bandwidth = bandwidth;
    roptions.environment = environment;
    roptions.includeReference = false;
    roptions.debug = debug;
    roptions.printTiming = debug;

    for (auto contigID : contigIDs)
    {
        seqan::AnsonMyersRealigner_<TFragmentStore> realigner(store, roptions);
        realigner.run(contigID);
    }
}

void realign(TFragmentStore & store, int bandwidth, int environment, bool debug)
{
    std::vector<unsigned> contigIDs;
    for (unsigned i = 0; i < length(store.contigStore); ++i)
        contigIDs.push_back(i);
    realign(store, contigIDs, bandwidth, environment, debug);
}
