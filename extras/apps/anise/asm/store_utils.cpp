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

#include "store_utils.h"

#include <type_traits>

#include <seqan/realign.h>

// ----------------------------------------------------------------------------
// Function printStore()
// ----------------------------------------------------------------------------

void printStore(std::ostream & out, TFragmentStore const & storeC, unsigned contigID)
{
    TFragmentStore store(storeC);

    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    for (unsigned i = 0; i < length(layout.contigRows); ++i)
    {
        if (contigID != (unsigned)-1 && contigID != i)
            continue;
        // Get coordinates to plot.
        __int64 l = 0;
        __int64 r = l;
        for (unsigned j = 0; j < length(layout.contigRows[i]); ++j)
        {
            unsigned id = back(layout.contigRows[i][j]);
            if (r < store.alignedReadStore[id].beginPos)
                r = store.alignedReadStore[id].beginPos;
            if (r < store.alignedReadStore[id].endPos)
                r = store.alignedReadStore[id].endPos;
        }

        out << ">multi-read-alignment: contig_" << i << "\n";
        printAlignment(out, seqan::Raw(), layout, store, i, l, r, 0, 1000);
    }
}

// ----------------------------------------------------------------------------
// Function permuteStoreContigs()
// ----------------------------------------------------------------------------

void permuteStoreContigs(TFragmentStore & store,
                         std::vector<unsigned> const & permutation)
{
    using std::swap;

    // Update contigStore.
    decltype(store.contigStore) tmp;
    resize(tmp, length(store.contigStore));
    for (unsigned i = 0; i < length(store.contigStore); ++i)
        swap(tmp[permutation[i]], store.contigStore[i]);
    swap(tmp, store.contigStore);

    // Update contig IDs.
    for (auto & el : store.alignedReadStore)
        el.contigId = permutation[el.contigId];
}
