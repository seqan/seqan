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

#include "frag_store_builder.h"

// --------------------------------------------------------------------------
// Class FragmentStoreBuilderImpl
// --------------------------------------------------------------------------

class FragmentStoreBuilderImpl
{
    // The records to build from.  Must be pairs with the left mate even and right mate in uneven positions.
    std::vector<seqan::BamAlignmentRecord> const & records;
    // The reference sequences to use.
    seqan::StringSet<seqan::Dna5String> const & refs;

public:

    FragmentStoreBuilderImpl(std::vector<seqan::BamAlignmentRecord> const & records,
                             seqan::StringSet<seqan::Dna5String> const & refs) :
        records(records), refs(refs)
    {}

    void run(TFragmentStore & store) const
    {
        SEQAN_ASSERT_EQ(records.size() % 2u, 0u);

        fillContigs(store);
        for (unsigned i = 0; i < records.size(); i += 2)
            appendAlignmentPair(store, records[i], records[i + 1], i / 2);

        // Convert the matches to global alignments.
        convertMatchesToGlobalAlignment(store, seqan::Score<int, seqan::Simple>(1, -1, -1), seqan::False());
    }

private:

    // Fill store.contigStore and store.contigNameStore.
    void fillContigs(TFragmentStore & store) const
    {
        resize(store.contigStore, length(refs));
        for (unsigned i = 0; i < length(store.contigStore); ++i)
            store.contigStore[i].seq = refs[i];

        for (unsigned i = 0; i < length(refs); ++i)  // fill contig names
        {
            std::stringstream ss;
            ss << "contig_" << i;
            appendValue(store.contigNameStore, ss.str());
        }
    }

    // Add an alignment pair.
    void appendAlignmentPair(TFragmentStore & store,
                             seqan::BamAlignmentRecord const & left,
                             seqan::BamAlignmentRecord const & right,
                             unsigned pairMatchID) const
    {
        SEQAN_CHECK(hasFlagFirst(left) || hasFlagFirst(right), "One must be the first.");
        bool flipped = !hasFlagFirst(left);  // need to switch first/last roles

        // Append the mate pair.
        auto matePairID = appendMatePair(store, left.seq, right.seq);
        auto leftID = store.matePairStore[matePairID].readId[0];
        auto rightID = store.matePairStore[matePairID].readId[1];
        // The first/last information is stored in store.matePairStore and we need to switch these roles.
        if (flipped)
            std::swap(store.matePairStore[matePairID].readId[0], store.matePairStore[matePairID].readId[1]);

        // std::cerr << "APPENDING PAIR\tleftID=" << leftID << "\trightID=" << rightID << "\n"
        //           << "  left\t" << left.qName << "\n";

        // Append left alignment.
        if (!hasFlagUnmapped(left))
        {
            // std::cerr << "    NOT UNMAPPED\n";
            appendAlignedRead(store, leftID, left.rID, (int)left.beginPos,
                              (int)left.beginPos + (int)getAlignmentLengthInRef(left), pairMatchID);
            if (hasFlagRC(left))
            {
                std::swap(back(store.alignedReadStore).beginPos, back(store.alignedReadStore).endPos);
                reverseComplement(store.readSeqStore[leftID]);
            }
            resize(store.alignQualityStore, length(store.alignedReadStore));
        }
        // std::cerr << "  right\t" << right.qName << "\tpairMatchID=" << pairMatchID << "\n";
        // Append right alignment.
        if (!hasFlagUnmapped(right))
        {
            // std::cerr << "    NOT UNMAPPED\n";
            appendAlignedRead(store, rightID, right.rID, (int)right.beginPos,
                              (int)right.beginPos + (int)getAlignmentLengthInRef(right), pairMatchID);
            if (hasFlagRC(right))
            {
                std::swap(back(store.alignedReadStore).beginPos, back(store.alignedReadStore).endPos);
                reverseComplement(store.readSeqStore[rightID]);
            }
            resize(store.alignQualityStore, length(store.alignedReadStore));
        }
    }

};

// --------------------------------------------------------------------------
// Class FragmentStoreBuilder
// --------------------------------------------------------------------------

FragmentStoreBuilder::FragmentStoreBuilder(std::vector<seqan::BamAlignmentRecord> const & records,
                                           seqan::StringSet<seqan::Dna5String> const & refs) :
        impl(new FragmentStoreBuilderImpl(records, refs))
{}

FragmentStoreBuilder::~FragmentStoreBuilder()
{}

void FragmentStoreBuilder::run(TFragmentStore & store) const
{
    impl->run(store);
}
