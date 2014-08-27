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

#include "filter_discordant.h"

#include "utils.h"

// #define BASIL_DEBUG

// ----------------------------------------------------------------------------
// Class DiscordantPairsFilter
// ----------------------------------------------------------------------------

void DiscordantPairsFilter::filter(std::vector<seqan::BamAlignmentRecord *> & out,
                                   std::vector<seqan::BamAlignmentRecord *> const & in)
{
    for (auto ptr : in)
    {
        seqan::BamAlignmentRecord const & record = *ptr;

        // Ignore singletons, mates mapping to different contigs, or records with too large fragment size.
        if (!hasFlagMultiple(record) || record.rID != record.rNextId || abs(record.tLen) > maxFragmentSize)
        {
#ifdef BASIL_DEBUG
            std::cerr << "DISCARDING (PRE #1)\t" << record.qName << "/" << (hasFlagLast(record) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG
            delete ptr;
            continue;
        }

        // Ignore alignments with too bad quality, except clipped ones that will become orphans later.
        if (!hasFlagUnmapped(record) && !hasClipping(record) &&
            (record.mapQ < minQuality || (minQuality > 0 && record.mapQ == 255)))
        {
#ifdef BASIL_DEBUG
            std::cerr << "DISCARDING (PRE #2)\t" << record.qName << "/" << (hasFlagLast(record) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG
            delete ptr;
            continue;
        }

        // Ignore alignments that are unmapped.
        if (hasFlagUnmapped(record) && hasFlagNextUnmapped(record))
        {
#ifdef BASIL_DEBUG
            std::cerr << "DISCARDING (PRE #3)\t" << record.qName << "/" << (hasFlagLast(record) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG
            delete ptr;
            continue;
        }

        // We now have a record that is part of a pair within the allowed insert size or a record that has an
        // anchor/shadow structure.  We will take care of missing mates in such pairs (e.g. if one aligned into a high
        // coverage region) and treat clipped records as unmapped in AntiClippingFilter.
#ifdef BASIL_DEBUG
        std::cerr << "WRITING OUT (PRE)\t" << record.qName << "/" << (hasFlagLast(record) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG
        out.push_back(ptr);
    }

}
