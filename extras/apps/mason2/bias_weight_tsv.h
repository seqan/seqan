// ==========================================================================
//                         Mason - A Read Simulator
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Author: Oliver Stolpe <oliver.stolpe@fu-berlin.de>
// ==========================================================================
// TSV format for exome bias. First column gives contig name, second column
// gives the bias weight.
// ==========================================================================

#ifndef EXTRAS_APPS_MASON2_BIAS_WEIGHT_TSV_H_
#define EXTRAS_APPS_MASON2_BIAS_WEIGHT_TSV_H_

#include <seqan/sequence.h>
#include <seqan/stream.h>

// TODO(holtgrew): Consistently use exceptions instead of return values in mason.

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct BiasWeightTsv_;
typedef seqan::Tag<BiasWeightTsv_> BiasWeightTsv;

class BiasWeightRecord
{
public:
    // The characters to insert into the genome.
    seqan::CharString contigId;

    // Size of the operation.  Always positive, also in case of deletions.
    float weight;

    BiasWeightRecord() : contigId(""), weight(1.0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TForwardIterator>
int readRecord(BiasWeightRecord & record, TForwardIterator & iter,
               BiasWeightTsv const & /*tag*/)
{
    seqan::CharString buffer, buffer2;

    // Read until tab, must not read file end.
    readUntil(buffer, iter, seqan::OrFunctor<seqan::IsTab, seqan::IsNewline>());

    // Read until tab or line break, may reach end of file.
    readUntil(buffer2, iter, seqan::IsNewline());

    record.contigId = buffer;
    record.weight = seqan::lexicalCast<float>(buffer2);

    // Make sure the weight greater than 0
    if (record.weight < 0.0)
        return 1;  // invalid weight

    // Skip rest of the line, may reach end of file.
    skipLine(iter);

    return 0;
}

#endif  // #ifndef EXTRAS_APPS_MASON2_BIAS_WEIGHT_TSV_H_
