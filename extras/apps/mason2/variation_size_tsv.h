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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// TSV format for variation sizes.  First column gives type, second column
// gives the variation size.
// ==========================================================================

#ifndef EXTRAS_APPS_MASON2_VARIATION_SIZE_TSV_H_
#define EXTRAS_APPS_MASON2_VARIATION_SIZE_TSV_H_

#include <seqan/sequence.h>
#include <seqan/stream.h>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct VariationSizeTsv_;
typedef seqan::Tag<VariationSizeTsv_> VariationSizeTsv;

class VariationSizeRecord
{
public:
    enum Kind
    {
        INVALID,
        INDEL,
        INVERSION,
        TRANSLOCATION,
        DUPLICATION
    };

    // Kind of the operation.
    Kind kind;

    // Size of the operation.  Always positive, also in case of deletions.
    int size;

    // The characters to insert into the genome.
    seqan::CharString seq;

    VariationSizeRecord() : kind(INVALID), size(-1)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TSpec>
int readRecord(VariationSizeRecord & record, seqan::RecordReader<TStream, TSpec> & reader,
               VariationSizeTsv const & /*tag*/)
{
    seqan::CharString buffer, buffer2;

    int res = 0;
    // Read until tab, must not read file end.
    if ((res = readUntilTabOrLineBreak(buffer, reader)) != 0)
        return res;

    // Read until tab, must not reach file end.
    if ((res = skipChar(reader, '\t')) != 0)
        return res;

    // Read until tab or line break, may reach end of file.
    if ((res = readUntilTabOrLineBreak(buffer2, reader)) != 0 && res != seqan::EOF_BEFORE_SUCCESS)
        return res;

    if (!lexicalCast2(record.size, buffer2))
        return 1;  // invalid number

    if (buffer == "INS" || buffer == "DEL")
        record.kind = VariationSizeRecord::INDEL;
    else if (buffer == "INV")
        record.kind = VariationSizeRecord::INVERSION;
    else if (buffer == "CTR")
        record.kind = VariationSizeRecord::TRANSLOCATION;
    else if (buffer == "DUP")
        record.kind = VariationSizeRecord::DUPLICATION;
    else
        return 1;  // Invalid record type.

    // Make sure that the record size is valid.
    if (record.size < 0)
        return 1;  // invalid size
    // Adjust record.size to type in case of indel.
    if (buffer == "DEL")
        record.size = -record.size;

    // In case of an insertion, check if we can read the sequence to be inserted.
    if (buffer == "INS" && !atEnd(reader) && value(reader) == '\t')
    {
        if ((res = skipChar(reader, '\t')) != 0)  // Skip TAB char.
            return res;

        clear(buffer2);
        if ((res = readLetters(buffer2, reader)) != 0 && res != seqan::EOF_BEFORE_SUCCESS)
            return res;
        record.seq = buffer2;
        record.size = length(record.seq);
    }

    // Skip rest of the line, may reach end of file.
    if ((res = skipLine(reader)) != 0 && res != seqan::EOF_BEFORE_SUCCESS)
        return res;

    return 0;
}

#endif  // #ifndef EXTRAS_APPS_MASON2_VARIATION_SIZE_TSV_H_
