// ==========================================================================
//                                  parse_lm
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

// SEQAN_NO_GENERATED_FORWARDS

#ifndef EXTRAS_INCLUDE_SEQAN_PARSE_LM_PARSE_BLASTN_TABULAR_H_
#define EXTRAS_INCLUDE_SEQAN_PARSE_LM_PARSE_BLASTN_TABULAR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct BlastnTabular_;
typedef Tag<BlastnTabular_> BlastnTabular;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

/**
.Function.LocalMatchStore#readRecord
..class:Class.LocalMatchStore
..cat:Local Match Store
..signature:readRecord(store, stream, BlastnTabular())
..param.store:@Class.LocalMatchStore@ object to read into.
...type:Class.LocalMatchStore
..returns:$int$, 0 on success, non-0 on errors and EOF
..include:seqan/parse_lm.h
 */

template <typename TLocalMatchStore, typename TStream, typename TPassSpec>
int
readRecord(TLocalMatchStore & store,
           RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
           BlastnTabular const & /*tag*/)
{
    typedef typename TLocalMatchStore::TPosition TPosition;
    //typedef typename TLocalMatchStore::TPosition TId;
    
    if (atEnd(recordReader))
        return 1;
    // Skip any comments.
    while (value(recordReader) == '#') {
        int res = skipLine(recordReader);
        if (res)
            return res;
        if (atEnd(recordReader))
            return 1;
    }

    SEQAN_ASSERT_NEQ(value(recordReader), '#');

    // Read line.
    CharString buffer;
    int res = 0;

    CharString subjectName;
    CharString queryName;
    TPosition subjectBeginPos = 0;
    TPosition subjectEndPos = 0;
    TPosition queryBeginPos = 0;
    TPosition queryEndPos = 0;

    // Field: query id
    res = readUntilChar(queryName, recordReader, '\t');
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: subject id
    res = readUntilChar(subjectName, recordReader, '\t');
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: % identity
    res = skipUntilChar(recordReader, '\t');
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: alignment length
    res = skipUntilChar(recordReader, '\t');
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: mismatches
    res = skipUntilChar(recordReader, '\t');
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: gap opens
    res = skipUntilChar(recordReader, '\t');
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: q. start
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    queryBeginPos = lexicalCast<TPosition>(buffer) - 1;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: q. end
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    queryEndPos = lexicalCast<TPosition>(buffer);
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: s. start
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    subjectBeginPos = lexicalCast<TPosition>(buffer) - 1;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: s. end
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    subjectEndPos = lexicalCast<TPosition>(buffer);
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: evalue
    res = skipUntilChar(recordReader, '\t');
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: bit score
    res = skipUntilWhitespace(recordReader);
    if (res) return res;

    // Finally, append the local match.
    appendLocalMatch(store, subjectName, subjectBeginPos, subjectEndPos, queryName, queryBeginPos, queryEndPos);

    return 0;
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_PARSE_LM_PARSE_BLASTN_TABULAR_H_
