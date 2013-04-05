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

#ifndef EXTRAS_INCLUDE_SEQAN_PARSE_LM_PARSE_STELLAR_GFF_H_
#define EXTRAS_INCLUDE_SEQAN_PARSE_LM_PARSE_STELLAR_GFF_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct StellarGff_;
typedef Tag<StellarGff_> StellarGff;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

/**
.Function.readRecord
..class:Class.LocalMatchStore
..cat:Local Match Store
..signature:readRecord(store, stream, StellarGff())
..param.store:@Class.LocalMatchStore@ object to read into.
...type:Class.LocalMatchStore
..returns:$int$, 0 on success, non-0 on errors and EOF
..include:seqan/parse_lm.h
 */

template <typename TLocalMatchStore, typename TStream, typename TPassSpec>
int
readRecord(TLocalMatchStore & store,
           RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
           StellarGff const & /*tag*/)
{
    typedef typename TLocalMatchStore::TPosition TPosition;
    //typedef typename TLocalMatchStore::TPosition TId;
    
    // Read line.
    CharString buffer;
    int res = 0;
    char subjectStrand = 'X';

    CharString subjectName;
    CharString queryName;
    TPosition subjectBeginPos = 0;
    TPosition subjectEndPos = 0;
    TPosition queryBeginPos = 0;
    TPosition queryEndPos = 0;

    // Field: SUBJECT
    res = readUntilChar(subjectName, recordReader, '\t');
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: SOURCE
    res = skipUntilChar(recordReader, '\t');
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: TYPE
    res = skipUntilChar(recordReader, '\t');
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: START
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    subjectBeginPos = lexicalCast<TPosition>(buffer) - 1;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: END
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    subjectEndPos = lexicalCast<TPosition>(buffer);
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: SCORE
    res = skipUntilChar(recordReader, '\t');
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: STRAND
    clear(buffer);
    res = readNChars(buffer, recordReader, 1);
    if (res) return res;
    subjectStrand = buffer[0];
    if (subjectStrand != '+' && subjectStrand != '-')
        return 1;  // FORMAT ERROR, should probably be a constant
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: FRAME
    res = skipNChars(recordReader, 1);
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;

    // The GROUP field contains the information about the query.
    // query sequence
    res = readUntilChar(queryName, recordReader, ';');
    if (res) return res;
    // Skip semicolon.
    res = skipChar(recordReader, ';');
    if (res) return res;
    // "seq2Range="
    clear(buffer);
    res = readUntilChar(buffer, recordReader, '=');
    if (res) return res;
    if (buffer != "seq2Range")
        return 1;  // FORMAT ERROR, should probably be a constant
    // Skip '='
    res = skipChar(recordReader, '=');
    if (res) return res;
    // query begin pos
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    queryBeginPos = lexicalCast<TPosition>(buffer) - 1;
    // skip comma
    res = skipChar(recordReader, ',');
    if (res) return res;
    // query end pos
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    queryEndPos = lexicalCast<TPosition>(buffer);
    // Skip semicolon.
    res = skipChar(recordReader, ';');
    if (res) return res;
    // Skip "eValue=*;"
    clear(buffer);
    res = readUntilChar(buffer, recordReader, '=');
    if (res) return res;
    if (buffer != "eValue")
        return 1;  // FORMAT ERROR, should probably be a constant
    // Skip '='
    res = skipChar(recordReader, '=');
    if (res) return res;
    // Skip until semicolon.
    res = skipUntilChar(recordReader, ';');
    if (res) return res;
    // And skip the semicolon.
    res = skipChar(recordReader, ';');
    if (res) return res;
    // Skip "cigar".
    clear(buffer);
    res = readUntilChar(buffer, recordReader, '=');
    if (res) return res;
    if (buffer != "cigar")
        return 1;  // FORMAT ERROR, should probably be a constant
    // Skip '='
    res = skipChar(recordReader, '=');
    if (res) return res;
    // Read cigar string.
    clear(buffer);
    res = readUntilChar(buffer, recordReader, ';');
    // Skip semicolon.
    res = skipChar(recordReader, ';');
    if (res) return res;
    // ignore rest of the field, skip to next line
    skipLine(recordReader);

	// Finally, append the local match.
    if (subjectStrand == '-')
        ::std::swap(subjectBeginPos, subjectEndPos);
    appendLocalMatch(store, subjectName, subjectBeginPos, subjectEndPos, queryName, queryBeginPos, queryEndPos, buffer);

    return 0;
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_PARSE_LM_PARSE_STELLAR_GFF_H_
