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
// Parsing for the "general" format of lastz.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef EXTRAS_INCLUDE_SEQAN_PARSE_LM_PARSE_LASTZ_GENERAL_H_
#define EXTRAS_INCLUDE_SEQAN_PARSE_LM_PARSE_LASTZ_GENERAL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct LastzGeneral_;
typedef Tag<LastzGeneral_> LastzGeneral;

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
..signature:readRecord(store, stream, LastzGeneral())
..summary:Read Lastz "general" format record.
..param.store:@Class.LocalMatchStore@ object to read into.
...type:Class.LocalMatchStore
..returns:$int$, 0 on success, non-0 on errors and EOF
..include:seqan/parse_lm.h
 */

template <typename TLocalMatchStore, typename TStream, typename TPassSpec>
int
readRecord(TLocalMatchStore & store,
           RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
           LastzGeneral const & /*tag*/)
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
    char subjectStrand = 'X';
    char queryStrand = 'X';
    TPosition subjectBeginPos = 0;
    TPosition subjectEndPos = 0;
    TPosition queryBeginPos = 0;
    TPosition queryEndPos = 0;

    // Field: score
    res = readDigits(buffer, recordReader);
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: name1
    readUntilChar(subjectName, recordReader, '\t');
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: strand1
    clear(buffer);
    res = readNChars(buffer, recordReader, 1);
    if (res) return res;
    subjectStrand = buffer[0];
    if (subjectStrand != '+' && subjectStrand != '-')
        return 1;  // FORMAT ERROR, should probably be a constant
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: size1
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    // NOTE: Result is ignored.
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: zstart1
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    subjectBeginPos = lexicalCast<TPosition>(buffer);
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: end1
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    subjectEndPos = lexicalCast<TPosition>(buffer);
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: name2
    readUntilChar(queryName, recordReader, '\t');
    if (res) return res;
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: strand2
    clear(buffer);
    res = readNChars(buffer, recordReader, 1);
    if (res) return res;
    queryStrand = buffer[0];
    if (queryStrand != '+' && queryStrand != '-')
        return 1;  // FORMAT ERROR, should probably be a constant
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: size2
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    // NOTE: Result is ignored.
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: zstart2
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    queryBeginPos = lexicalCast<TPosition>(buffer);
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: end2
    clear(buffer);
    res = readDigits(buffer, recordReader);
    if (res) return res;
    queryEndPos = lexicalCast<TPosition>(buffer);
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: identity
    clear(buffer);
    res = readUntilWhitespace(buffer, recordReader);
    if (res) return res;
    // NOTE: Result is ignored.
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: idPct
    clear(buffer);
    res = readUntilWhitespace(buffer, recordReader);
    if (res) return res;
    // NOTE: Result is ignored.
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: coverage
    clear(buffer);
    res = readUntilWhitespace(buffer, recordReader);
    if (res) return res;
    // NOTE: Result is ignored.
    // Skip TAB.
    res = skipChar(recordReader, '\t');
    if (res) return res;
    // Field: covPct
    clear(buffer);
    res = readUntilWhitespace(buffer, recordReader);
    if (res) return res;
    // NOTE: Result is ignored.
    // Skip to next line, ignore EOF here.
    skipLine(recordReader);

    // Finally, append the local match.
    if (subjectStrand == '-')
        ::std::swap(subjectBeginPos, subjectEndPos);
    if (queryStrand == '-')
        ::std::swap(queryBeginPos, queryEndPos);
    appendLocalMatch(store, subjectName, subjectBeginPos, subjectEndPos, queryName, queryBeginPos, queryEndPos);

    return 0;
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_PARSE_LM_PARSE_LASTZ_GENERAL_H_
