// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_JOURNALED_STRING_TREE_TEST_CONFIG_READER_H_
#define TESTS_JOURNALED_STRING_TREE_TEST_CONFIG_READER_H_

#include <seqan/stream.h>

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct TestDeltaConfig__;
typedef Tag<TestDeltaConfig__> TestDeltaConfig_;

typedef FormattedFile<TestDeltaConfig_, Input>     TestConfigFileIn_;

template <typename TSequence>
struct TestConfigHeader_
{
    TSequence   ref;
};

template <typename TSnp, typename TDel>
struct TestConfigRecord_
{
    __int32                     pos;
    DeltaType                   deltaType;
    String<__int32>             coverage;
    TSnp                        snp;
    TDel                        del;
    String<TSnp>                ins;
    Pair<TDel, String<TSnp> >   stv;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TSpec, typename TDirection, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<TestDeltaConfig_, TDirection, TSpec>, TStorageSpec>
{
    typedef StringSet<CharString>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TSpec>
struct FileFormat<FormattedFile<TestDeltaConfig_, Input, TSpec> >
{
    typedef TestDeltaConfig_ Type;
};

template <typename T>
struct MagicHeader<TestDeltaConfig_, T>
{
    static unsigned char const VALUE[11];
};

template <typename T>
unsigned char const MagicHeader<TestDeltaConfig_, T>::VALUE[11] =
{
    '!', 'C', 'o', 'n', 'f', 'i', 'g', 'F', 'i', 'l', 'e'  // magic header
};

template <typename T>
struct FileExtensions<TestDeltaConfig_, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<TestDeltaConfig_, T>::VALUE[1] =
{
    ".txt"     // default output extension
};

// ============================================================================
// Functions
// ============================================================================

template <typename TSnp, typename TDel, typename TFwdIter, typename TContext>
inline void
readRecord(TestConfigRecord_<TSnp, TDel> & record,
           TFwdIter & fwdIter,
           TContext const & context,
           TestDeltaConfig_ const & /*tag*/)
{
//    typedef EqualsChar<'='>     TKeyValueSeparator;

    CharString buffer;
    readUntil(buffer, fwdIter, IsWhitespace());
    lexicalCast(record.pos, buffer);  // Read pos.

    skipOne(fwdIter);  // Skip tab.
    clear(buffer);
    readUntil(buffer, fwdIter, IsWhitespace());
    if (buffer == "SNP")
        record.deltaType = DELTA_TYPE_SNP;
    else if (buffer == "DEL")
        record.deltaType = DELTA_TYPE_DEL;
    else if (buffer == "INS")
        record.deltaType = DELTA_TYPE_INS;
    else if (buffer == "STV")
        record.deltaType = DELTA_TYPE_SV;
    else
        SEQAN_ASSERT_FAIL("Wrong delta type");

    skipOne(fwdIter);
    unsigned count = 0;
    while (count != length(context))
    {
        clear(buffer);
        readUntil(buffer, fwdIter, IsWhitespace());
        if (buffer == "1")
            appendValue(record.coverage, count);
        ++count;
        skipOne(fwdIter);
    }

    clear(buffer);
    if (record.deltaType == DELTA_TYPE_SNP)
    {
        readLine(buffer, fwdIter);
        record.snp = buffer[0];
    }
    else if (record.deltaType == DELTA_TYPE_DEL)
    {
        readLine(buffer, fwdIter);
        lexicalCast(record.del, buffer);
    }
    else if (record.deltaType == DELTA_TYPE_INS)
    {
        readLine(buffer, fwdIter);
        record.ins = buffer;
    }
    else
    {
        readUntil(buffer, fwdIter, EqualsChar<','>());
        lexicalCast(record.stv.i1, buffer);
        skipOne(fwdIter);
        clear(buffer);
        readLine(buffer, fwdIter);
        record.stv.i2 = buffer;
    }
}

// ----------------------------------------------------------------------------
// Function readRecord();
// ----------------------------------------------------------------------------

template <typename TSpec, typename TSnp, typename TDel>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename FormattedFile<TestDeltaConfig_, Input, TSpec>::TStream> >, void)
readRecord(TestConfigRecord_<TSnp, TDel> & record,
           FormattedFile<TestDeltaConfig_, Input, TSpec> & file)
{
    readRecord(record, file.iter, context(file), file.format);
}

// ----------------------------------------------------------------------------
// Function readHeader();
// ----------------------------------------------------------------------------

template <typename TSequence, typename TFwdIter>
inline void
readHeader(TestConfigHeader_<TSequence> & header,
           StringSet<CharString> & context,
           TFwdIter & fwdIter,
           TestDeltaConfig_ const & /*tag*/)
{
    skipLine(fwdIter);      // Skip magic header.
    while (value(fwdIter) == '#')
        skipLine(fwdIter);

    CharString buffer;
    readUntil(buffer, fwdIter, EqualsChar<'='>());
    SEQAN_ASSERT(buffer == "REF");
    skipOne(fwdIter);  // Skip separator.
    readLine(header.ref, fwdIter);

    clear(buffer);
    skipUntil(fwdIter, IsWhitespace());
    skipOne(fwdIter);  // Skip tab.
    skipUntil(fwdIter, IsWhitespace());
    skipOne(fwdIter);  // Skip tab.

    while (true)
    {
        clear(buffer);
        readUntil(buffer, fwdIter, IsWhitespace());
        if (buffer == "MOD")
        {
            skipLine(fwdIter);
            break;
        }
        appendValue(context, buffer);
        skipOne(fwdIter);  // Skip tab.
    }
}

template <typename TSequence, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename FormattedFile<TestDeltaConfig_, Input, TSpec>::TStream> >, void)
readHeader(TestConfigHeader_<TSequence> & header,
           FormattedFile<TestDeltaConfig_, Input, TSpec> & file)
{
    readHeader(header, context(file), file.iter, file.format);
}

}

#endif  // #ifndef TESTS_JOURNALED_STRING_TREE_TEST_CONFIG_READER_H_
