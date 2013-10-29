// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Read FASTA and FASTQ files.
// ==========================================================================

#ifndef SEQAN_SEQ_IO_READ_FASTA_FASTQ_H_
#define SEQAN_SEQ_IO_READ_FASTA_FASTQ_H_

namespace seqan {

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord(Fasta)
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TFwdIterator>
inline void
readRecord(TIdString & meta, TSeqString & seq, TFwdIterator & iter, Fasta)
{
    EqualsChar<'>'> fastaBegin;
    IgnoreOrAssertFunctor<IsWhitespace, IsInAlphabet<typename Value<TSeqString>::Type>, std::runtime_error>
        ignoreWhiteSpaceAndAssertAlphabet("Invalid character in Fasta sequence!");

    clear(meta);
    clear(seq);

    skipUntil(iter, fastaBegin);    // forward to the next '>'
    ++iter;                         // skip '>'
    readLine(meta, iter);           // read Fasta id
    readUntil(seq, iter, fastaBegin, ignoreWhiteSpaceAndAssertAlphabet);    // read Fasta sequence
}

// ----------------------------------------------------------------------------
// Function readRecord(Fastq)
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TQualString, typename TFwdIterator>
inline void
readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TFwdIterator & iter, Fastq)
{
    EqualsChar<'@'> fastqBegin;
    EqualsChar<'+'> qualsBegin;

    IgnoreOrAssertFunctor<IsWhitespace, IsInAlphabet<typename Value<TSeqString>::Type>, std::runtime_error>
        ignoreWhiteSpaceAndAssertAlphabet("Invalid sequence character in Fastq sequence!");

    IgnoreOrAssertFunctor<IsBlank, IsInAlphabet<typename Value<TQualString>::Type>, std::runtime_error>
        ignoreBlankAssertQuality("Invalid quality character in Fastq sequence!");

    clear(meta);
    clear(seq);
    clear(qual);

    skipUntil(iter, fastqBegin);    // forward to the next '@'
    ++iter;                         // skip '@'

    readLine(meta, iter);           // read Fastq id

    readUntil(seq, iter, qualsBegin, ignoreWhiteSpaceAndAssertAlphabet);    // read Fastq sequence
    skipLine(iter);                 // skip '+' and 2nd Fastq id

    readUntil(qual, iter, IsNewline(), ignoreBlankAssertQuality);           // read Fastq qualities
    skipUntil(iter, fastqBegin);    // forward to the next '@'
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQ_IO_READ_FASTA_FASTQ_H_
