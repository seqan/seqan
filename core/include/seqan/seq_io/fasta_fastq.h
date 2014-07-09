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
// Author: David Weese <david.weese@fu-berlin.de>
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Input/Output on FASTA and FASTQ files.
// ==========================================================================

#ifndef SEQAN_SEQ_IO_FASTA_FASTQ_H_
#define SEQAN_SEQ_IO_FASTA_FASTQ_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class SequenceOutputOptions
// ----------------------------------------------------------------------------

/*!
 * @class SequenceOutputOptions
 * @headerfile <seqan/seq_io.h>
 * @brief Configuration for writing sequence (FASTA/FASTQ) files.
 * 
 * This struct is used for the configuration of writing out FASTA and FASTQ files.
 * 
 * @var int SequenceOutputOptions::lineLength;
 * @brief Length of the lines when writing out.
 * 
 * Set to <tt>-1</tt> for default behaviour (no line break for FASTQ, line length of 70 for FASTA) and <tt>0</tt> for
 * disabling line breaks.
 * 
 * @var bool SequenceOutputOptions::qualMeta;
 * @brief Whether or not to write the meta information into the <tt>"+"</tt> line before the qualities (interpreted for
 *        FASTQ only). Default is <tt>false</tt>.
 */

/**
.Class.SequenceOutputOptions
..cat:Input/Output
..summary:Configuration for writing sequence (FASTA/FASTQ) files.
..description:
This $struct$ is used for the configuration of writing out FASTA and FASTQ files.
..include:seqan/seq_io.h

.Memvar.SequenceOutputOptions#lineLength
..class:Class.SequenceOutputOptions
..type:nolink:$int$
..summary:Length of the lines when writing out.
..description:Set to $-1$ for default behaviour (no line break for FASTQ, line length of 70 for FASTA) and $0$ for disabling line breaks.

.Memvar.SequenceOutputOptions#qualMeta
..class:Class.SequenceOutputOptions
..type:nolink:$bool$
..summary:Whether or not to write the meta information into the $"+"$ line before the qualities (interpreted for FASTQ only). Default is $false$.
*/

// TODO(holtgrew): Would it be worth having two/three shortcuts for "short reads" and "genomic sequence" and faster or can the compiler optimize the creation away?

struct SequenceOutputOptions
{
    int lineLength;
    bool qualMeta;

    explicit
    SequenceOutputOptions(int lineLength = -1, bool qualMeta = false) :
        lineLength(lineLength),
        qualMeta(qualMeta)
    {}
};

// ----------------------------------------------------------------------------
// Class QualityExtractor
// ----------------------------------------------------------------------------

template <typename TValue>
struct QualityExtractor : public std::unary_function<TValue, char>
{
    inline char operator()(TValue const & x) const
    {
        return '!' + static_cast<char>(getQualityValue(x));
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function writeWrappedString()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TSequence, typename TSize>
inline void writeWrappedString(TTarget & target, TSequence const & seq, TSize lineLength)
{
    typedef typename Size<TSequence>::Type TSeqSize;
    typedef typename Iterator<TSequence const, Rooted>::Type TIter;

    TIter iter = begin(seq, Rooted());
    TSeqSize charsLeft = length(seq);
    TSeqSize charsPerLine;
    TSeqSize lineLength_ = (lineLength == 0)? maxValue<TSeqSize>() : lineLength;

    for (; charsLeft != 0; charsLeft -= charsPerLine)
    {
        charsPerLine = std::min(charsLeft, lineLength_);
        write(target, iter, charsPerLine);
        writeValue(target, '\n');
    }
}

// ----------------------------------------------------------------------------
// Function readRecord(TagSelector); Qualities inside seq
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TFwdIterator>
inline void
readRecord(TIdString & /* meta */, TSeqString & /* seq */, TFwdIterator & /* iter */,
           TagSelector<> const & /* format */)
{}

template <typename TIdString, typename TSeqString, typename TFwdIterator, typename TTagList>
inline void
readRecord(TIdString & meta, TSeqString & seq, TFwdIterator & iter, TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormatTag;

    if (value(format) == LENGTH<TTagList>::VALUE - 1)
        readRecord(meta, seq, iter, TFormatTag());
    else
        readRecord(meta, seq, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// ----------------------------------------------------------------------------
// Function readRecord(TagSelector); Qualities inside qual
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TQualString, typename TFwdIterator>
inline void
readRecord(TIdString & /* meta */, TSeqString & /* seq */, TQualString & /* qual */, TFwdIterator & /* iter */,
           TagSelector<> const & /* format */)
{}

template <typename TIdString, typename TSeqString, typename TQualString, typename TFwdIterator, typename TTagList>
inline void
readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TFwdIterator & iter, TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormatTag;

    if (value(format) == LENGTH<TTagList>::VALUE - 1)
        readRecord(meta, seq, qual, iter, TFormatTag());
    else
        readRecord(meta, seq, qual, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// ----------------------------------------------------------------------------
// Function readRecord(Raw); Qualities inside seq
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TFwdIterator>
inline void readRecord(TIdString & meta, TSeqString & seq, TFwdIterator & iter, Raw)
{
    typedef typename Value<TSeqString>::Type                            TAlphabet;
    typedef AssertFunctor<IsInAlphabet<TAlphabet>, ParseError, Fasta>   TAsserter;
    typedef OrFunctor<IsWhitespace, TAsserter>                          TIgnoreOrAssert;
    typedef EqualsChar<'>'>                                             TFastaBegin;

    clear(meta);
    clear(seq);

    skipUntil(iter, TFastaBegin());     // forward to the next '>'
    skipLine(iter);                      // assert and skip '>'
    readUntil(seq, iter, TFastaBegin(), TIgnoreOrAssert()); // read Fasta sequence
}

// ----------------------------------------------------------------------------
// Function readRecord(Raw); Qualities inside qual
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TQualString, typename TFwdIterator>
inline void readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TFwdIterator & iter, Raw const & raw)
{
    clear(qual);
    readRecord(meta, seq, iter, raw);
}

// ----------------------------------------------------------------------------
// Function readRecord(Fasta); Qualities inside seq
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TFwdIterator>
inline void readRecord(TIdString & meta, TSeqString & seq, TFwdIterator & iter, Fasta)
{
    typedef typename Value<TSeqString>::Type                            TAlphabet;
    typedef AssertFunctor<IsInAlphabet<TAlphabet>, ParseError, Fasta>   TAsserter;
    typedef OrFunctor<IsWhitespace, TAsserter>                          TIgnoreOrAssert;
    typedef EqualsChar<'>'>                                             TFastaBegin;

    clear(meta);
    clear(seq);

    skipUntil(iter, TFastaBegin());     // forward to the next '>'
    skipOne(iter);                      // assert and skip '>'

    readLine(meta, iter);               // read Fasta id
    readUntil(seq, iter, TFastaBegin(), TIgnoreOrAssert()); // read Fasta sequence
}

// ----------------------------------------------------------------------------
// Function readRecord(Fasta); Qualities inside qual
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TQualString, typename TFwdIterator>
inline void readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TFwdIterator & iter, Fasta)
{
    clear(qual);
    readRecord(meta, seq, iter, Fasta());
}

// ----------------------------------------------------------------------------
// Function readRecord(Fastq); Qualities inside seq
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TFwdIterator>
inline void readRecord(TIdString & meta, TSeqString & seq, TFwdIterator & iter, Fastq)
{
    typedef typename Value<TSeqString>::Type                                TSeqAlphabet;
    typedef AssertFunctor<IsInAlphabet<TSeqAlphabet>, ParseError, Fastq>    TSeqAsserter;
    typedef OrFunctor<IsWhitespace, TSeqAsserter>                           TSeqIgnoreOrAssert;
    typedef EqualsChar<'@'>                                                 TFastqBegin;
    typedef EqualsChar<'+'>                                                 TQualsBegin;

    clear(meta);
    clear(seq);

    skipUntil(iter, TFastqBegin());     // forward to the next '@'
    skipOne(iter);                      // skip '@'

    readLine(meta, iter);               // read Fastq id

    readUntil(seq, iter, TQualsBegin(), TSeqIgnoreOrAssert());  // read Fastq sequence
    skipOne(iter, TQualsBegin());       // assert and skip '+'
    skipLine(iter);                     // skip optional 2nd Fastq id

    skipUntil(iter, TFastqBegin());     // forward to the next '@'
}

// ----------------------------------------------------------------------------
// Function readRecord(Fastq); Qualities inside qual
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TQualString, typename TFwdIterator>
inline void readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TFwdIterator & iter, Fastq)
{
    typedef typename Value<TSeqString>::Type                                TSeqAlphabet;
    typedef typename Value<TQualString>::Type                               TQualAlphabet;
    typedef AssertFunctor<IsInAlphabet<TSeqAlphabet>, ParseError, Fastq>    TSeqAsserter;
    typedef AssertFunctor<IsInAlphabet<TQualAlphabet>, ParseError, Fastq>   TQualAsserter;
    typedef OrFunctor<IsWhitespace, TSeqAsserter>                           TSeqIgnoreOrAssert;
    typedef OrFunctor<IsBlank, TQualAsserter>                               TQualIgnoreOrAssert;
    typedef EqualsChar<'@'>                                                 TFastqBegin;
    typedef EqualsChar<'+'>                                                 TQualsBegin;

    clear(meta);
    clear(seq);
    clear(qual);

    skipUntil(iter, TFastqBegin());     // forward to the next '@'
    skipOne(iter);                      // skip '@'

    readLine(meta, iter);               // read Fastq id

    readUntil(seq, iter, TQualsBegin(), TSeqIgnoreOrAssert());  // read Fastq sequence
    skipOne(iter, TQualsBegin());       // assert and skip '+'
    skipLine(iter);                     // skip optional 2nd Fastq id

    readUntil(qual, iter, IsNewline(), TQualIgnoreOrAssert());  // read Fastq qualities
    skipUntil(iter, TFastqBegin());     // forward to the next '@'
}

// ----------------------------------------------------------------------------
// Function writeRecord(Fasta); Qualities inside seq
// ----------------------------------------------------------------------------

/*!
 * @fn FastaFastqIO#writeRecord
 * @headerfile <seqan/seq_io.h>
 * @brief Write one FASTA or FASTQ record.
 * 
 * @signature int writeRecord(target, id, seq, tag[, options]);
 * @signature int writeRecord(target, id, seq, quals, tag[, options]);
 * 
 * @param[in,out] target  The target to write to.  Type: StreamConcept
 * @param[in]     id      ID/Meta information line to write out. Types: SequenceConcept
 * @param[in]     seq     Sequence to write out.  Type: SequenceConcept
 * @param[in]     quals   ASCII quality characters to write out.  Types: SequenceConcept
 * @param[in]     tag     The format selector. Types: nolink:<tt>Fasta</tt>, <tt>Fastq</tt>
 * @param[in]     options if not supplied, defaults are chosen.  Types: SequenceOutputOptions
 *
 * @return int 0 on success, non-0 value on errors.
 */

/**
.Function.FASTA/FASTQ I/O#writeRecord
..summary:Write one FASTA or FASTQ record.
..signature:int writeRecord(target, id, seq, tag[, options])
..signature:int writeRecord(target, id, seq, quals, tag[, options])
..param.target:The target to write to.
...type:Concept.StreamConcept
..param.id:ID/Meta information line to write out.
...type:Concept.SequenceConcept
..param.seq:Sequence to write out.
...type:Concept.SequenceConcept
..param.quals:ASCII quality characters to write out.
...type:Concept.SequenceConcept
..param.tag:The format selector.
...type:nolink:$Fasta$, $Fastq$
..param.options:if not supplied defaults are chosen.
...type:Class.SequenceOutputOptions
..include:seqan/seq_io.h
*/


template <typename TTarget, typename TIdString, typename TSeqString>
inline void writeRecord(TTarget & target,
                        TIdString const & meta,
                        TSeqString const & seq,
                        Fasta const & /*tag*/,
                        SequenceOutputOptions const & options = SequenceOutputOptions())
{
    writeValue(target, '>');
    write(target, meta);
    writeValue(target, '\n');

    writeWrappedString(target, seq, (options.lineLength < 0)? 70 : options.lineLength); // 70bp wrapping, by default
}

// ----------------------------------------------------------------------------
// Function writeRecord(Fastq); Qualities inside qual
// ----------------------------------------------------------------------------

template <typename TTarget, typename TIdString, typename TSeqString, typename TQualString>
inline void writeRecord(TTarget & target,
                        TIdString const & meta,
                        TSeqString const & seq,
                        TQualString const & qual,
                        Fastq const & /*tag*/,
                        SequenceOutputOptions const & options = SequenceOutputOptions())
{
    writeValue(target, '@');
    write(target, meta);
    writeValue(target, '\n');

    int lineLength = (options.lineLength < 0)? 0 : options.lineLength;  // no wrapping, by default
    writeWrappedString(target, seq, lineLength);

    writeValue(target, '+');
    if (options.qualMeta)
        write(target, meta);
    writeValue(target, '\n');

    writeWrappedString(target, qual, lineLength);
}

// ----------------------------------------------------------------------------
// Function writeRecord(Fastq); Qualities inside seq
// ----------------------------------------------------------------------------

template <typename TTarget, typename TIdString, typename TSeqString>
inline void writeRecord(TTarget & target,
                        TIdString const & meta,
                        TSeqString const & seq,
                        Fastq const & tag,
                        SequenceOutputOptions const & options = SequenceOutputOptions())
{
    typedef QualityExtractor<typename Value<TSeqString>::Type> TQualityExtractor;
    ModifiedString<TSeqString const, ModView<TQualityExtractor> > quals(seq);
    writeRecord(target, meta, seq, quals, tag, options);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQ_IO_FASTA_FASTQ_H_
