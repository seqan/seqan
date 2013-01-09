// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Code for writing FASTA and FASTQ files
// ==========================================================================

#ifndef SEQAN_SEQ_IO_WRITE_H_
#define SEQAN_SEQ_IO_WRITE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Enum.FastAQOutputOptions
..cat:Input / Output
..summary:Enum with flags for writing Fasta and Fastq files
..value.NO_FLAGS:No Flags are set.
..value.LINEBREAKS:Introduce linebreaks in sequences at column 70
..value.WRITEQUALITIESMETA:Write the meta information into the "+" line in addition to writing it to the "\at" line [FASTQ only]
..value.DEFAULT_FASTA: currently == LINEBREAKS
..value.DEFAULT_FASTQ: currently == NO_FLAGS
..include:seqan/stream.h
 */

// TODO(holtgrew): Rename enums?

enum FastAQOutputOptions
{
    NO_FLAGS = 0,
    LINEBREAKS = 1,           // introduce linebreaks into sequence at col 70
    WRITEQUALITIESMETA = 2,    // for FastQ write meta line before qualities too
//  FOOFOO = 4,
//  BARBAR = 8
    DEFAULT_FASTA = 1,
    DEFAULT_FASTQ = 0
};

// ============================================================================
// Metafunctions
// ============================================================================


// ============================================================================
// Functions
// ============================================================================


// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------


/**
.Function.writeRecord
..signature:writeRecord(TStream & stream, TIdString const & meta, TSeqString const & seq, Fasta const &[, FastAQOutputOptions const options])
..param.options:if not supplied defaults are chosen, see @Enum.FastAQOutputOptions@
...type:Enum.FastAQOutputOptions
*/
template <typename TStream,
          typename TIdString,
          typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            Fasta const & /*tag*/,
            unsigned const options)
{
    int res = streamPut(stream, '>');
    if (res)
        return res;

    res = streamPut(stream, meta);
    if (res)
        return res;

    res = streamPut(stream, '\n');
    if (res)
        return res;

    if (options & LINEBREAKS)
    {
        // write stream character by character
        typename Iterator<TSeqString const, Standard>::Type it = begin(seq);
        typename Iterator<TSeqString const, Standard>::Type it_end = end(seq);
        for (unsigned long l = 0; it < it_end; ++it)
        {
            res = streamPut(stream, *it);
            if (res)
                return res;
            if (++l == 70)
            {
                res = streamPut(stream, '\n');
                l = 0;
                if (res)
                    return res;
            }
        }
        if (res)
            return res;
    } else
    {
        res = streamPut(stream, seq);
        if (res)
            return res;
    }

    res = streamPut(stream, '\n');
    return res;
}

// FASTA
template <typename TStream,
          typename TIdString,
          typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            Fasta const & /*tag*/)
{
    return writeRecord(stream, meta, seq, Fasta(), DEFAULT_FASTA);
}


/**
.Function.writeRecord
..signature:writeRecord(TStream & stream, TIdString const & meta, TSeqString const & seq, TQualString const & qual, Fastq const &[, FastAQOutputOptions const options])
*/

template <typename TStream, typename TSequence, typename TQualString>
inline int _writeRecordFastq(TStream & stream, TSequence const & seq, TQualString const &, unsigned const options, True const & /*HasQualities<Value<TSequence>::Type>::VALUE*/)
{
    int res = 0;

    if (options & LINEBREAKS)
    {
        // write stream character by character
        typename Iterator<TSequence const>::Type it = begin(seq);
        typename Iterator<TSequence const>::Type it_end = end(seq);
        for (unsigned long l = 0; it < it_end; ++it)
        {
            res = streamPut(stream, static_cast<char>('!' + getQualityValue(*it)));
            if (res)
                return res;
            if (++l == 70)
            {
                res = streamPut(stream, '\n');
                l = 0;
                if (res)
                    return res;
            }
        }
        if (res)
            return res;
    }
    else
    {
        typename Iterator<TSequence const>::Type it = begin(seq);
        typename Iterator<TSequence const>::Type it_end = end(seq);
        for (; it < it_end; ++it)
        {
            res = streamPut(stream, static_cast<char>('!' + getQualityValue(*it)));
            if (res)
                return res;
        }
    }

    return 0;
}


template <typename TStream, typename TSequence, typename TQualString>
inline int _writeRecordFastq(TStream & stream, TSequence const & seq, TQualString const & qual, unsigned const options, False const & /*HasQualities<Value<TSequence>::Type>::VALUE*/)
{
    int res = 0;

    if (qual == "") // we don't actually have qualities
    {
        if (options & LINEBREAKS)
        {
            for (unsigned long i = 0, l = 0;
                 i < length(seq);
                 ++i)
            {
                res = streamPut(stream, char(126));
                if (res)
                    return res;
                if (++l == 70)
                {
                    res = streamPut(stream, '\n');
                    l = 0;
                    if (res)
                        return res;
                }
            }
        } else
        {
            for (unsigned long i = 0; i < length(seq); ++i)
            {
                res = streamPut(stream, char(33 + 40));
                if (res)
                    return res;
            }
        }
    } else
    {
        if (options & LINEBREAKS)
        {
            // write stream character by character
            typename Iterator<TQualString const>::Type it = begin(qual);
            typename Iterator<TQualString const>::Type it_end = end(qual);
            for (unsigned long l = 0; it < it_end; ++it)
            {
                res = streamPut(stream, *it);
                if (res)
                    return res;
                if (++l == 70)
                {
                    res = streamPut(stream, '\n');
                    l = 0;
                    if (res)
                        return res;
                }
            }
            if (res)
                return res;
        } else
        {
            res = streamPut(stream, qual);
            if (res)
                return res;
        }
    }
    return 0;
}

// FASTQ
template <typename TStream,
          typename TIdString,
          typename TQualString,
          typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            TQualString const & qual,
            Fastq const & /*tag*/,
            unsigned const options)
{
    int res = streamPut(stream, '@');
    if (res)
        return res;

    res = streamPut(stream, meta);
    if (res)
        return res;

    res = streamPut(stream, '\n');
    if (res)
        return res;

    if (options & LINEBREAKS)
    {
        // write stream character by character
        typename Iterator<TSeqString const>::Type it = begin(seq);
        typename Iterator<TSeqString const>::Type it_end = end(seq);
        for (unsigned long l = 0; it < it_end; ++it)
        {
            res = streamPut(stream, *it);
            if (res)
                return res;
            if (++l == 70)
            {
                res = streamPut(stream, '\n');
                l = 0;
                if (res)
                    return res;
            }
        }
        if (res)
            return res;
    } else
    {
        res = streamPut(stream, seq);
        if (res)
            return res;
    }

    res = streamPut(stream, "\n+");
    if (res)
        return res;

    if (options & WRITEQUALITIESMETA)
    {
        res = streamPut(stream, meta);
        if (res)
            return res;
    }

    res = streamPut(stream, "\n");
    if (res)
        return res;

    res = _writeRecordFastq(stream, seq, qual, options, typename HasQualities<typename Value<TSeqString>::Type>::Type());
    if (res)
        return res;
    res = streamPut(stream, '\n');
    return res;
}

// FASTQ and we have no qualities
template <typename TStream,
          typename TIdString,
          typename TQualString,
          typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            TQualString const & qual,
            Fastq const & /*tag*/)
{
    return writeRecord(stream, meta, seq, qual, Fastq(), DEFAULT_FASTQ);
}

/**
.Function.writeRecord
..signature:writeRecord(TStream & stream, TIdString const & meta, TSeqString const & seq, Fastq const &[, FastAQOutputOptions const options])
..remarks:Writing to FASTQ without specifying qualities currently will result in all qualities being set to maximum ('\126')
*/
// FASTQ and we don't have the qualities
template <typename TStream,
          typename TIdString,
          typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            Fastq const & /*tag*/,
            unsigned const options)
{
    return writeRecord(stream, meta, seq, CharString(), Fastq(), options);
}
// FASTQ and we don't have the qualities
template <typename TStream,
          typename TIdString,
          typename TSeqString>
inline int
writeRecord(TStream & stream,
            TIdString const & meta,
            TSeqString const & seq,
            Fastq const & /*tag*/)
{
    return writeRecord(stream, meta, seq, CharString(), Fastq(), DEFAULT_FASTQ);
}



// ----------------------------------------------------------------------------
// Function write2()
// ----------------------------------------------------------------------------


/**
.Function.write2
..signature:write2(TStream & stream, StringSet<TIdString, TIdSpec> & sequenceIds, StringSet<TSeqString, TSeqSpec> & sequences, Fasta const &[, FastAQOutputOptions const options])
..param.options:if not supplied defaults are chosen, see @Enum.FastAQOutputOptions@
...type:Enum.FastAQOutputOptions
*/
// FASTA
template <typename TStream,
          typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> const & sequenceIds,
         StringSet<TSeqString, TSeqSpec> const & sequences,
         Fasta const & /*tag*/,
         FastAQOutputOptions const options)
{
    if (length(sequenceIds) != length(sequences))
        return -1;

    typedef StringSet<TIdString, TIdSpec> const TIdSet;
    typedef StringSet<TSeqString, TSeqSpec> const TSeqSet;

    typename Iterator<TIdSet>::Type  itMeta     = begin(sequenceIds);
    typename Iterator<TIdSet>::Type  itMeta_end = end(sequenceIds);
    typename Iterator<TSeqSet>::Type itSeq      = begin(sequences);
//    typename Iterator<TSeqSet>::Type itSeq_end  = end(sequences);

    for (; itMeta != itMeta_end; ++itMeta, ++itSeq)
    {
        int res = writeRecord(stream, *itMeta, *itSeq, Fasta(), options);
        if (res)
            return res;
    }
    return 0;
}

template <typename TStream,
          typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> const & sequenceIds,
         StringSet<TSeqString, TSeqSpec> const & sequences,
         Fasta const & /*tag*/)
{
    return write2(stream, sequenceIds, sequences, Fasta(), DEFAULT_FASTA);

}

/**
.Function.write2
..signature:write2(TStream & stream, StringSet<TIdString, TIdSpec> & sequenceIds, StringSet<TSeqString, TSeqSpec> & sequences, [StringSet<TQualString, TQualSpec> & qualities, ]Fastq const &[, FastAQOutputOptions const options])
*/
// FASTQ
template <typename TStream,
          typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TQualString, typename TQualSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> const & sequenceIds,
         StringSet<TSeqString, TSeqSpec> const & sequences,
         StringSet<TQualString, TQualSpec> const & qualities,
         Fastq const & /*tag*/,
         FastAQOutputOptions const options)
{
    if (length(sequenceIds) != length(sequences) ||
        length(qualities) != length(sequences))
        return -1;

    typedef StringSet<TIdString, TIdSpec> const TIdSet;
    typedef StringSet<TSeqString, TSeqSpec> const TSeqSet;
    typedef StringSet<TQualString, TSeqSpec> const TQualSet;

    typename Iterator<TIdSet>::Type   itMeta      = begin(sequenceIds);
    typename Iterator<TIdSet>::Type   itMeta_end  = end(sequenceIds);
    typename Iterator<TSeqSet>::Type  itSeq       = begin(sequences);
//    typename Iterator<TSeqSet>::Type  itSeq_end   = end(sequences);
    typename Iterator<TQualSet>::Type itQual      = begin(qualities);
    // typename Iterator<TQualSet>::Type itQual_end  = end(qualities);

    for (; itMeta != itMeta_end; ++itMeta, ++itSeq, ++itQual)
    {
        int res = writeRecord(stream,
                              *itMeta, *itSeq, *itQual,
                              Fastq(),
                              options);
        if (res)
            return res;
    }
    return 0;
}

template <typename TStream,
          typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TQualString, typename TQualSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> const & sequenceIds,
         StringSet<TSeqString, TSeqSpec> const & sequences,
         StringSet<TQualString, TQualSpec> const & qualities,
         Fastq const & /*tag*/)
{
    return write2(stream,
                  sequenceIds, sequences, qualities,
                  Fastq(),
                  DEFAULT_FASTQ);
}

template <typename TStream,
          typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> const & sequenceIds,
         StringSet<TSeqString, TSeqSpec> const & sequences,
         Fastq const & /*tag*/,
         FastAQOutputOptions const options)
{
    typedef StringSet<TIdString, TIdSpec> const TIdSet;
    typedef StringSet<TSeqString, TSeqSpec> const TSeqSet;

    typename Iterator<TIdSet>::Type   itMeta      = begin(sequenceIds);
    typename Iterator<TIdSet>::Type   itMeta_end  = end(sequenceIds);
    typename Iterator<TSeqSet>::Type  itSeq       = begin(sequences);
//    typename Iterator<TSeqSet>::Type  itSeq_end   = end(sequences);

    for (; itMeta != itMeta_end; ++itMeta, ++itSeq)
    {
        int res = writeRecord(stream,
                              *itMeta, *itSeq,
                              Fastq(),
                              options);
        if (res)
            return res;
    }
    return 0;
}

template <typename TStream,
          typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec>
int write2(TStream & stream,
         StringSet<TIdString, TIdSpec> const & sequenceIds,
         StringSet<TSeqString, TSeqSpec> const & sequences,
         Fastq const & /*tag*/)
{
    return write2(stream,
                  sequenceIds, sequences,
                  Fastq(),
                  DEFAULT_FASTQ);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQ_IO_WRITE_H_
