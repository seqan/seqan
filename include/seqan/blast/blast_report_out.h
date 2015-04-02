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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// This file contains routines to generate BLAST default output
// ==========================================================================

#ifndef SEQAN_EXTRAS_BLAST_WRITE_BLAST_REPORT_H_
#define SEQAN_EXTRAS_BLAST_WRITE_BLAST_REPORT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================
// ----------------------------------------------------------------------------
// Type BlastReportOut
// ----------------------------------------------------------------------------

/*!
 * @tag BlastReport
 * @signature typedef Tag<BlastReport_> BlastReport;
 * @headerfile <seqan/blast.h>
 * @brief Support for Blast default file format
 *
 * This tag represents support for Blast's default file format (<tt>blastall -m 0</tt> / <tt>blast* -outfmt 0</tt>).
 * Only support for reading is available, see @link BlastReportOut @endlink for more details.
 *
 * The reference Blast implementation used for developing the SeqAn support is NCBI Blast+ 2.2.26. In contrast to the
 * tabular format their is no support for writing legacy files (without the +).
 *
 * SeqAn also supports writing the default blast output format, see @link BlastReport @endlink.
 */
struct BlastReport_;
typedef Tag<BlastReport_> BlastReport;



/*!
 * @class BlastReportOut
 * @signature template <typename TBlastIOContext>
 * using BlastReportOut = FormattedFile<BlastReport, Output, TBlastIOContext>;
 * @extends FormattedFileOut
 * @headerfile <seqan/blast.h>
 * @brief FormattedFileOut abstraction for @link BlastReport @endlink
 *
 * @remarks
 *
 * This is a @link FormattedFile @endlink specialization for writing @link BlastReport @endlink formats. For details
 * on how to influence the writing of files , see @link BlastIOContext @endlink.
 * Please note that you have to  the type of the context as a template parameter to BlastReportOut, see the example
 * below.
 *
 *
 * TODO update example
 * @example
 * @code{.cpp}
 * BlastReportOut out("/tmp/example.blast");
 *
 * context(out).dbName = "Legendary Nucleotide Database";
 *
 * BlastRecord<> r;
 * r.qId = "FIRSTREAD abcdefg";
 *
 * for (...)
 * {
 *     BlastMatch<> m;
 *
 *     // "fill" the match object
 *
 *     appendValue(r.matches, m);
 * }
 *
 * writeRecord(out, r);
 * @endcode
 *
 * @see BlastRecord
 */

template <typename TBlastIOContext = BlastIOContext<> >
using BlastReportOut = FormattedFile<BlastReport, Output, TBlastIOContext>;


// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

// TODO adapt magicheader?
template <typename T>
struct MagicHeader<BlastReport, T> :
    public MagicHeader<Nothing, T> {};

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<BlastReport, T>
{
    static constexpr char const * VALUE[5] =
    {
        ".blast",
        ".m0",
        ".bm0"
    };
};

template <typename T>
constexpr char const * FileExtensions<BlastReport, T>::VALUE[5];

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TContext, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<BlastReport, Output, TContext>, TStorageSpec>
{
    typedef TContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TSpec>
struct FileFormat<FormattedFile<BlastReport, Output, TSpec> >
{
    typedef BlastReport Type;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function guessFormat()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool guessFormat(FormattedFile<BlastReport, Output, TSpec> &)
{
    return true;
}

// ----------------------------------------------------------------------------
// Function _numberOfDigits
// ----------------------------------------------------------------------------

template <typename T>
constexpr typename std::enable_if<std::is_integral<T>::value, T>::type
_numberOfDigits(T const number)
{
    return (number == 0) ? 1 : static_cast<T>(std::floor(std::log10(number))+1);
}

// ----------------------------------------------------------------------------
// some string constants
// ----------------------------------------------------------------------------

constexpr const char *
_blastReference()
{
    return "Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro "
    "A. Schaffer,\nJinghui Zhang, Zheng Zhang, Webb Miller, and David J. "
    "Lipman (1997),\n\"Gapped BLAST and PSI-BLAST: a new generation of protein "
    "database search\nprograms\",  Nucleic Acids Res. 25:3389-3402.\n";
}

constexpr const char *
_seqanReference()
{
    return "Reference for SeqAn: Döring, A., D. Weese, T. Rausch, K. Reinert (2008): "
    "SeqAn --\nAn efficient, generic C++ library for sequence analysis. BMC "
    "Bioinformatics,\n9(1), 11. BioMed Central Ltd."
    " doi:10.1186/1471-2105-9-11\n";
}

template <typename T>
constexpr const char *
_matrixName(T const & /**/)
{
    return "'_matrixName not implemented for this matrix'";
}

constexpr const char *
_matrixName(Blosum62 const & /**/)
{
    return "BLOSUM62";
}

//TODO add more matrices

// ----------------------------------------------------------------------------
// Function _statsBlock
// ----------------------------------------------------------------------------

template <typename TMatch>
inline void
_statsBlock(char * buffer,
            TMatch const & m,
            BlastProgram const,
            BlastProgramTagBlastN const &,
            BlastReport const &)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%), Gaps = %d/%d (%d%%)\n"
                   " Strand=", // no spaces here for whatever reason
            m.bitScore, unsigned(m.alignStats.alignmentScore), m.eValue,
            m.alignStats.numMatches, m.alignLength, int(ROUND(m.alignStats.alignmentIdentity)),
            m.alignStats.numGapOpens + m.alignStats.numGapExtensions, m.alignLength,
            int(ROUND(double(m.alignStats.numGapOpens + m.alignStats.numGapExtensions) * 100 / m.alignLength)));
    if (m.qFrameShift == 1)
        strcat(buffer, "Plus/");
    else
        strcat(buffer, "Minus/");

    if (m.sFrameShift == 1)
        strcat(buffer, "Plus\n\n");
    else
        strcat(buffer, "Minus\n\n");
}

template <typename TMatch>
inline void
_statsBlock(char * buffer,
            TMatch const & m,
            BlastProgram const,
            BlastProgramTagBlastP const &,
            BlastReport const &)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Positives = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n\n",
            m.bitScore, unsigned(m.alignStats.alignmentScore), m.eValue,
            m.alignStats.numMatches, m.alignLength,
            int(ROUND(m.alignStats.alignmentIdentity)),
            m.alignStats.numPositiveScores, m.alignLength,
            int(ROUND(m.alignStats.alignmentSimilarity)),
            m.alignStats.numGapOpens + m.alignStats.numGapExtensions, m.alignLength,
            int(ROUND(double(m.alignStats.numGapOpens + m.alignStats.numGapExtensions) * 100 / m.alignLength)));
}

template <typename TMatch>
inline void
_statsBlock(char * buffer,
            TMatch const & m,
            BlastProgram const,
            BlastProgramTagBlastX const &,
            BlastReport const &)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Positives = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n"
                   " Frame = %+d\n\n",
            m.bitScore, unsigned(m.alignStats.alignmentScore), m.eValue,
            m.alignStats.numMatches, m.alignLength,
            int(ROUND(m.alignStats.alignmentIdentity)),
            m.alignStats.numPositiveScores, m.alignLength,
            int(ROUND(m.alignStats.alignmentSimilarity)),
            m.alignStats.numGapOpens + m.alignStats.numGapExtensions, m.alignLength,
            int(ROUND(double(m.alignStats.numGapOpens + m.alignStats.numGapExtensions) * 100 / m.alignLength)),
            m.qFrameShift);
}

template <typename TMatch>
inline void
_statsBlock(char * buffer,
            TMatch const & m,
            BlastProgram const,
            BlastProgramTagTBlastN const &,
            BlastReport const &)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Positives = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n"
                   " Frame = %+d\n\n",
            m.bitScore, unsigned(m.alignStats.alignmentScore), m.eValue,
            m.alignStats.numMatches, m.alignLength,
            int(ROUND(m.alignStats.alignmentIdentity)),
            m.alignStats.numPositiveScores, m.alignLength,
            int(ROUND(m.alignStats.alignmentSimilarity)),
            m.alignStats.numGapOpens + m.alignStats.numGapExtensions, m.alignLength,
            int(ROUND(double(m.alignStats.numGapOpens + m.alignStats.numGapExtensions) * 100 / m.alignLength)),
            m.sFrameShift);
}

template <typename TMatch>
inline void
_statsBlock(char * buffer,
            TMatch const & m,
            BlastProgram const,
            BlastProgramTagTBlastX const &,
            BlastReport const &)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Positives = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n"
                   " Frame = %+d/%+d\n\n",
            m.bitScore, unsigned(m.alignStats.alignmentScore), m.eValue,
            m.alignStats.numMatches, m.alignLength,
            int(ROUND(m.alignStats.alignmentIdentity)),
            m.alignStats.numPositiveScores, m.alignLength,
            int(ROUND(m.alignStats.alignmentSimilarity)),
            m.alignStats.numGapOpens + m.alignStats.numGapExtensions, m.alignLength,
            int(ROUND(double(m.alignStats.numGapOpens + m.alignStats.numGapExtensions) * 100 / m.alignLength)),
            m.qFrameShift, m.sFrameShift);
    //TODO there is an N-column beside e-value here, whats that?
}

template <typename TMatch>
inline void
_statsBlock(char * buffer,
            TMatch const & m,
            BlastProgram const p,
            BlastProgramTagUnknown const &,
            BlastReport const &)
{
    // selection at run-time
    switch(p)
    {
        case BlastProgram::BLASTN:
        {
            typedef BlastProgramTag<BlastProgram::BLASTN> TNewTag;
            _statsBlock(buffer, m, p, TNewTag(), BlastReport());
        } break;
        case BlastProgram::BLASTP:
        {
            typedef BlastProgramTag<BlastProgram::BLASTP> TNewTag;
            _statsBlock(buffer, m, p, TNewTag(), BlastReport());
        } break;
        case BlastProgram::BLASTX:
        {
            typedef BlastProgramTag<BlastProgram::BLASTX> TNewTag;
            _statsBlock(buffer, m, p, TNewTag(), BlastReport());
        } break;
        case BlastProgram::TBLASTN:
        {
            typedef BlastProgramTag<BlastProgram::TBLASTN> TNewTag;
            _statsBlock(buffer, m, p, TNewTag(), BlastReport());
        } break;
        case BlastProgram::TBLASTX:
        {
            typedef BlastProgramTag<BlastProgram::TBLASTX> TNewTag;
            _statsBlock(buffer, m, p, TNewTag(), BlastReport());
        } break;
        case BlastProgram::INVALID:
        case BlastProgram::UNKNOWN:
            SEQAN_FAIL("Invalid or unkown BlastProgram specified. Don't know how to print it.");
            break;
    }
}

// ----------------------------------------------------------------------------
// Function _writeAlignmentBlockIntermediateChar
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TScore,
          typename TConString,
          typename TChar1,
          typename TChar2,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeAlignmentBlockIntermediateChar(TStream & stream,
                                     BlastIOContext<TScore, TConString, p, h> & context,
                                     TChar1 const & char1,
                                     TChar2 const & char2,
                                     BlastReport const & /*tag*/)
{
    if ((char1 == '-') || (char2 == '-'))
        write(stream, ' ');
    else if (char1 == char2)
        write(stream, char1);
    else if (score(context.scoringAdapter.scheme, char1, char2) > 0)
        write(stream, '+');
    else
        write(stream, ' ');
}

template <typename TStream,
          typename TScore,
          typename TConString,
          typename TChar1,
          typename TChar2,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeAlignmentBlockIntermediateChar(TStream & stream,
                                     BlastIOContext<TScore, TConString, p, h> &,
                                     TChar1 & char1,
                                     TChar2 & char2,
                                     BlastReport const & /*tag*/)
{
    if (char1 == '-' || char2 == '-')
        write(stream, ' ');
    else if (char1 == char2)
        write(stream, '|');
    else
        write(stream, ' ');
}

// ----------------------------------------------------------------------------
// Function _writeAlignmentBlock
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TScore,
          typename TConString,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeAlignmentBlock(TStream & stream,
                     BlastIOContext<TScore, TConString, p, h> & context,
                     BlastMatch<TQId, TSId, TPos, TAlign> const & m,
                     BlastReport const & /*tag*/)
{
    TPos const   windowSize  = 60;

    char            buffer[40]  = "";

    TPos            aPos        = 0; // position in alignment
    uint32_t        qPos        = 0; // position in query (without gaps)
    uint32_t        sPos        = 0; // position in subject (without gaps)
    // the latter two can become negative

    TPos            effQStart   = m.qStart;
    TPos            effQEnd     = m.qEnd;
    TPos            effSStart   = m.sStart;
    TPos            effSEnd     = m.sEnd;

    _untranslateQPositions(effQStart, effQEnd, m.qFrameShift, m.qLength, context.blastProgram, BlastProgramTag<p>());
    _untranslateSPositions(effSStart, effSEnd, m.sFrameShift, m.sLength, context.blastProgram, BlastProgramTag<p>());

    int8_t const  qStepOne = (m.qFrameShift < 0) ?  -1 : 1;
    int8_t const  sStepOne = (m.sFrameShift < 0) ?  -1 : 1;
    int8_t const     qStep = qIsTranslated(context.blastProgram, BlastProgramTag<p>()) ? qStepOne * 3 : qStepOne;
    int8_t const     sStep = sIsTranslated(context.blastProgram, BlastProgramTag<p>()) ? sStepOne * 3 : sStepOne;

    auto const & row0        = row(m.align, 0);
    auto const & row1        = row(m.align, 1);

    TPos const   maxPos      = std::max({effQStart, effQEnd, effSStart, effSEnd});
    // max # digits in pos's
    unsigned char const numberWidth = _numberOfDigits(maxPos);

//     std::cout << "m.alignLength: " << m.alignLength
//               << "\t length(row0): " << length(row0)
//               << "\t length(row1): " << length(row1)
//               << "\n";qStepOne

    while (aPos < m.alignLength)
    {
        // Query line
        sprintf(buffer, "Query  %-*d  ", numberWidth, qPos + effQStart);
        write(stream, buffer);

        TPos const end = _min(aPos + windowSize, m.alignLength);
        for (TPos i = aPos; i < end; ++i)
        {
            if (!isGap(row0, i))
                qPos += qStep;
            write(stream, value(row0, i));
        }
        sprintf(buffer, "  %-*d", numberWidth, (qPos + effQStart) - qStepOne);
        write(stream, buffer);

        // intermediate line
        write(stream, "\n         ");
        for (unsigned i = 0; i < numberWidth; ++i)
            write(stream, ' ');

        for (TPos i = aPos; i < end; ++i)
            _writeAlignmentBlockIntermediateChar(stream, context, value(row0,i), value(row1,i), BlastReport());

        // Subject line
        sprintf(buffer, "\nSbjct  %-*d  ", numberWidth, sPos + effSStart);
        write(stream, buffer);

        for (TPos i = aPos; i < end; ++i)
        {
            if (!isGap(row1, i))
                sPos += sStep;
            write(stream, value(row1, i));
        }
        sprintf(buffer, "  %-*d\n\n", numberWidth, (sPos + effSStart) - sStepOne);
        write(stream, buffer);

        aPos = end;
    }
    write(stream, '\n');
}

// ----------------------------------------------------------------------------
// Function _writeFullMatch
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TScore,
          typename TConString,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeFullMatch(TStream & stream,
                BlastIOContext<TScore, TConString, p, h> & context,
                BlastMatch<TQId, TSId, TPos, TAlign> const & m,
                BlastReport const & /*tag*/)
{
    write(stream, "> ");
    for (unsigned beg = 0, end = 0; end < length(m.sId);)
    {
        if (beg == 0)
            end += 64;
        else
            end += 60;

        if (end >= length(m.sId))
            end = length(m.sId);

        write(stream, infix(m.sId, beg, end));
        write(stream, '\n');//            ");

        beg = end;
    }
    write(stream, "Length=");
    write(stream, m.sLength);
    write(stream, "\n\n");

    char buffer[512] = "";
    _statsBlock(buffer, m, context.blastProgram, BlastProgramTag<p>(), BlastReport());
    write(stream, buffer);

    _writeAlignmentBlock(stream, context, m, BlastReport());
}

// ----------------------------------------------------------------------------
// Function _writeMatchOneLiner
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TScore,
          typename TConString,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeMatchOneLiner(TStream & stream,
                    BlastIOContext<TScore, TConString, p, h> &,
                    BlastMatch<TQId, TSId, TPos, TAlign> const & m,
                    BlastReport const & /*tag*/)
{
    if (length(m.sId) == 66) // it fits
    {
        write(stream, m.sId);
    }
    else if (length(m.sId) < 66) // needs to be padded with ' '
    {
        write(stream, m.sId);
        for (unsigned char i = 0; i < 66 -length(m.sId); ++i)
            write(stream, ' ');
    }
    else // needs to be truncated
    {
        write(stream, prefix(m.sId, 63));
        write(stream, "...");
    }
    write(stream, ' ');

    char buffer[20] = "";
    sprintf(buffer, "%4li  %.1g\n", long(m.bitScore), m.eValue);
    write(stream, buffer);
}

// constexpr
// const char * _uint_label(unsigned short)
// {
//     return "%uh";
// }
// 
// constexpr
// const char * _uint_label(unsigned int)
// {
//     return "%u";
// }
// 
// constexpr
// const char * _uint_label(unsigned long)
// {
//     return "%ul";
// }
// 
// constexpr
// const char * _uint_label(unsigned long long)
// {
//     return "%ull"; // requires C99
// }

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TScore,
          typename TConString,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeRecordHeader(TStream & stream,
                   BlastIOContext<TScore, TConString, p, h> &,
                   BlastRecord<TQId, TSId, TPos, TAlign> const & record,
                   BlastReport const & /*tag*/)
{
    // write query header
    write(stream, "\nQuery= ");
    write(stream, record.qId);
    write(stream, "\n\nLength=");
    write(stream, record.qLength);
    write(stream, "\n");
}

template <typename TStream,
          typename TScore,
          typename TConString,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeRecordFooter(TStream & stream,
                   BlastIOContext<TScore, TConString, p, h> & context,
                   BlastRecord<TQId, TSId, TPos, TAlign> const & record,
                   BlastReport const & /*tag*/)
{
    write(stream, "\n"
                  "Lambda     K      H\n"
                  "   ");
    write(stream, FormattedNumber<double>("%-4.3f", getLambda(context.scoringAdapter)));
    write(stream, "   ");
    write(stream, FormattedNumber<double>("%-5.4f", getKappa(context.scoringAdapter)));
    write(stream, "   ");
    write(stream, FormattedNumber<double>("%-5.4f", getH(context.scoringAdapter)));
    write(stream, "\n\n"
                  "Gapped\n"
                  "Lambda     K      H\n");
    write(stream, "   ");
    write(stream, FormattedNumber<double>("%-4.3f", getLambda(context.scoringAdapter)));
    write(stream, "   ");
    write(stream, FormattedNumber<double>("%-5.4f", getKappa(context.scoringAdapter)));
    write(stream, "   ");
    write(stream, FormattedNumber<double>("%-5.4f", getH(context.scoringAdapter)));
    write(stream, "\n\n"
                  "Effective search space used: ");
    write(stream, record.qLength * context.dbTotalLength);
    write(stream, "\n\n");
}

// DOX for this in blast_tabular_out
template <typename TStream,
          typename TScore,
          typename TConString,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastProgram p,
          BlastTabularSpec h>
inline void
writeRecord(TStream & stream,
            BlastIOContext<TScore, TConString, p, h> & context,
            BlastRecord<TQId, TSId, TPos, TAlign> const & record,
            BlastReport const & /*tag*/)
{

    #ifdef DEBUG
    {
        CharString str1(record.qId);
        for (auto const & m : record.matches)
        {
            CharString str2(m.qId);
            SEQAN_ASSERT_EQ(str1, str2);
        }
    }
    #endif// DEBUG

    _writeRecordHeader(stream, context, record, BlastReport());

    if (!empty(record.matches))
    {
            write(stream, "                                                                   Score     E\n"
                          "Sequences producing significant alignments:                       (Bits)  Value\n\n");
        // match one-liners
        for (auto const & m : record.matches)
            _writeMatchOneLiner(stream, context, m, BlastReport());

        write(stream, "\nALIGNMENTS\n");
        // full matches
        for (auto const & m : record.matches)
            _writeFullMatch(stream, context, m, BlastReport());
    }
    else
    {
        write(stream, "\n\n***** No hits found *****\n\n\n");
    }

    _writeRecordFooter(stream, context, record, BlastReport());
}

/*!
 * @fn BlastReportOut#writeRecord
 * @headerfile seqan/blast.h
 * @brief write a @link BlastRecord @endlink including it's @link BlastMatch @endlinkes and possible headers to a file.
 * @signature void writeRecord(blastReportOut, blastRecord);
 *
 * @param[in,out] blastReportOut A @link BlastReportOut @endlink formattedFile.
 * @param[in]     blastRecord     The @link BlastRecord @endlink you wish to print.
 *
 * Modifiy the formattedFile's @link BlastIOContext @endlink to set some properties of the output.
 * Note also that this will effect downstream functions like @link BlastRecord#writeRecordHeader @endlink and
 * @link BlastMatch#writeMatch @endlink!
 *
 * @throw IOError On low-level I/O errors.
 */

template <typename TContext,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign>
inline void
writeRecord(BlastReportOut<TContext> & formattedFile,
            BlastRecord<TQId, TSId, TPos, TAlign> const & r)
{
    writeRecord(formattedFile.iter, context(formattedFile), r, BlastReport());
}

// ----------------------------------------------------------------------------
// Function writeHeader()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastReport#writeHeader
 * @headerfile seqan/blast.h
 * @brief write the header (top-most section) of a BlastReport file
 * @signature void writeHeader(stream, context, blastReport);
 *
 * @param[in,out] stream         The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in,out] context        A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     blastReport   The @link BlastReport @endlink tag.
 */

template <typename TStream,
          typename TScore,
          typename TConString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
writeHeader(TStream & stream,
            BlastIOContext<TScore, TConString, p, h> & context,
            BlastReport const & /*tag*/)
{

    if (empty(context.versionString))
        context._setDefaultVersionString();
    write(stream, context.versionString);
    write(stream, "\n\n\n");

    // write references
    write(stream, _blastReference());
    write(stream, "\n\n\n");
    write(stream, _seqanReference());

    write(stream, "\n\n\nDatabase: ");
    write(stream, context.dbName);
    write(stream, "\n           ");
    char buffer[40] = "";
//     sprintf(buffer,
//             _uint_label(context.dbNumberOfSeqs),
//             context.dbNumberOfSeqs); //TODO insert commata
//     write(stream, buffer);
    write(stream, context.dbNumberOfSeqs);
    write(stream, " sequences; ");
        clear(buffer);
//     sprintf(buffer,
//             _uint_label(context.dbTotalLength),
//             context.dbTotalLength); //TODO insert commata
//     write(stream, buffer);
    write(stream, context.dbTotalLength);
    write(stream, " total letters\n\n");
}

/*!
 * @fn BlastReportOut#writeHeader
 * @headerfile seqan/blast.h
 * @brief write the header (top-most section) of a BlastReport file
 * @signature void writeHeader(blastReportOut);
 *
 * @param[in,out] blastReportOut A @link BlastReportOut @endlink formattedFile.
 */

template <typename TContext>
inline void
writeHeader(BlastReportOut<TContext> & formattedFile)
{
    writeHeader(formattedFile.iter, context(formattedFile), BlastReport());
}

// ----------------------------------------------------------------------------
// Function writeFooter()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastReport#writeFooter
 * @headerfile seqan/blast.h
 * @brief write the footer of a BlastReport file
 * @signature void writeFooter(stream, context, blastReport);
 *
 * @param[in,out] stream         The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in,out] context        A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     blastReport   The @link BlastReport @endlink tag.
 */

template <typename TStream,
          typename TScore,
          typename TConString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
writeFooter(TStream & stream,
            BlastIOContext<TScore, TConString, p, h> & context,
            BlastReport const & /*tag*/)
{
    TScore scheme(getScoringScheme(context));

    write(stream, "\n  Database: ");
    write(stream, context.dbName);
    write(stream, "\n  Number of letters in database: ");
    write(stream, context.dbTotalLength);
    write(stream, "\n  Number of sequences in database:  ");
    write(stream, context.dbNumberOfSeqs);
    write(stream, "\n\n\n\n");

    write(stream, "Matrix:");
    write(stream, _matrixName(scheme));
    write(stream, "\nGap Penalties: Existence: ");
    write(stream, -scoreGapOpen(scheme)); // convert scores to penalties
    write(stream, ", Extension: ");
    write(stream, -scoreGapExtend(scheme)); // convert scores to penalties
    //TODO possibly add more parameter information
    write(stream, "\n\n");
}

/*!
 * @fn BlastReportOut#writeFooter
 * @headerfile seqan/blast.h
 * @brief write the footer of a BlastReport
 * @signature void writeFooter(blastReportOut);
 *
 * @param[in,out] blastReportOut A @link BlastReportOut @endlink formattedFile.
 */

template <typename TContext>
inline void
writeFooter(BlastReportOut<TContext> & formattedFile)
{
    writeFooter(formattedFile.iter, context(formattedFile), BlastReport());
}

//TODO check behaviour regarding empty records and truncating of IDs

} // namespace seqan

#endif // header guard
