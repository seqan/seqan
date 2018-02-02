// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Tag BlastReport
// ----------------------------------------------------------------------------

/*!
 * @class BlastReport
 * @signature typedef Tag<BlastReport_> BlastReport;
 * @headerfile <seqan/blast.h>
 * @brief Support for Blast default file format
 *
 * There are three blast format related tags in SeqAn:
 *
 * <li> @link BlastReport @endlink with the FormattedFile output specialization @link BlastReportFileOut @endlink</li>
 * <li> @link BlastTabular @endlink with the FormattedFile output and input specializations
 * @link BlastTabularFileOut @endlink and @link BlastTabularFileIn @endlink</li>
 * <li> @link BlastTabularLL @endlink which provides light-weight, but very basic tabular IO </li>
 *
 * This is the first tag, it represents support for Blast's default file format (<tt>blastall -m 0</tt> / <tt>blast*
 * -outfmt 0</tt>).
 * Only support for writing is available, see the example below.
 *
 * The reference Blast implementation used for developing the SeqAn support is NCBI Blast+ 2.2.26. In contrast to the
 * tabular format there is no support for writing legacy files (without the +).
 *
 * See @link BlastReportFileOut @endlink for more details.
 */

struct BlastReport_;
typedef Tag<BlastReport_> BlastReport;

// ----------------------------------------------------------------------------
// Tag BlastReportFileOut
// ----------------------------------------------------------------------------

/*!
 * @class BlastReportFileOut
 * @signature template <typename TBlastIOContext = BlastIOContext<>>
 * using BlastReportFileOut = FormattedFile<BlastReport, Output, TBlastIOContext>;
 * @extends FormattedFileOut
 * @headerfile <seqan/blast.h>
 * @brief FormattedFileOut abstraction for @link BlastReport @endlink
 *
 * This is a @link FormattedFile @endlink specialization for writing @link BlastReport @endlink formats. For details
 * on how to influence the writing of files, see @link BlastIOContext @endlink.
 * Please note that you have to specify the type of the context as a template parameter to BlastReportFileOut.
 *
 * @section Overview
 *
 * <ul>
 * <li> open @link BlastReportFileOut @endlink,</li>
 * <li> configure the @link BlastIOContext @endlink </li>
 * <li> @link BlastReportFileOut#writeHeader @endlink </li>
 * <li> @link BlastReportFileOut#writeRecord @endlink repeated up to n times</li>
 * <li> @link BlastReportFileOut#writeFooter @endlink </li>
 * </ul>
 *
 * The following members of the context have to be defined before writing:
 * <ul>
 * <li> @link BlastIOContext::blastProgram @endlink (unless fixed at compile-time)</li>
 * <li> @link BlastIOContext::scoringScheme @endlink</li>
 * <li> @link BlastIOContext::dbName @endlink</li>
 * <li> @link BlastIOContext::dbTotalLength @endlink</li>
 * <li> @link BlastIOContext::dbNumberOfSeqs @endlink</li>
 * </ul>
 *
 * For a detailed example have a look at the
 * <a href="http://seqan.readthedocs.io/en/develop/Tutorial/InputOutput/BlastIO.html">Blast IO tutorial</a>.
 *
 * @see BlastRecord
 */

template <typename TBlastIOContext = BlastIOContext<> >
using BlastReportFileOut = FormattedFile<BlastReport, Output, TBlastIOContext>;

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<BlastReport, T> :
    public MagicHeader<Nothing, T> {};

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<BlastReport, T>
{
    static constexpr char const * VALUE[1] =
    {
        ".m0"
    };
};

template <typename T>
constexpr char const * FileExtensions<BlastReport, T>::VALUE[1];

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
    return "Reference for SeqAn: Doering, A., D. Weese, T. Rausch, K. Reinert (2008): "
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
_matrixName(Blosum45 const & /**/)
{
    return "BLOSUM45";
}

constexpr const char *
_matrixName(Blosum62 const & /**/)
{
    return "BLOSUM62";
}

constexpr const char *
_matrixName(Blosum80 const & /**/)
{
    return "BLOSUM80";
}

constexpr const char *
_matrixName(Pam250 const & /**/)
{
    return "PAM250";
}

template <typename TCurTag>
constexpr const char *
_matrixNameTagDispatch(TagList<TCurTag, void> const &,
                       AminoAcidScoreMatrixID const m)
{
    using TUndType = std::underlying_type_t<AminoAcidScoreMatrixID>;
    return (Find<impl::score::MatrixTags, TCurTag>::VALUE == static_cast<TUndType>(m))
            ? _matrixName(Score<int, ScoreMatrix<AminoAcid, TCurTag>>())
            : "ERROR: Matrix name deduction failed, report this as a bug!";
}

template <typename TCurTag,
          typename TRestList>
constexpr const char *
_matrixNameTagDispatch(TagList<TCurTag, TRestList> const &,
                       AminoAcidScoreMatrixID const m)
{
    using TUndType = std::underlying_type_t<AminoAcidScoreMatrixID>;
    return (Find<impl::score::MatrixTags, TCurTag>::VALUE == static_cast<TUndType>(m))
            ? _matrixName(Score<int, ScoreMatrix<AminoAcid, TCurTag>>())
            : _matrixNameTagDispatch(TRestList(), m);
}

template <typename TValue>
constexpr const char *
_matrixName(Score<TValue, ScoreMatrix<AminoAcid, ScoreSpecSelectable> > const & sc)
{
    return _matrixNameTagDispatch(impl::score::MatrixTags(), getScoreMatrixId(sc));
}

// ----------------------------------------------------------------------------
// Function _writeStatsBlock
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TScore,
          typename TMatch,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeStatsBlock(TStream & stream,
                 BlastIOContext<TScore, p, h> const &,
                 TMatch const & m,
                 BlastReport const &)
{
    write(stream, " Score =  ");
    write(stream, FormattedNumber<double>("%.1f", m.bitScore));
    write(stream, " bits (");
    write(stream, unsigned(m.alignStats.alignmentScore));
    write(stream, "), Expect =  ");
    write(stream, FormattedNumber<double>("%.1g", m.eValue));
    write(stream, "\n Identities = ");
    write(stream, m.alignStats.numMatches);
    write(stream, '/');
    write(stream, m.alignStats.alignmentLength);
    write(stream, " (");
    write(stream, int(std::lround(m.alignStats.alignmentIdentity)));
    write(stream, "%), Gaps = ");
    write(stream, m.alignStats.numGaps);
    write(stream, '/');
    write(stream, m.alignStats.alignmentLength);
    write(stream, " (");
    write(stream, int(std::lround(double(m.alignStats.numGaps) * 100 / m.alignStats.alignmentLength)));
    write(stream, "%)\n Strand=");

    if (m.qFrameShift == 1)
        write(stream, "Plus/");
    else
        write(stream, "Minus/");

    if (m.sFrameShift == 1)
        write(stream, "Plus\n\n");
    else
        write(stream, "Minus\n\n");
}

template <typename TStream,
          typename TScoreSpec,
          typename TMatch,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeStatsBlock(TStream & stream,
                 BlastIOContext<Score<int, ScoreMatrix<AminoAcid, TScoreSpec> >, p, h> const & context,
                 TMatch const & m,
                 BlastReport const &)
{
    write(stream, " Score =  ");
    write(stream, FormattedNumber<double>("%.1f", m.bitScore));
    write(stream, " bits (");
    write(stream, unsigned(m.alignStats.alignmentScore));
    write(stream, "), Expect =  ");
    write(stream, FormattedNumber<double>("%.1g", m.eValue));
    write(stream, "\n Identities = ");
    write(stream, m.alignStats.numMatches);
    write(stream, '/');
    write(stream, m.alignStats.alignmentLength);
    write(stream, " (");
    write(stream, int(std::lround(m.alignStats.alignmentIdentity)));
    write(stream, "%), Positives = ");
    write(stream, m.alignStats.numPositiveScores);
    write(stream, '/');
    write(stream, m.alignStats.alignmentLength);
    write(stream, " (");
    write(stream, int(std::lround(m.alignStats.alignmentSimilarity)));
    write(stream, "%), Gaps = ");
    write(stream, m.alignStats.numGaps);
    write(stream, '/');
    write(stream, m.alignStats.alignmentLength);
    write(stream, " (");
    write(stream, int(std::lround(double(m.alignStats.numGaps) * 100 / m.alignStats.alignmentLength)));
    write(stream, "%)");

    if (context.blastProgram != BlastProgram::BLASTP)
        write(stream, "\n Frame = ");

    if (context.blastProgram == BlastProgram::BLASTX)
    {
        write(stream, FormattedNumber<int8_t>("%+d", m.qFrameShift));
    }
    else if (context.blastProgram == BlastProgram::TBLASTN)
    {
        write(stream, FormattedNumber<int8_t>("%+d", m.sFrameShift));
    }
    else if (context.blastProgram == BlastProgram::TBLASTX)
    {
        write(stream, FormattedNumber<int8_t>("%+d", m.qFrameShift));
        write(stream, "/");
        write(stream, FormattedNumber<int8_t>("%+d", m.sFrameShift));
    }
    write(stream, "\n\n");
}

// ----------------------------------------------------------------------------
// Function _writeAlignmentBlockIntermediateChar
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TScore,
          typename TChar1,
          typename TChar2,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeAlignmentBlockIntermediateChar(TStream & stream,
                                     BlastIOContext<TScore, p, h> const &,
                                     TChar1 const & char1,
                                     TChar2 const & char2,
                                     BlastReport const & /*tag*/)
{
    if (char1 == char2)
        write(stream, '|');
    else
        write(stream, ' ');
}

template <typename TStream,
          typename TScoreSpec,
          typename TChar1,
          typename TChar2,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeAlignmentBlockIntermediateChar(TStream & stream,
                                     BlastIOContext<Score<int, ScoreMatrix<AminoAcid, TScoreSpec> >,
                                                    p,
                                                    h> const & context,
                                     TChar1 const & char1,
                                     TChar2 const & char2,
                                     BlastReport const & /*tag*/)
{
    // TODO(h4nn3s): this function could be replaced by a matrix, which would be a lot faster
    // but it would also require a matrix for every scoringScheme

    if (char1 == char2)
        write(stream, char1);
    else if ((char1 == '-') || (char2 == '-'))
        write(stream, ' ');
    else if (score(seqanScheme(context.scoringScheme), char1, char2) > 0)
        write(stream, '+');
    else
        write(stream, ' ');
}

// ----------------------------------------------------------------------------
// Function _writeAlignmentBlock
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TScore,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlignRow0,
          typename TAlignRow1,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeAlignmentBlock(TStream & stream,
                     BlastIOContext<TScore, p, h> & context,
                     BlastMatch<TAlignRow0, TAlignRow1, TPos, TQId, TSId> const & m,
                     BlastReport const & /*tag*/)
{
    TPos const   windowSize  = 60;

    char            buffer[40]  = "";

    TPos            aPos        = 0; // position in alignment
    TPos            qPos        = 0; // position in query (without gaps)
    TPos            sPos        = 0; // position in subject (without gaps)
    // the latter two might become negative in below calculation, but it is safe even for unsigned types

    TPos            effQStart   = m.qStart;
    TPos            effQEnd     = m.qEnd;
    TPos            effSStart   = m.sStart;
    TPos            effSEnd     = m.sEnd;

    _untranslateQPositions(effQStart, effQEnd, m.qFrameShift, m.qLength, context.blastProgram);
    _untranslateSPositions(effSStart, effSEnd, m.sFrameShift, m.sLength, context.blastProgram);

    int8_t const  qStepOne = (m.qFrameShift < 0) ?  -1 : 1;
    int8_t const  sStepOne = (m.sFrameShift < 0) ?  -1 : 1;
    int8_t const     qStep = qIsTranslated(context.blastProgram) ? qStepOne * 3 : qStepOne;
    int8_t const     sStep = sIsTranslated(context.blastProgram) ? sStepOne * 3 : sStepOne;

    TPos const   maxPos      = std::max({effQStart, effQEnd, effSStart, effSEnd});
    // max # digits in pos's
    unsigned char const numberWidth = _numberOfDigits(maxPos);

    while (aPos < m.alignStats.alignmentLength)
    {
        // Query line
        snprintf(buffer, 40, "Query  %-*d  ", numberWidth, qPos + effQStart);
        write(stream, buffer);

        TPos const end = std::min(static_cast<TPos>(aPos + windowSize), m.alignStats.alignmentLength);
        for (TPos i = aPos; i < end; ++i)
        {
            if (!isGap(m.alignRow0, i))
                qPos += qStep;
            write(stream, getValue(m.alignRow0, i));
        }
        snprintf(buffer, 40, "  %-*d", numberWidth, (qPos + effQStart) - qStepOne);
        write(stream, buffer);

        // intermediate line
        write(stream, "\n         ");
        for (unsigned i = 0; i < numberWidth; ++i)
            write(stream, ' ');

        for (TPos i = aPos; i < end; ++i)
            _writeAlignmentBlockIntermediateChar(stream, context, getValue(m.alignRow0,i), getValue(m.alignRow1,i), BlastReport());

        // Subject line
        snprintf(buffer, 40, "\nSbjct  %-*d  ", numberWidth, sPos + effSStart);
        write(stream, buffer);

        for (TPos i = aPos; i < end; ++i)
        {
            if (!isGap(m.alignRow1, i))
                sPos += sStep;
            write(stream, getValue(m.alignRow1, i));
        }
        snprintf(buffer, 40, "  %-*d\n\n", numberWidth, (sPos + effSStart) - sStepOne);
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
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlignRow0,
          typename TAlignRow1,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeFullMatch(TStream & stream,
                BlastIOContext<TScore, p, h> & context,
                BlastMatch<TAlignRow0, TAlignRow1, TPos, TQId, TSId> const & m,
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

    _writeStatsBlock(stream, context, m, BlastReport());

    _writeAlignmentBlock(stream, context, m, BlastReport());
}

// ----------------------------------------------------------------------------
// Function _writeMatchOneLiner
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TScore,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlignRow0,
          typename TAlignRow1,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeMatchOneLiner(TStream & stream,
                    BlastIOContext<TScore, p, h> &,
                    BlastMatch<TAlignRow0, TAlignRow1, TPos, TQId, TSId> const & m,
                    BlastReport const & /*tag*/)
{
    if (length(m.sId) == 66) // it fits
    {
        write(stream, m.sId);
    }
    else if (length(m.sId) < 66) // needs to be padded with ' '
    {
        write(stream, m.sId);
        for (unsigned char i = 0; i < 66 - length(m.sId); ++i)
            write(stream, ' ');
    }
    else // needs to be truncated
    {
        write(stream, prefix(m.sId, 63));
        write(stream, "...");
    }
    write(stream, ' ');

    write(stream, FormattedNumber<long>("%4li", m.bitScore));
    write(stream, "  ");
    write(stream, FormattedNumber<double>("%.1g", m.eValue));
    write(stream, '\n');
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TScore,
          typename ... TSpecs,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeRecordTop(TStream & stream,
                   BlastIOContext<TScore, p, h> &,
                   BlastRecord<TSpecs...> const & record,
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
          typename ... TSpecs,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeRecordFooter(TStream & stream,
                   BlastIOContext<TScore, p, h> & context,
                   BlastRecord<TSpecs...> const & record,
                   BlastReport const & /*tag*/)
{
    write(stream, "\n"
                  "Lambda     K      H\n"
                  "   ");
    write(stream, FormattedNumber<double>("%-4.3f", getLambda(context.scoringScheme)));
    write(stream, "   ");
    write(stream, FormattedNumber<double>("%-5.4f", getKappa(context.scoringScheme)));
    write(stream, "   ");
    write(stream, FormattedNumber<double>("%-5.4f", getH(context.scoringScheme)));
    write(stream, "\n\n"
                  "Gapped\n"
                  "Lambda     K      H\n");
    write(stream, "   ");
    write(stream, FormattedNumber<double>("%-4.3f", getLambda(context.scoringScheme)));
    write(stream, "   ");
    write(stream, FormattedNumber<double>("%-5.4f", getKappa(context.scoringScheme)));
    write(stream, "   ");
    write(stream, FormattedNumber<double>("%-5.4f", getH(context.scoringScheme)));
    write(stream, "\n\n"
                  "Effective search space used: ");
    write(stream, record.qLength * context.dbTotalLength);
    write(stream, "\n\n");
}

/*!
 * @fn BlastReportFileOut#writeRecord
 * @headerfile seqan/blast.h
 * @brief write a @link BlastRecord @endlink including it's @link BlastMatch @endlinkes and possibly comment lines to a file.
 * @signature void writeRecord(blastReportOut, blastRecord);
 *
 * @param[in,out] blastReportOut A @link BlastReportFileOut @endlink formattedFile.
 * @param[in]     blastRecord     The @link BlastRecord @endlink you wish to print.
 *
 * @section Remarks
 *
 * Modifiy the formattedFile's @link BlastIOContext @endlink to set some properties of the output, like the
 * @link BlastIOContext::versionString @endlink or @link BlastIOContext::dbName @endlink.
 *
 * @throw IOError On low-level I/O errors.
 *
 * @see BlastRecord
 * @see BlastIOContext
 */

template <typename TStream,
          typename TScore,
          typename ... TSpecs,
          BlastProgram p,
          BlastTabularSpec h>
inline void
writeRecord(TStream & stream,
            BlastIOContext<TScore, p, h> & context,
            BlastRecord<TSpecs...> const & record,
            BlastReport const & /*tag*/)
{

    #ifdef DEBUG
    {
        CharString str1(record.qId);
        for (auto const & m : record.matches)
        {
            CharString str2(m.qId);
            SEQAN_ASSERT(startsWith(str1, str2));
        }
    }
    #endif// DEBUG

    _writeRecordTop(stream, context, record, BlastReport());

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

template <typename TContext,
          typename ... TSpecs>
inline void
writeRecord(BlastReportFileOut<TContext> & formattedFile,
            BlastRecord<TSpecs...> const & r)
{
    writeRecord(formattedFile.iter, context(formattedFile), r, BlastReport());
}

// ----------------------------------------------------------------------------
// Function writeHeader()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastReportFileOut#writeHeader
 * @headerfile seqan/blast.h
 * @brief Write the header (top-most section) of a BlastReport file.
 * @signature void writeHeader(blastReportOut);
 *
 * @param[in,out] blastReportOut A @link BlastReportFileOut @endlink formattedFile.
 *
 * @throw IOError On low-level I/O errors.
 *
 * @see BlastRecord
 * @see BlastIOContext
 */

template <typename TStream,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline void
writeHeader(TStream & stream,
            BlastIOContext<TScore, p, h> & context,
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
    write(stream, context.dbNumberOfSeqs); //TODO insert commata
    write(stream, " sequences; ");
    write(stream, context.dbTotalLength); //TODO insert commata
    write(stream, " total letters\n\n");
}

template <typename TContext>
inline void
writeHeader(BlastReportFileOut<TContext> & formattedFile)
{
    writeHeader(formattedFile.iter, context(formattedFile), BlastReport());
}

// ----------------------------------------------------------------------------
// Function writeFooter()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastReportFileOut#writeFooter
 * @headerfile seqan/blast.h
 * @brief Write the footer of a BlastReport.
 * @signature void writeFooter(blastReportOut);
 *
 * @param[in,out] blastReportOut A @link BlastReportFileOut @endlink formattedFile.
 *
 * @throw IOError On low-level I/O errors.
 *
 * @see BlastRecord
 * @see BlastIOContext
 */

template <typename TStream>
inline void
writeMatrixName(TStream & stream, SimpleScore const & scheme)
{
    write(stream, " blastn matrix:");
    write(stream, scoreMatch(scheme));
    write(stream, ' ');
    write(stream, scoreMismatch(scheme));
}

template <typename TStream, typename TScheme>
inline void
writeMatrixName(TStream & stream, TScheme const & scheme)
{
    // see top of file
    write(stream, _matrixName(scheme));
}

template <typename TStream,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline void
writeFooter(TStream & stream,
            BlastIOContext<TScore, p, h> & context,
            BlastReport const & /*tag*/)
{
    write(stream, "\n  Database: ");
    write(stream, context.dbName);
    write(stream, "\n  Number of letters in database: ");
    write(stream, context.dbTotalLength);
    write(stream, "\n  Number of sequences in database:  ");
    write(stream, context.dbNumberOfSeqs);
    write(stream, "\n\n\n\n");

    write(stream, "Matrix: ");
    writeMatrixName(stream, seqanScheme(context.scoringScheme));
    write(stream, "\nGap Penalties: Existence: ");
    write(stream, -scoreGapOpenBlast(context.scoringScheme)); // convert scores to penalties
    write(stream, ", Extension: ");
    write(stream, -scoreGapExtend(context.scoringScheme)); // convert scores to penalties
    //TODO possibly add more parameter information
    write(stream, "\n\n");
}

template <typename TContext>
inline void
writeFooter(BlastReportFileOut<TContext> & formattedFile)
{
    writeFooter(formattedFile.iter, context(formattedFile), BlastReport());
}

} // namespace seqan

#endif // header guard
