// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013, Hannes Hauswedell, FU Berlin
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
// This file contains routines to generate BLAST default output
// ==========================================================================

#ifndef SEQAN_EXTRAS_BLAST_WRITE_BLAST_REPORT_H_
#define SEQAN_EXTRAS_BLAST_WRITE_BLAST_REPORT_H_

#include <cstdio>

#include <seqan/basic.h>
#include <seqan/blast/blast_base.h>
#include <seqan/score.h>


namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

constexpr const char *
_blastReference()
{
    return "Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro "
    "A. Schaffer,\nJinghui Zhang, Zheng Zhang, Webb Miller, and David J. "
    "Lipman (1997),\n\"Gapped BLAST and PSI-BLAST: a new generation of protein "
    "database search\nprograms\",  Nucleic Acids Res. 25:3389-3402.\n\n";
}

constexpr const char *
_seqanReference()
{
    return "Reference for SeqAn: Döring, A., D. Weese, T. Rausch, K. Reinert (2008): "
    "SeqAn --\nAn efficient, generic C++ library for sequence analysis. BMC "
    "Bioinformatics,\n9(1), 11. BioMed Central Ltd."
    " doi:10.1186/1471-2105-9-11\n\n";
}

template <typename T>
constexpr const char *
_matrixName(T const & /**/)
{
    return "UNSCPECIFIED";
}


constexpr const char *
_matrixName(Blosum62 const & /**/)
{
    return "BLOSUM62";
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TMatch,
          BlastFormatOptions::Generation g>
inline void
_statsBlock(char                      * buffer,
            TMatch              const & m,
            BlastFormat<BlastFormatOptions::Pairwise,
                        BlastFormatOptions::BlastN,
                        g>      const & /*tag*/)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n"
                   " Strand=", // no spaces here for whatever reason
            m.bitScore, unsigned(m.score), m.eVal,
            m.identities, m.aliLength, m.identities * 100 / m.aliLength,
            m.gaps,       m.aliLength, m.gaps       * 100 / m.aliLength);
    if (m.qFrameShift == 1)
        strcat(buffer, "Plus/");
    else
        strcat(buffer, "Minus/");

    if (m.sFrameShift == 1)
        strcat(buffer, "Plus\n\n");
    else
        strcat(buffer, "Minus\n\n");
}


template <typename TMatch,
          BlastFormatOptions::Generation g>
inline void
_statsBlock(char                      * buffer,
            TMatch              const & m,
            BlastFormat<BlastFormatOptions::Pairwise,
                        BlastFormatOptions::BlastP,
                        g>      const & /*tag*/)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Positives = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n\n",
            m.bitScore, unsigned(m.score), m.eVal,
            m.identities, m.aliLength, m.identities * 100 / m.aliLength,
            m.positives,  m.aliLength, m.positives  * 100 / m.aliLength,
            m.gaps,       m.aliLength, m.gaps       * 100 / m.aliLength);
}

template <typename TMatch,
          BlastFormatOptions::Generation g>
inline void
_statsBlock(char                     * buffer,
            TMatch             const & m,
            BlastFormat<BlastFormatOptions::Pairwise,
                        BlastFormatOptions::BlastX,
                        g>     const & /*tag*/)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Positives = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n"
                   " Frame = %+d\n\n",
            m.bitScore, unsigned(m.score), m.eVal,
            m.identities, m.aliLength, m.identities * 100 / m.aliLength,
            m.positives,  m.aliLength, m.positives  * 100 / m.aliLength,
            m.gaps,       m.aliLength, m.gaps       * 100 / m.aliLength,
            m.qFrameShift);
}

template <typename TMatch,
          BlastFormatOptions::Generation g>
inline void
_statsBlock(char                    * buffer,
            TMatch            const & m,
            BlastFormat<BlastFormatOptions::Pairwise,
                        BlastFormatOptions::TBlastN,
                        g>    const & /*tag*/)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Positives = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n"
                   " Frame = %+d\n\n",
            m.bitScore, unsigned(m.score), m.eVal,
            m.identities, m.aliLength, m.identities * 100 / m.aliLength,
            m.positives,  m.aliLength, m.positives  * 100 / m.aliLength,
            m.gaps,       m.aliLength, m.gaps       * 100 / m.aliLength,
            m.sFrameShift);
}

template <typename TMatch,
          BlastFormatOptions::Generation g>
inline void
_statsBlock(char                    * buffer,
            TMatch            const & m,
            BlastFormat<BlastFormatOptions::Pairwise,
                        BlastFormatOptions::TBlastX,
                        g>    const & /*tag*/)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Positives = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n"
                   " Frame = %+d/%+d\n\n",
            m.bitScore, unsigned(m.score), m.eVal,
            m.identities, m.aliLength, m.identities * 100 / m.aliLength,
            m.positives,  m.aliLength, m.positives  * 100 / m.aliLength,
            m.gaps,       m.aliLength, m.gaps       * 100 / m.aliLength,
            m.qFrameShift, m.sFrameShift);
            //TODO verify that the order is actually qFS/sFS and
            // not the other way around
    //TODO there is an N-column beside e-value here, whats that?
}


template <typename TStream,
          typename TChar1,
          typename TChar2,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
_writeAlignmentBlockIntermediateChar(TStream                 & stream,
                                     TChar1            const & char1,
                                     TChar2            const & char2,
                                     BlastFormat<BlastFormatOptions::Pairwise,
                                                p,
                                                g>     const & /*tag*/)
{
    int ret = 0;
    if ((char1 == '-') || (char2 == '-'))
        ret = streamPut(stream, ' ');
    else if (char1 == char2)
        ret = streamPut(stream, char1);
    //TODO softcode scoring scheme
    else if (score(Blosum62(), char1, char2) > 0)
        ret = streamPut(stream, '+');
    else
        ret = streamPut(stream, ' ');
    return ret;
}

template <typename TStream,
          typename TChar1,
          typename TChar2,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
_writeAlignmentBlockIntermediateChar(TStream                 & stream,
                                     TChar1                  & char1,
                                     TChar2                  & char2,
                                     BlastFormat<BlastFormatOptions::Pairwise,
                                                 BlastFormatOptions::BlastN,
                                                 g>    const & /*tag*/)
{
    int ret = 0;
    if (char1 == '-' || char2 == '-')
        ret = streamPut(stream, ' ');
    else if (char1 == char2)
        ret = streamPut(stream, '|');
    else
        ret = streamPut(stream, ' ');
    return ret;
}


template <typename TStream,
          typename TMatch,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
_writeAlignmentBlock(TStream                 & stream,
                     TMatch            const & m,
                     BlastFormat<BlastFormatOptions::Pairwise,
                                 p,
                                 g>    const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Pairwise,p,g> TFormat;
    typedef decltype(m.qStart) TPos;

    int             ret         = 0;
    TPos    const   windowSize  = 60;

    char            buffer[40]  = "";

    TPos            aPos        = 0; // position in alignment
    TPos            qPos        = 0; // position in query (without gaps)
    TPos            sPos        = 0; // position in subject (without gaps)

    auto    const & row0        = row(m.align, 0);
    auto    const & row1        = row(m.align, 1);

    unsigned char   numberWidth = 5; //TODO get biggest number of digits in pos's

//     std::cout << "m.aliLength: " << m.aliLength
//               << "\t length(row0): " << length(row0)
//               << "\t length(row1): " << length(row1)
//               << "\n";

    while (aPos < m.aliLength)
    {
        // Query line
        sprintf(buffer, "Query  %-*d  ", numberWidth, qPos + m.qStart);
        ret = streamPut(stream, buffer);
        if (ret)
            return ret;

        TPos const end = _min(aPos + windowSize,TPos(m.aliLength));
        for (TPos i = aPos; i < end; ++i)
        {
            if (!isGap(row0, i))
                ++qPos;
            ret = streamPut(stream, value(row0, i));
            if (ret)
                return ret;
        }
        sprintf(buffer, "  %-*d", numberWidth, qPos + m.qStart - 1);
        ret = streamPut(stream, buffer);
        if (ret)
            return ret;

        // intermediate line
        ret = streamPut(stream, "\n         ");
        if (ret)
            return ret;
        for (unsigned i = 0; i < numberWidth; ++i)
        {
            ret = streamPut(stream, ' ');
            if (ret)
                return ret;
        }
        for (TPos i = aPos; i < end; ++i)
        {
            ret = _writeAlignmentBlockIntermediateChar(stream,
                                                       value(row0,i),
                                                       value(row1,i),
                                                       TFormat());
            if (ret)
                return ret;
        }

        // Subject line
        sprintf(buffer, "\nSbjct  %-*d  ", numberWidth, sPos + m.sStart);
        ret = streamPut(stream, buffer);
        if (ret)
            return ret;

        for (TPos i = aPos; i < end; ++i)
        {
            if (!isGap(row1, i))
                ++sPos;
            ret = streamPut(stream, value(row1, i));
            if (ret)
                return ret;
        }
        sprintf(buffer, "  %-*d\n\n", numberWidth, sPos + m.sStart - 1);
        ret = streamPut(stream, buffer);
        if (ret)
            return ret;

        aPos = end;
    }
    return 0;
}


template <typename TStream,
          typename TMatch,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
_writeFullMatch(TStream             & stream,
                TMatch        const & m,
                BlastFormat<BlastFormatOptions::Pairwise,
                            p,
                            g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Pairwise,p,g> TFormat;
    int ret = streamPut(stream, "> ");
    if (ret)
        return ret;
    for (unsigned beg = 0, end = 0; end < length(m.sId);)
    {
        if (beg == 0)
            end += 64;
        else
            end += 60;

        if (end >= length(m.sId))
            end = length(m.sId);

        ret = streamPut(stream, infix(m.sId, beg, end));
        if (ret)
            return ret;
        ret = streamPut(stream, '\n');//            ");
        if (ret)
            return ret;

        beg = end;
    }
    ret = streamPut(stream, "Length=");
    if (ret)
        return ret;
    ret = streamPut(stream, m.sLength);
    if (ret)
        return ret;
    ret = streamPut(stream, "\n\n");
    if (ret)
        return ret;

    char buffer[512] = "";
    _statsBlock(buffer, m, TFormat());
    ret = streamPut(stream, buffer);
    if (ret)
        return ret;

    ret = _writeAlignmentBlock(stream, m, TFormat());

    return ret;
}

template <typename TStream,
          typename TMatch,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
_writeMatchOneLiner(TStream             & stream,
                   TMatch        const & m,
                   BlastFormat<BlastFormatOptions::Pairwise,
                               p,
                               g> const & /*tag*/)
{
    int ret = 0;
    if (length(m.sId) == 66) // it fits
    {
        ret = streamPut(stream, m.sId);
        if (ret)
            return ret;
    }
    else if (length(m.sId) < 66) // needs to be padded with ' '
    {
        ret = streamPut(stream, m.sId);
        if (ret)
            return ret;
        for (unsigned char i = 0; i < 66 -length(m.sId); ++i)
        {
            ret = streamPut(stream, ' ');
            if (ret)
                return ret;
        }
    }
    else // needs to be truncated
    {
        ret = streamPut(stream, prefix(m.sId, 63));
        if (ret)
            return ret;
        ret = streamPut(stream, "...");
        if (ret)
            return ret;
    }
    ret = streamPut(stream, ' ');
    if (ret)
        return ret;

    char buffer[20] = "";
    sprintf(buffer, "%4li  %.1g\n", long(m.bitScore), m.eVal);
    ret = streamPut(stream, buffer);

    return ret;
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TString,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writeTop(TStream             & stream,
         TString       const & dbName,
         unsigned long const   dbNumSeqs,
         unsigned long const   dbTotalLength,
            BlastFormat<BlastFormatOptions::Pairwise,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Pairwise, p, g> TFormat;

    // write TOP
    int ret = streamPut(stream, _programTagToString(TFormat()));
    if (ret)
        return ret;
    ret = streamPut(stream, " I/O Module of SeqAn-");
    if (ret)
        return ret;
    ret = streamPut(stream, SEQAN_VERSION_MAJOR);
    if (ret)
        return ret;
    ret = streamPut(stream, '.');
    if (ret)
        return ret;
    ret = streamPut(stream, SEQAN_VERSION_MINOR);
    if (ret)
        return ret;
    ret = streamPut(stream, '.');
    if (ret)
        return ret;
    ret = streamPut(stream, SEQAN_VERSION_PATCH);
    if (ret)
        return ret;
    ret = streamPut(stream, " (http://www.seqan.de)\n\n");
    if (ret)
        return ret;

    // write references
    ret = streamPut(stream, _blastReference());
    if (ret)
        return ret;
    ret = streamPut(stream, _seqanReference());
    if (ret)
        return ret;

    ret = streamPut(stream, "\n\nDatabase: ");
    if (ret)
        return ret;
    ret = streamPut(stream, dbName);
    if (ret)
        return ret;
    ret = streamPut(stream, "\n           ");
    if (ret)
        return ret;
    char buffer[40] = "";
    sprintf(buffer, "%lu", dbNumSeqs); //TODO insert commata
    ret = streamPut(stream, buffer);
    if (ret)
        return ret;
    ret = streamPut(stream, " sequences; ");
    if (ret)
        return ret;
        clear(buffer);
    sprintf(buffer, "%lu", dbTotalLength); //TODO insert commata
    ret = streamPut(stream, buffer);
    if (ret)
        return ret;
    ret = streamPut(stream, " total letters\n\n");
    return ret;
}

template <typename TStream,
          typename TRecord,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
_writeRecordHeader(TStream               & stream,
                    TRecord        const & record,
                    BlastFormat<BlastFormatOptions::Pairwise,
                                p,
                                g> const & /*tag*/)
{
    // write query header
    int ret = streamPut(stream, "\nQuery= ");
    if (ret)
        return ret;
    ret = streamPut(stream, record.qId);
    if (ret)
        return ret;
    ret = streamPut(stream, "\n\nLength=");
    if (ret)
        return ret;
    ret = streamPut(stream, record.qLength);
    if (ret)
        return ret;

    ret = streamPut(stream,
    "\n\n\n                                                                   "
    "Score     E\n"
    "Sequences producing significant alignments:                       "
    "(Bits)  Value\n\n");

    return ret;
}

template <typename TStream,
          typename TRecord,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writeRecord(TStream             & stream,
            TRecord       const & record,
            BlastFormat<BlastFormatOptions::Pairwise,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Pairwise, p, g> TFormat;

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

    int ret = _writeRecordHeader(stream, record, TFormat());
    if (ret)
        return ret;

    // match one-liners
    for (auto const & m : record.matches)
    {
        ret = _writeMatchOneLiner(stream, m, TFormat());
        if (ret)
            return ret;
    }
    ret = streamPut(stream, "\nALIGNMENTS\n");
    if (ret)
        return ret;
    // full matches
    for (auto const & m : record.matches)
    {
        ret = _writeFullMatch(stream, m, TFormat());
        if (ret)
            return ret;
    }

    return 0;
}

template <typename TStream,
          typename TValue, typename TSpec,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writeBottom(TStream             & stream,
            Score<TValue, TSpec>   const & score,
            BlastFormat<BlastFormatOptions::Pairwise,
                        p,
                        g> const & /*tag*/)
{
    //TODO add database specs

    int ret = streamPut(stream, "\nMatrix:");
    if (ret)
        return ret;
    ret = streamPut(stream, _matrixName(score));
    if (ret)
        return ret;
    ret = streamPut(stream, "\nGap Penalties: Existence: ");
    if (ret)
        return ret;
    ret = streamPut(stream, scoreGapOpen(score));
    if (ret)
        return ret;
    ret = streamPut(stream, ", Extension: ");
    if (ret)
        return ret;
    ret = streamPut(stream, scoreGapExtend(score));
    if (ret)
        return ret;
    ret = streamPut(stream, "\n\n");
    return ret;
}

} // namespace seqan

#endif // header guard