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
    return "'_matrixName not implemented for this matrix'";
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
          BlastFormatGeneration g>
inline void
_statsBlock(char                      * buffer,
            TMatch              const & m,
            BlastFormat<BlastFormatFile::PAIRWISE,
                        BlastFormatProgram::BLASTN,
                        g>      const & /*tag*/)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n"
                   " Strand=", // no spaces here for whatever reason
            m.bitScore, unsigned(m.score), m.eValue,
                        m.identities, m.aliLength,
            int(ROUND(double(m.identities) * 100 / m.aliLength)),
            m.gaps,       m.aliLength,
            int(ROUND(double(m.gaps)       * 100 / m.aliLength)));
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
          BlastFormatGeneration g>
inline void
_statsBlock(char                      * buffer,
            TMatch              const & m,
            BlastFormat<BlastFormatFile::PAIRWISE,
                        BlastFormatProgram::BLASTP,
                        g>      const & /*tag*/)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Positives = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n\n",
            m.bitScore, unsigned(m.score), m.eValue,
            m.identities, m.aliLength,
            int(ROUND(double(m.identities) * 100 / m.aliLength)),
            m.positives,  m.aliLength,
            int(ROUND(double(m.positives)  * 100 / m.aliLength)),
            m.gaps,       m.aliLength,
            int(ROUND(double(m.gaps)       * 100 / m.aliLength)));
}

template <typename TMatch,
          BlastFormatGeneration g>
inline void
_statsBlock(char                     * buffer,
            TMatch             const & m,
            BlastFormat<BlastFormatFile::PAIRWISE,
                        BlastFormatProgram::BLASTX,
                        g>     const & /*tag*/)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Positives = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n"
                   " Frame = %+d\n\n",
            m.bitScore, unsigned(m.score), m.eValue,
            m.identities, m.aliLength,
            int(ROUND(double(m.identities) * 100 / m.aliLength)),
            m.positives,  m.aliLength,
            int(ROUND(double(m.positives)  * 100 / m.aliLength)),
            m.gaps,       m.aliLength,
            int(ROUND(double(m.gaps)       * 100 / m.aliLength)),
            m.qFrameShift);
}

template <typename TMatch,
          BlastFormatGeneration g>
inline void
_statsBlock(char                    * buffer,
            TMatch            const & m,
            BlastFormat<BlastFormatFile::PAIRWISE,
                        BlastFormatProgram::TBLASTN,
                        g>    const & /*tag*/)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Positives = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n"
                   " Frame = %+d\n\n",
            m.bitScore, unsigned(m.score), m.eValue,
            m.identities, m.aliLength,
            int(ROUND(double(m.identities) * 100 / m.aliLength)),
            m.positives,  m.aliLength,
            int(ROUND(double(m.positives)  * 100 / m.aliLength)),
            m.gaps,       m.aliLength,
            int(ROUND(double(m.gaps)       * 100 / m.aliLength)),
            m.sFrameShift);
}

template <typename TMatch,
          BlastFormatGeneration g>
inline void
_statsBlock(char                    * buffer,
            TMatch            const & m,
            BlastFormat<BlastFormatFile::PAIRWISE,
                        BlastFormatProgram::TBLASTX,
                        g>    const & /*tag*/)
{
    sprintf(buffer," Score =  %.1f bits (%d), Expect =  %.1g\n"
                   " Identities = %d/%d (%d%%),"
                   " Positives = %d/%d (%d%%),"
                   " Gaps = %d/%d (%d%%)\n"
                   " Frame = %+d/%+d\n\n",
            m.bitScore, unsigned(m.score), m.eValue,
            m.identities, m.aliLength,
            int(ROUND(double(m.identities) * 100 / m.aliLength)),
            m.positives,  m.aliLength,
            int(ROUND(double(m.positives)  * 100 / m.aliLength)),
            m.gaps,       m.aliLength,
            int(ROUND(double(m.gaps)       * 100 / m.aliLength)),
            m.qFrameShift, m.sFrameShift);
    //TODO there is an N-column beside e-value here, whats that?
}

template <typename TStream,
          typename TChar1,
          typename TChar2,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeAlignmentBlockIntermediateChar(TStream                 & stream,
                                     TChar1            const & char1,
                                     TChar2            const & char2,
                                     BlastFormat<BlastFormatFile::PAIRWISE,
                                                p,
                                                g>     const & /*tag*/)
{
    if ((char1 == '-') || (char2 == '-'))
        write(stream, ' ');
    else if (char1 == char2)
        write(stream, char1);
    //TODO softcode scoring scheme
    else if (score(Blosum62(), char1, char2) > 0)
        write(stream, '+');
    else
        write(stream, ' ');
}

template <typename TStream,
          typename TChar1,
          typename TChar2,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeAlignmentBlockIntermediateChar(TStream                 & stream,
                                     TChar1                  & char1,
                                     TChar2                  & char2,
                                     BlastFormat<BlastFormatFile::PAIRWISE,
                                                 BlastFormatProgram::BLASTN,
                                                 g>    const & /*tag*/)
{
    if (char1 == '-' || char2 == '-')
        write(stream, ' ');
    else if (char1 == char2)
        write(stream, '|');
    else
        write(stream, ' ');
}

template <typename T>
constexpr typename std::enable_if<std::is_integral<T>::value, T>::type
_numberOfDigits(T const number)
{
    return (number == 0) ? 1 : static_cast<T>(std::floor(std::log10(number))+1);
}


template <typename TStream,
          typename TMatch,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeAlignmentBlock(TStream                 & stream,
                     TMatch            const & m,
                     BlastFormat<BlastFormatFile::PAIRWISE,
                                 p,
                                 g>    const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::PAIRWISE,p,g> TFormat;
    typedef decltype(m.qStart) TPos;

    TPos    const   windowSize  = 60;

    char            buffer[40]  = "";

    TPos            aPos        = 0; // position in alignment
    uint32_t        qPos        = 0; // position in query (without gaps)
    uint32_t        sPos        = 0; // position in subject (without gaps)
    // the latter two can become negative

    TPos            effQStart   = m.qStart;
    TPos            effQEnd     = m.qEnd;
    TPos            effSStart   = m.sStart;
    TPos            effSEnd     = m.sEnd;

    _untranslatePositions(effQStart, effQEnd, m.qFrameShift, m.qLength,
                          QHasRevComp<TFormat>(), QIsTranslated<TFormat>());
    _untranslatePositions(effSStart, effSEnd, m.sFrameShift, m.sLength,
                          SHasRevComp<TFormat>(), SIsTranslated<TFormat>());

    int8_t const     qStep = _step(m.qFrameShift,
                                   QHasRevComp<TFormat>(),
                                   QIsTranslated<TFormat>());
    int8_t const     sStep = _step(m.sFrameShift,
                                   SHasRevComp<TFormat>(),
                                   SIsTranslated<TFormat>());
    int8_t const  qStepOne = _step(m.qFrameShift,
                                   QHasRevComp<TFormat>(),
                                   False());
    int8_t const  sStepOne = _step(m.sFrameShift,
                                   SHasRevComp<TFormat>(),
                                   False());

    auto    const & row0        = row(m.align, 0);
    auto    const & row1        = row(m.align, 1);

    TPos    const   maxPos      = std::max({effQStart, effQEnd, effSStart,
                                            effSEnd});
    // max # digits in pos's
    unsigned char const numberWidth = _numberOfDigits(maxPos);

//     std::cout << "m.aliLength: " << m.aliLength
//               << "\t length(row0): " << length(row0)
//               << "\t length(row1): " << length(row1)
//               << "\n";

    while (aPos < m.aliLength)
    {
        // Query line
        sprintf(buffer, "Query  %-*d  ", numberWidth, qPos + effQStart);
        write(stream, buffer);

        TPos const end = _min(aPos + windowSize, m.aliLength);
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
            _writeAlignmentBlockIntermediateChar(stream,
                                                 value(row0,i),
                                                 value(row1,i),
                                                 TFormat());

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
}

template <typename TStream,
          typename TMatch,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeFullMatch(TStream             & stream,
                TMatch        const & m,
                BlastFormat<BlastFormatFile::PAIRWISE,
                            p,
                            g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::PAIRWISE,p,g> TFormat;
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
    _statsBlock(buffer, m, TFormat());
    write(stream, buffer);

    _writeAlignmentBlock(stream, m, TFormat());
}

template <typename TStream,
          typename TMatch,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeMatchOneLiner(TStream             & stream,
                   TMatch        const & m,
                   BlastFormat<BlastFormatFile::PAIRWISE,
                               p,
                               g> const & /*tag*/)
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

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

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

template <typename TStream,
          typename TDbSpecs,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeTop(TStream                                            & stream,
         TDbSpecs                                     const & dbSpecs,
         BlastFormat<BlastFormatFile::PAIRWISE, p, g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::PAIRWISE, p, g> TFormat;

    // write TOP
    write(stream, _programTagToString(TFormat()));
    write(stream, " I/O Module of SeqAn-");
    write(stream, SEQAN_VERSION_MAJOR);
    write(stream, '.');
    write(stream, SEQAN_VERSION_MINOR);
    write(stream, '.');
    write(stream, SEQAN_VERSION_PATCH);
    write(stream, " (http://www.seqan.de)\n\n");

    // write references
    write(stream, _blastReference());
    write(stream, _seqanReference());

    write(stream, "\n\nDatabase: ");
    write(stream, dbSpecs.dbName);
    write(stream, "\n           ");
    char buffer[40] = "";
//     sprintf(buffer,
//             _uint_label(dbSpecs.dbNumberOfSeqs),
//             dbSpecs.dbNumberOfSeqs); //TODO insert commata
//     write(stream, buffer);
    write(stream, dbSpecs.dbNumberOfSeqs);
    write(stream, " sequences; ");
        clear(buffer);
//     sprintf(buffer,
//             _uint_label(dbSpecs.dbTotalLength),
//             dbSpecs.dbTotalLength); //TODO insert commata
//     write(stream, buffer);
    write(stream, dbSpecs.dbTotalLength);
    write(stream, " total letters\n\n");
}

template <typename TStream,
          typename TRecord,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeRecordHeader(TStream               & stream,
                    TRecord        const & record,
                    BlastFormat<BlastFormatFile::PAIRWISE,
                                p,
                                g> const & /*tag*/)
{
    // write query header
    write(stream, "\nQuery= ");
    write(stream, record.qId);
    write(stream, "\n\nLength=");
    write(stream, record.qLength);

    write(stream,
    "\n\n\n                                                                   "
    "Score     E\n"
    "Sequences producing significant alignments:                       "
    "(Bits)  Value\n\n");
}

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeRecord(TStream                                            & stream,
            TRecord                                      const & record,
            TDbSpecs                                     const & /**/,
            BlastFormat<BlastFormatFile::PAIRWISE, p, g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::PAIRWISE, p, g> TFormat;

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

    _writeRecordHeader(stream, record, TFormat());

    // match one-liners
    for (auto const & m : record.matches)
        _writeMatchOneLiner(stream, m, TFormat());

    write(stream, "\nALIGNMENTS\n");
    // full matches
    for (auto const & m : record.matches)
        _writeFullMatch(stream, m, TFormat());

}

template <typename TStream,
          typename TDbSpecs,
          typename TScore,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeBottom(TStream                                           & stream,
            TDbSpecs                                    const & dbSpecs,
            BlastScoringAdapter<TScore>                 const & adapter,
            BlastFormat<BlastFormatFile::PAIRWISE, p,g> const & /*tag*/)
{
    (void)dbSpecs; //TODO add database specs

    TScore scheme(getScoreScheme(adapter));
    seqanScoringScheme2blastScoringScheme(scheme);
    write(stream, "\nMatrix:");
    write(stream, _matrixName(scheme));
    write(stream, "\nGap Penalties: Existence: ");
    write(stream, scoreGapOpen(scheme));
    write(stream, ", Extension: ");
    write(stream, scoreGapExtend(scheme));
    write(stream, "\n\n");
}

} // namespace seqan

#endif // header guard
