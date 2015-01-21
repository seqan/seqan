// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2014, Knut Reinert, FU Berlin
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
// This file contains routines to generate BLAST tab-seperated output
// ==========================================================================

#ifndef SEQAN_EXTRAS_BLAST_WRITE_BLAST_TABULAR_H_
#define SEQAN_EXTRAS_BLAST_WRITE_BLAST_TABULAR_H_

#include <cinttypes>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// BlastMatchField
// ----------------------------------------------------------------------------

/*!
 * @class BlastMatchField BlastMatchField
 * @brief A "meta" datastructure that contains information about members of @link BlastMatch @endlink
 * @headerfile seqan/blast.h
 *
 * @signature template <BlastFormatGeneration g = BlastFormatGeneration::BLAST_PLUS, typename TVoidSpec = void>
 * struct BlastMatchField;
 * @tparam g            The @link BlastFormatGeneration @endlink you are using.
 * @tparam TVoidSpec    An extra spec to prevent global inclusion (you can safely ignore this)
 *
 * This data structure conveniently gives access to all possible fields used in
 * BLAST-compatabile tabular output formats. @link BlastMatchField::Enum @endlink is needed to
 * specify a custom field composition for @link BlastRecord#writeHeader @endlink,
 * @link BlastRecord#writeRecord @endlink and @link BlastMatch#writeMatch @endlink.
 * The member variables offer the correct labels for the tabular formats' I/O and strings
 * to interact with the user on the command line.
 *
 * <h3>Legacy Format</h3>
 *
 * For g == BlastFormatGeneration::@link BlastFormatGeneration::BLAST_LEGACY @endlink
 * only the default columns are supported and the functions mentioned above do
 * not offer an interface with the <tt>fieldList</tt> parameter. This equivalent
 * to the following values for the member variables:
 * @htmlonly
 * <span style="font-size:90%"><table>
 * <tr><th>#</th><th>Enum</th><th>optionLabels</th><th>columnLabels</th><th>descriptions</th><th>implemented</th></tr>
 * <tr><td>0</td><td>STD</td><td>std</td><td>Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score</td><td>Default 12 columns (Query Seq-id, Subject Seq-id, Percentage of identical matches, Alignment length, Number of mismatches, Number of gap openings, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Expect value, Bit score)</td><td>&#9745;</td></tr>
 * </table></span>
 * @endhtmlonly
 *
 * <h3>Current Format</h3>
 *
 * Value Overview for g == BlastFormatGeneration::@link BlastFormatGeneration::BLAST_PLUS @endlink
 *
 * More fields will likely be implemented in the future.
 *
 * @htmlonly
 * <span style="font-size:90%"><table>
 * <tr><th>#</th><th>Enum</th><th>optionLabels</th><th>columnLabels</th><th>descriptions</th><th>implemented</th></tr>
 * <tr><td>0</td><td>STD</td><td>std</td><td>query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score</td><td>Default 12 columns (Query Seq-id, Subject Seq-id, Percentage of identical matches, Alignment length, Number of mismatches, Number of gap openings, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Expect value, Bit score)</td><td>&#9745;</td></tr>
 * <tr><td>1</td><td>Q_SEQ_ID</td><td>qseqid</td><td>query id</td><td>Query Seq-id</td><td>&#9745;</td></tr>
 * <tr><td>2</td><td>Q_GI</td><td>qgi</td><td>query gi</td><td>Query GI</td><td>&#9744;</td></tr>
 * <tr><td>3</td><td>Q_ACC</td><td>qacc</td><td>query acc.</td><td>Query accesion</td><td>&#9744;</td></tr>
 * <tr><td>4</td><td>Q_ACCVER</td><td>qaccver</td><td>query acc.ver</td><td>Query accesion.version</td><td>&#9744;</td></tr>
 * <tr><td>5</td><td>Q_LEN</td><td>qlen</td><td>query length</td><td>Query sequence length</td><td>&#9745;</td></tr>
 * <tr><td>6</td><td>S_SEQ_ID</td><td>sseqid</td><td>subject id</td><td>Subject Seq-id</td><td>&#9745;</td></tr>
 * <tr><td>7</td><td>S_ALL_SEQ_ID</td><td>sallseqid</td><td>subject ids</td><td>All subject Seq-id(s), separated by a ';'</td><td>&#9744;</td></tr>
 * <tr><td>8</td><td>S_GI</td><td>sgi</td><td>subject gi</td><td>Subject GI</td><td>&#9744;</td></tr>
 * <tr><td>9</td><td>S_ALL_GI</td><td>sallgi</td><td>subject gis</td><td>All subject GIs</td><td>&#9744;</td></tr>
 * <tr><td>10</td><td>S_ACC</td><td>sacc</td><td>subject acc.</td><td>Subject accession</td><td>&#9744;</td></tr>
 * <tr><td>11</td><td>S_ACCVER</td><td>saccver</td><td>subject acc.ver</td><td>Subject accession.version</td><td>&#9744;</td></tr>
 * <tr><td>12</td><td>S_ALLACC</td><td>sallacc</td><td>subject accs.</td><td>All subject accessions</td><td>&#9744;</td></tr>
 * <tr><td>13</td><td>S_LEN</td><td>slen</td><td>subject length</td><td>Subject sequence length</td><td>&#9745;</td></tr>
 * <tr><td>14</td><td>Q_START</td><td>qstart</td><td>q. start</td><td>Start of alignment in query</td><td>&#9745;</td></tr>
 * <tr><td>15</td><td>Q_END</td><td>qend</td><td>q. end</td><td>End of alignment in query</td><td>&#9745;</td></tr>
 * <tr><td>16</td><td>S_START</td><td>sstart</td><td>s. start</td><td>Start of alignment in subject</td><td>&#9745;</td></tr>
 * <tr><td>17</td><td>S_END</td><td>send</td><td>s. end</td><td>End of alignment in subject</td><td>&#9745;</td></tr>
 * <tr><td>18</td><td>Q_SEQ</td><td>qseq</td><td>query seq</td><td>Aligned part of query sequence</td><td>&#9744;</td></tr>
 * <tr><td>19</td><td>S_SEQ</td><td>sseq</td><td>subject seq</td><td>Aligned part of subject sequence</td><td>&#9744;</td></tr>
 * <tr><td>20</td><td>E_VALUE</td><td>evalue</td><td>evalue</td><td>Expect value</td><td>&#9745;</td></tr>
 * <tr><td>21</td><td>BIT_SCORE</td><td>bitscore</td><td>bit score</td><td>Bit score</td><td>&#9745;</td></tr>
 * <tr><td>22</td><td>SCORE</td><td>score</td><td>score</td><td>Raw score</td><td>&#9745;</td></tr>
 * <tr><td>23</td><td>LENGTH</td><td>length</td><td>alignment length</td><td>Alignment length</td><td>&#9745;</td></tr>
 * <tr><td>24</td><td>P_IDENT</td><td>pident</td><td>% identity</td><td>Percentage of identical matches</td><td>&#9745;</td></tr>
 * <tr><td>25</td><td>N_IDENT</td><td>nident</td><td>identical</td><td>Number of identical matches</td><td>&#9745;</td></tr>
 * <tr><td>26</td><td>MISMATCH</td><td>mismatch</td><td>mismatches</td><td>Number of mismatches</td><td>&#9745;</td></tr>
 * <tr><td>27</td><td>POSITIVE</td><td>positive</td><td>positives</td><td>Number of positive-scoring matches</td><td>&#9745;</td></tr>
 * <tr><td>28</td><td>GAP_OPEN</td><td>gapopen</td><td>gap opens</td><td>Number of gap openings</td><td>&#9745;</td></tr>
 * <tr><td>29</td><td>GAPS</td><td>gaps</td><td>gaps</td><td>Total number of gaps</td><td>&#9745;</td></tr>
 * <tr><td>30</td><td>P_POS</td><td>ppos</td><td>% positives</td><td>Percentage of positive-scoring matches</td><td>&#9745;</td></tr>
 * <tr><td>31</td><td>FRAMES</td><td>frames</td><td>query/sbjct frames</td><td>Query and subject frames separated by a '/'</td><td>&#9745;</td></tr>
 * <tr><td>32</td><td>Q_FRAME</td><td>qframe</td><td>query frame</td><td>Query frame</td><td>&#9745;</td></tr>
 * <tr><td>33</td><td>S_FRAME</td><td>sframe</td><td>sbjct frame</td><td>Subject frame</td><td>&#9745;</td></tr>
 * <tr><td>34</td><td>BTOP</td><td>btop</td><td>BTOP</td><td>Blast traceback operations (BTOP)</td><td>&#9744;</td></tr>
 * <tr><td>35</td><td>S_TAX_IDS</td><td>staxids</td><td>subject tax ids</td><td>unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)</td><td>&#9744;</td></tr>
 * <tr><td>36</td><td>S_SCI_NAMES</td><td>sscinames</td><td>subject sci names</td><td>unique Subject Scientific Name(s), separated by a ';'</td><td>&#9744;</td></tr>
 * <tr><td>37</td><td>S_COM_NAMES</td><td>scomnames</td><td>subject com names</td><td>unique Subject Common Name(s), separated by a ';'</td><td>&#9744;</td></tr>
 * <tr><td>38</td><td>S_BLAST_NAMES</td><td>sblastnames</td><td>subject blast names</td><td>unique Subject Blast Name(s), separated by a ';' (in alphabetical order)</td><td>&#9744;</td></tr>
 * <tr><td>39</td><td>S_S_KINGDOMS</td><td>sskingdoms</td><td>subject super kingdoms</td><td>unique Subject Super Kingdom(s), separated by a ';' (in alphabetical order)</td><td>&#9744;</td></tr>
 * <tr><td>40</td><td>S_TITLE</td><td>stitle</td><td>subject title</td><td>Subject Title</td><td>&#9744;</td></tr>
 * <tr><td>41</td><td>S_ALL_TITLES</td><td>salltitles</td><td>subject titles</td><td>All Subject Title(s), separated by a '<>'</td><td>&#9744;</td></tr>
 * <tr><td>42</td><td>S_STRAND</td><td>sstrand</td><td>subject strand</td><td>Subject Strand</td><td>&#9744;</td></tr>
 * <tr><td>43</td><td>Q_COV_S</td><td>qcovs</td><td>% subject coverage</td><td>Query Coverage Per Subject</td><td>&#9744;</td></tr>
 * <tr><td>45</td><td>Q_COV_HSP</td><td>qcovhsp</td><td>% hsp coverage</td><td>Query Coverage Per HSP</td><td>&#9744;</td></tr>
 * </table></span>
 * @endhtmlonly
 *
 * @var static_constexpr_char_const_*_const BlastMatchField::columnLabels[]
 * @brief An array of CStrings representing the <b>column label</b> of each field in tabular @link BlastFormatFile @endlink
 *
 * @var static_constexpr_char_const_*_const BlastMatchField::optionLabels[]
 * @brief An array of CStrings representing the command line parameter name of each field
 *
 * @var static_constexpr_char_const_*_const BlastMatchField::descriptions[]
 * @brief An array of CStrings representing the human-readable descriptions of each field
 *
 * @var static_constexpr_bool_const BlastMatchField::implemented[]
 * @brief An array of bools revealing whether the Blast I/O module supports printing this field
 *
 * @var static_constexpr_std::array<Enum_const,12> BlastMatchField::defaults
 * @brief An std::array of @link BlastMatchField::Enum @endlink indicating the fields that are default
 *
 * @code{.cpp}
 * static constexpr std::array<Enum const, 12> defaults
 * {
 *     {
 *         Enum::Q_SEQ_ID,
 *         Enum::S_SEQ_ID,
 *         Enum::P_IDENT,
 *         Enum::LENGTH,
 *         Enum::MISMATCH,
 *         Enum::GAP_OPEN,
 *         Enum::Q_START,
 *         Enum::Q_END,
 *         Enum::S_START,
 *         Enum::S_END,
 *         Enum::E_VALUE,
 *         Enum::BIT_SCORE
 *     }
 * };
 * @endcode
 */

/*!
 * @enum BlastMatchField::Enum
 * @headerfile seqan/blast.h
 * @signature enum class BlastMatchField<g, TVoidSpec>::Enum : uint8_t { ... };
 * @brief A strongly typed enum mapping all fields supported by NCBI Blast to an interact
 *
 * The available values are visible in detailed description of @link BlastMatchField @endlink.
 */

template <BlastFormatGeneration g = BlastFormatGeneration::BLAST_PLUS,
          typename TVoidSpec = void>
struct BlastMatchField;

// Only canonical 12 columns supported for legacy blast
template <typename TVoidSpec>
struct BlastMatchField<BlastFormatGeneration::BLAST_LEGACY, TVoidSpec>
{
    enum class Enum : uint8_t
    {
        STD
    };

    // this is what Enum::STD stands for
    static constexpr std::array<Enum const, 1> defaults
    {
        {
            Enum::STD
        }
    };

    static constexpr char const * const columnLabels [] =
    {
        "Query id, Subject id, % identity, alignment length,"
         " mismatches, gap openings, q. start, q. end, s. start, s."
         " end, e-value, bit score"
    };

    static constexpr char const * const optionLabels [] =
    {
        "std"
    };

    static constexpr char const * const descriptions [] =
    {
        "Default 12 columns (Query Seq-id, Subject Seq-id, Percentage of "
         "identical matches, Alignment length, Number of mismatches, Number of "
         "gap openings, Start of alignment in query, End of alignment in query,"
         " Start of alignment in subject, End of alignment in subject, Expect "
         "value, Bit score)"
    };

    static constexpr bool const implemented [] =
    {
        true
    };
};

// declarations
template <typename TVoidSpec>
constexpr char const * const
BlastMatchField<BlastFormatGeneration::BLAST_LEGACY, TVoidSpec>::optionLabels[1];

template <typename TVoidSpec>
constexpr char const * const
BlastMatchField<BlastFormatGeneration::BLAST_LEGACY, TVoidSpec>::columnLabels[1];

template <typename TVoidSpec>
constexpr char const * const
BlastMatchField<BlastFormatGeneration::BLAST_LEGACY, TVoidSpec>::descriptions[1];

template <typename TVoidSpec>
constexpr bool const
BlastMatchField<BlastFormatGeneration::BLAST_LEGACY, TVoidSpec>::implemented[1];

template <typename TVoidSpec>
constexpr std::array<typename BlastMatchField<BlastFormatGeneration::BLAST_LEGACY, TVoidSpec>::Enum const, 1>
BlastMatchField<BlastFormatGeneration::BLAST_LEGACY, TVoidSpec>::defaults;


template <typename TVoidSpec>
struct BlastMatchField<BlastFormatGeneration::BLAST_PLUS, TVoidSpec>
{
    enum class Enum : uint8_t
    {
        STD,
        Q_SEQ_ID,
        Q_GI,
        Q_ACC,
        Q_ACCVER,
        Q_LEN,
        S_SEQ_ID,
        S_ALL_SEQ_ID,
        S_GI,
        S_ALL_GI,
        S_ACC,
        S_ACCVER,
        S_ALLACC,
        S_LEN,
        Q_START,
        Q_END,
        S_START,
        S_END,
        Q_SEQ,
        S_SEQ,
        E_VALUE,
        BIT_SCORE,
        SCORE,
        LENGTH,
        P_IDENT,
        N_IDENT,
        MISMATCH,
        POSITIVE,
        GAP_OPEN,
        GAPS,
        P_POS,
        FRAMES,
        Q_FRAME,
        S_FRAME,
        BTOP,
        S_TAX_IDS,
        S_SCI_NAMES,
        S_COM_NAMES,
        S_BLAST_NAMES,
        S_S_KINGDOMS,
        S_TITLE,
        S_ALL_TITLES,
        S_STRAND,
        Q_COV_S,
        Q_COV_HSP
    };

    // this is what Enum::STD stands for
    static constexpr std::array<Enum const, 12> defaults
    {
        {
            Enum::Q_SEQ_ID,
            Enum::S_SEQ_ID,
            Enum::P_IDENT,
            Enum::LENGTH,
            Enum::MISMATCH,
            Enum::GAP_OPEN,
            Enum::Q_START,
            Enum::Q_END,
            Enum::S_START,
            Enum::S_END,
            Enum::E_VALUE,
            Enum::BIT_SCORE
        }
    };

    static constexpr char const * const optionLabels [] =
    {
        "std",
        "qseqid",
        "qgi",
        "qacc",
        "qaccver",
        "qlen",
        "sseqid",
        "sallseqid",
        "sgi",
        "sallgi",
        "sacc",
        "saccver",
        "sallacc",
        "slen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "qseq",
        "sseq",
        "evalue",
        "bitscore",
        "score",
        "length",
        "pident",
        "nident",
        "mismatch",
        "positive",
        "gapopen",
        "gaps",
        "ppos",
        "frames",
        "qframe",
        "sframe",
        "btop",
        "staxids",
        "sscinames",
        "scomnames",
        "sblastnames",
        "sskingdoms",
        "stitle",
        "salltitles",
        "sstrand",
        "qcovs",
        "qcovhsp"
    };

    static constexpr char const * const columnLabels [] =
    {
        "query id, subject id, % identity, alignment "
         "length, mismatches, gap opens, q. start, q. end, s. "
         "start, s. end, evalue, bit score",
        "query id",
        "query gi",
        "query acc.",
        "query acc.ver",
        "query length",
        "subject id",
        "subject ids",
        "subject gi",
        "subject gis",
        "subject acc.",
        "subject acc.ver",
        "subject accs.",
        "subject length",
        "q. start",
        "q. end",
        "s. start",
        "s. end",
        "query seq",
        "subject seq",
        "evalue",
        "bit score",
        "score",
        "alignment length",
        "% identity",
        "identical",
        "mismatches",
        "positives",
        "gap opens",
        "gaps",
        "% positives",
        "query/sbjct frames",
        "query frame",
        "sbjct frame",
        "BTOP",
        "subject tax ids",
        "subject sci names",
        "subject com names",
        "subject blast names",
        "subject super kingdoms",
        "subject title",
        "subject titles",
        "subject strand",
        "% subject coverage",
        "% hsp coverage"
    };

    static constexpr char const * const descriptions [] =
    {
        "Default 12 columns (Query Seq-id, Subject Seq-id, Percentage of "
         "identical matches, Alignment length, Number of mismatches, Number of "
         "gap openings, Start of alignment in query, End of alignment in query,"
         " Start of alignment in subject, End of alignment in subject, Expect "
         "value, Bit score",
        "Query Seq-id",
        "Query GI",
        "Query accesion",
        "Query accesion.version",
        "Query sequence length",
        "Subject Seq-id",
        "All subject Seq-id(s), separated by a ';'",
        "Subject GI",
        "All subject GIs",
        "Subject accession",
        "Subject accession.version",
        "All subject accessions",
        "Subject sequence length",
        "Start of alignment in query",
        "End of alignment in query",
        "Start of alignment in subject",
        "End of alignment in subject",
        "Aligned part of query sequence",
        "Aligned part of subject sequence",
        "Expect value",
        "Bit score",
        "Raw score",
        "Alignment length",
        "Percentage of identical matches",
        "Number of identical matches",
        "Number of mismatches",
        "Number of positive-scoring matches",
        "Number of gap openings",
        "Total number of gaps",
        "Percentage of positive-scoring matches",
        "Query and subject frames separated by a '/'",
        "Query frame",
        "Subject frame",
        "Blast traceback operations (BTOP)",
        "unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)",
        "unique Subject Scientific Name(s), separated by a ';'",
        "unique Subject Common Name(s), separated by a ';'",
        "unique Subject Blast Name(s), separated by a ';' (in alphabetical order)",
        "unique Subject Super Kingdom(s), separated by a ';' (in alphabetical order)",
        "Subject Title",
        "All Subject Title(s), separated by a '<>'",
        "Subject Strand",
        "Query Coverage Per Subject",
        "Query Coverage Per HSP"
    };

    //TODO(c++14): change to std::bitset that is initialized with binary literal
    static constexpr bool const implemented [] =
    {
        true,
        true,
        false,
        false,
        false,
        true,
        true,
        false,
        false,
        false,
        false,
        false,
        false,
        true,
        true,
        true,
        true,
        true,
        false,
        false,
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false
    };
};

template <typename TVoidSpec>
constexpr char const * const
BlastMatchField<BlastFormatGeneration::BLAST_PLUS, TVoidSpec>::optionLabels[45];

template <typename TVoidSpec>
constexpr char const * const
BlastMatchField<BlastFormatGeneration::BLAST_PLUS, TVoidSpec>::columnLabels[45];

template <typename TVoidSpec>
constexpr char const * const
BlastMatchField<BlastFormatGeneration::BLAST_PLUS, TVoidSpec>::descriptions[45];

template <typename TVoidSpec>
constexpr bool const
BlastMatchField<BlastFormatGeneration::BLAST_PLUS, TVoidSpec>::implemented[45];

template <typename TVoidSpec>
constexpr std::array<typename BlastMatchField<BlastFormatGeneration::BLAST_PLUS, TVoidSpec>::Enum const, 12>
BlastMatchField<BlastFormatGeneration::BLAST_PLUS, TVoidSpec>::defaults;

// ============================================================================
// Metafunctions and global const-expressions
// ============================================================================

// ----------------------------------------------------------------------------
// _seperatorString
// ----------------------------------------------------------------------------


template <BlastFormatProgram p, BlastFormatGeneration g>
constexpr
const char *
_seperatorString(BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                              p,
                              g> const & /*tag*/)
{
    return ", ";
}

template <BlastFormatProgram p, BlastFormatGeneration g>
constexpr
const char *
_seperatorString(BlastFormat<BlastFormatFile::TABULAR,
                              p,
                              g> const & /*tag*/)
{
    return "\t";
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// _firstOcc
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Size<TString>::Type
_firstOcc(TString const & str, typename Value<TString>::Type const & val)
{
    typedef typename Size<TString>::Type S;
    for (S s = 0; s < length(str); ++s)
        if (value(str, s) == val)
            return s;
    return length(str);
}

// ----------------------------------------------------------------------------
// Function _writeField() [match object given]
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeField(TStream & s,
            BlastMatch<TQId, TSId, TPos, TAlign> const & match,
            typename BlastMatchField<g>::Enum const fieldId,
            BlastFormat<f, p, g> const &)
{
    typedef BlastFormat<f, p, g> TFormat;
    typedef typename BlastMatchField<g>::Enum ENUM;
    switch (fieldId)
    {
        case ENUM::STD:
             _writeFields(s,
                          match,
                          BlastMatchField<g>::defaults,
                          TFormat(),
                          false /*appendNewline*/);
            break;
        case ENUM::Q_SEQ_ID:
            write(s, prefix(match.qId, _firstOcc(match.qId, ' ')));
            break;
//         case ENUM::Q_GI: write(s,  * ); break;
//         case ENUM::Q_ACC: write(s,  * ); break;
//         case ENUM::Q_ACCVER: write(s,  * ); break;
        case ENUM::Q_LEN:
            write(s, match.qLength);
            break;
        case ENUM::S_SEQ_ID:
            write(s, prefix(match.sId, _firstOcc(match.sId, ' ')));
            break;
//         case ENUM::S_ALL_SEQ_ID: write(s,  * ); break;
//         case ENUM::S_GI: write(s,  * ); break;
//         case ENUM::S_ALL_GI: write(s,  * ); break;
//         case ENUM::S_ACC: write(s,  * ); break;
//         case ENUM::S_ACCVER: write(s,  * ); break;
//         case ENUM::S_ALLACC: write(s,  * ); break;
        case ENUM::S_LEN:
            write(s, match.sLength);
            break;
        case ENUM::Q_START:
        {
            TPos effectiveQStart    = match.qStart;
            TPos effectiveQEnd      = match.qEnd;
            _untranslatePositions(effectiveQStart, effectiveQEnd,
                                  match.qFrameShift, match.qLength,
                                  QHasRevComp<TFormat>(),
                                  QIsTranslated<TFormat>());
            write(s, effectiveQStart);
        } break;
        case ENUM::Q_END:
        {
            TPos effectiveQStart    = match.qStart;
            TPos effectiveQEnd      = match.qEnd;
            _untranslatePositions(effectiveQStart, effectiveQEnd,
                                  match.qFrameShift, match.qLength,
                                  QHasRevComp<TFormat>(),
                                  QIsTranslated<TFormat>());
            write(s, effectiveQEnd);
        } break;
        case ENUM::S_START:
        {
            TPos effectiveSStart    = match.sStart;
            TPos effectiveSEnd      = match.sEnd;

            _untranslatePositions(effectiveSStart, effectiveSEnd,
                                  match.sFrameShift, match.sLength,
                                  SHasRevComp<TFormat>(),
                                  SIsTranslated<TFormat>());
            write(s, effectiveSStart);
        } break;
        case ENUM::S_END:
        {
            TPos effectiveSStart    = match.sStart;
            TPos effectiveSEnd      = match.sEnd;

            _untranslatePositions(effectiveSStart, effectiveSEnd,
                                  match.sFrameShift, match.sLength,
                                  SHasRevComp<TFormat>(),
                                  SIsTranslated<TFormat>());
            write(s, effectiveSEnd);
        } break;
//         case ENUM::Q_SEQ: write(s,  * ); break;
//         case ENUM::S_SEQ: write(s,  * ); break;
        case ENUM::E_VALUE:
        {
            std::string formatString;
            // imported from NCBI code
            if (match.eValue < 1.0e-180)
                formatString = "%3.1lf";
            else if (match.eValue < 1.0e-99)
                formatString = "%2.0le";
            else if (match.eValue < 0.0009)
                formatString = "%3.0le";
            else if (match.eValue < 0.1)
                formatString = "%4.3lf";
            else if (match.eValue < 1.0)
                formatString = "%3.2lf";
            else if (match.eValue < 10.0)
                formatString = "%2.1lf";
            else
                formatString = "%5.0lf";

            write(s, FormattedNumber<double>(formatString.c_str(),
                                             match.eValue));
        } break;
        case ENUM::BIT_SCORE:
        {
            std::string formatString;
            // imported from NCBI code
            if (match.bitScore > 9999)
                formatString = "%4.3le";
            else if (match.bitScore > 99.9)
                formatString = "%4.0ld";
            else
                formatString = "%4.1lf";

            write(s, FormattedNumber<double>(formatString.c_str(),
                                             match.bitScore));
        } break;
        case ENUM::SCORE:
            write(s, match.score);
            break;
        case ENUM::LENGTH:
            write(s, match.aliLength);
            break;
        case ENUM::P_IDENT:
            write(s, FormattedNumber<double>("%.2f",
                                             double(match.identities) * 100
                                             / match.aliLength));
            break;
        case ENUM::N_IDENT:
            write(s, match.identities);
            break;
        case ENUM::MISMATCH:
            write(s, match.mismatches);
            break;
        case ENUM::POSITIVE:
            write(s, match.positives);
            break;
        case ENUM::GAP_OPEN:
            write(s, match.gapOpenings);
            break;
        case ENUM::GAPS:
            write(s, match.gaps);
            break;
        case ENUM::P_POS:
            write(s, FormattedNumber<double>("%.2f",
                                             double(match.positives) * 100
                                             / match.aliLength));
            break;
        case ENUM::FRAMES:
//             if (match.qFrameShift > 0)
//                 write(s, '+');
            write(s, FormattedNumber<signed char>("%+i", match.qFrameShift));
            write(s, '/');
//             if (match.sFrameShift > 0)
//                 write(s, '+');
            write(s, FormattedNumber<signed char>("%+i", match.sFrameShift));
            break;
        case ENUM::Q_FRAME:
            write(s, FormattedNumber<signed char>("%+i", match.qFrameShift));
            break;
        case ENUM::S_FRAME:
            write(s, FormattedNumber<signed char>("%+i", match.sFrameShift));
            break;
//         case ENUM::BTOP: write( * ); break;
//         case ENUM::S_TAX_IDS: write( * ); break;
//         case ENUM::S_SCI_NAMES: write( * ); break;
//         case ENUM::S_COM_NAMES: write( * ); break;
//         case ENUM::S_BLAST_NAMES: write( * ); break;
//         case ENUM::S_S_KINGDOMS: write( * ); break;
//         case ENUM::S_TITLE: write( * ); break;
//         case ENUM::S_ALL_TITLES: write( * ); break;
//         case ENUM::S_STRAND: write( * ); break;
//         case ENUM::Q_COV_S: write( * ); break;
//         case ENUM::Q_COV_HSP:
        default:
            write(s, "n/i"); // not implemented
    };
}

// ----------------------------------------------------------------------------
// Function _writeFields() [match object given]
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TFieldList,
          BlastFormatFile f,
          BlastFormatProgram p>
inline void
_writeFields(TStream & stream,
             BlastMatch<TQId, TSId, TPos, TAlign> const & match,
             TFieldList const & fields,
             BlastFormat<f, p, BlastFormatGeneration::BLAST_PLUS> const &,
             bool const appendNewline = true)
{
    typedef BlastFormat<f, p, BlastFormatGeneration::BLAST_PLUS> TFormat;

    for (auto it = seqan::begin(fields), itB = it, itEnd = seqan::end(fields);
         it != itEnd;
         ++it)
    {
        if (it != itB)
            write(stream, _seperatorString(TFormat()));

        _writeField(stream, match, *it, TFormat());
    }
    if (appendNewline)
        writeValue(stream, '\n');
}

template <typename TStream,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TFieldList,
          BlastFormatFile f,
          BlastFormatProgram p>
inline void
_writeFields(TStream & stream,
             BlastMatch<TQId, TSId, TPos, TAlign> const & match,
             TFieldList const &,
             BlastFormat<f, p, BlastFormatGeneration::BLAST_LEGACY> const &)
{
    // we override the Generation here and discard the list, since LEGACY
    // only supports defaults anyway
    typedef BlastFormat<f, p, BlastFormatGeneration::BLAST_PLUS> TFormat;
    _writeFields(stream,
                 match,
                 BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::defaults,
                 TFormat());
}

// ----------------------------------------------------------------------------
// Function _writeFieldLabels()
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TFieldList,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeFieldLabels(TStream & stream,
                  TFieldList const & fields,
                  BlastFormat<f, p, g> const &)
{
    typedef BlastFormat<f, p, g> TFormat;
    for (auto it = seqan::begin(fields), itB = it, itEnd = seqan::end(fields);
         it != itEnd;
         ++it)
    {
        if (it != itB)
            write(stream, _seperatorString(TFormat()));

        write(stream, BlastMatchField<g>::columnLabels[static_cast<uint8_t>(*it)]);
    }

    writeValue(stream, '\n');
}

// ----------------------------------------------------------------------------
// Function _writeFields() or labels [no match object given]
// ----------------------------------------------------------------------------


template <typename TStream, BlastFormatFile f, BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeFields(TStream & /**/,
             BlastFormat<f, p, g> const & /*tag*/)
{
}

template <typename TStream, typename TField, typename... TFields,
          BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
inline void
_writeFields(TStream & stream,
             BlastFormat<f, p, g> const & /*tag*/,
             TField const & field1, TFields const & ... fields)
{
    typedef BlastFormat<f, p, g> TFormat;

    write(stream, _seperatorString(TFormat()));
    write(stream, field1);
    _writeFields(stream, TFormat(), fields... );
}

// ----------------------------------------------------------------------------
// Function writeHeader()
// ----------------------------------------------------------------------------

template <typename TStream, typename TString1, typename TString2,
          BlastFormatProgram p, BlastFormatGeneration g>
inline void
_writeHeaderWithoutColumnLabels(TStream & stream,
                                TString1 const & qryId,
                                TString2 const & dbName,
                                BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                            p,
                                            g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,p,g> TFormat;

    write(stream, "# ");
    write(stream, _programTagToString(TFormat()));
    write(stream, " I/O Module of SeqAn-");
    write(stream, SEQAN_VERSION_MAJOR);
    write(stream, '.');
    write(stream, SEQAN_VERSION_MINOR);
    write(stream, '.');
    write(stream, SEQAN_VERSION_PATCH);
    write(stream, " (http://www.seqan.de)\n# Query: ");
    write(stream, qryId);
    write(stream, "\n# Database: ");
    write(stream, dbName);
    write(stream, '\n');
}

/*!
 * @fn BlastFormat#writeHeader
 * @headerfile seqan/blast.h
 * @brief write the header for @link BlastFormatFile::TABULAR_WITH_HEADER @endlink without @link BlastRecord @endlink-object
 * @signature writeHeader(stream, qryId, dbname, [matchCount,] blastFormatTag[, labels...])
 *
 * This function writes the header of a record if @link BlastFormatFile @endlink
 * is TABULAR_WITH_HEADER (for TABULAR this is a no-op). In contrast to
 * @link BlastRecord#writeHeader @endlink no @link BlastRecord @endlink-object
 * is required. Custom column labels can be specified as variadic
 * columnlabel arguments (just pass a printable argument for each label).
 *
 * @param[in,out] stream    The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in] qryId     The full ID of the query sequence (e.g. FASTA
 * one-line description)
 * @param[in] dbName    The name of the database (NOT the ID of a subject sequence,
 * but the name of the database, e.g. "NCBI NR" or path to a file)
 * @param[in] matchCount The amount of matches in the record. Mandatory parameter
 * for @link BlastFormatGeneration::BLAST_PLUS @endlink, optional for
 * @link BlastFormatGeneration::BLAST_LEGACY @endlink
 * @param[in] blastFormatTag The @link BlastFormat @endlink specifier.
 * @param[in] labels...   Optional custom column labels
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastRecord#writeRecord
 */

// default fields
template <typename TStream,
          typename TString,
          typename TString2,
          typename TNumeric,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & stream,
            TString const & qryId,
            TString2 const & dbName,
            TNumeric const hitCount,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> TFormat;

    _writeHeaderWithoutColumnLabels(stream, qryId, dbName, TFormat());

    write(stream, "# Fields: ");
    write(stream,
          BlastMatchField<g>::columnLabels
          [static_cast<uint8_t>(BlastMatchField<g>::Enum::STD)/* == 0 */]);
    write(stream, '\n');

    if (g == BlastFormatGeneration::BLAST_PLUS)
    {
        write(stream, "# ");
        write(stream, hitCount);
        write(stream, " hits found\n");
    }
}

// custom fields
// this is not exposed by the documentation above, but only 
template <typename TStream,
          typename TString,
          typename TString2,
          typename TNumeric,
          typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & stream,
            TString const & qryId,
            TString2 const & dbName,
            TNumeric const hitCount,
            TFieldList const & fieldList,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    static_assert(g == BlastFormatGeneration::BLAST_PLUS,
                  "writeMatch() with custom fields only supported with "
                  "BlastFormatGeneration::BLAST_PLUS; use default fields or "
                  "the third signature (arbitrary columns) instead.");
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;

    _writeHeaderWithoutColumnLabels(stream, qryId, dbName, TFormat());

    write(stream, "# Fields: ");
    _writeFieldLabels(stream, fieldList, TFormat());

    write(stream, "# ");
    write(stream, hitCount);
    write(stream, " hits found\n");
}

// custom columns
template <typename TStream,
          typename TString,
          typename TString2,
          typename TNumeric,
          typename TField,
          typename... TFields,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & stream,
            TString const & qryId,
            TString2 const & dbName,
            TNumeric const hitCount,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /**/,
            TField const & field1,
            TFields const & ... fields)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> TFormat;

    _writeHeaderWithoutColumnLabels(stream, qryId, dbName, TFormat());

    write(stream, "# Fields: ");
    write(stream, field1);
    _writeFields(stream, TFormat(), fields...);
    write(stream, '\n');

    if (g == BlastFormatGeneration::BLAST_PLUS)
    {
        write(stream, "# ");
        write(stream, hitCount);
        write(stream, " hits found\n");
    }
}

// BlastFormatGeneration::BLAST_LEGACY works without hitCount aswell
template <typename TStream,
          typename TString1,
          typename TString2,
          BlastFormatProgram p,
          typename std::enable_if<IsSequence<TString1>::VALUE>::type = 0>
inline void
writeHeader(TStream & stream,
            TString1 const & qryId,
            TString2 const & dbName,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_LEGACY> const & /**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_LEGACY> TFormat;

    writeHeader(stream, qryId, dbName, 0, TFormat());
}

// BlastFormatGeneration::BLAST_LEGACY works without hitCount aswell, custom fields
template <typename TStream,
          typename TString,
          typename TString2,
          typename TField,
          typename... TFields,
          BlastFormatProgram p>
inline void
writeHeader(TStream & stream,
            TString const & qryId,
            TString2 const & dbName,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_LEGACY> const & /**/,
            TField const & field1,
            TFields const & ... fields)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_LEGACY> TFormat;

    writeHeader(stream, qryId, dbName, 0, TFormat(), field1, fields...);
}

/*!
 * @fn BlastRecord#writeHeader
 * @headerfile seqan/blast.h
 * @brief write the header of a @link BlastRecord @endlink to file
 * @signature writeHeader(stream, blastRecord, [fieldList,] blastFormatTag)
 * 
 * This function writes the header of a record if @link BlastFormatFile @endlink
 * is TABULAR_WITH_HEADER (for TABULAR this is a no-op). Custom column labels
 * can be specified by passing a sequence of
 * @link BlastMatchField::Enum @endlink. Use this in conjunction with the same 
 * options of @link BlastMatch#writeMatch @endlink .
 *
 * @param[in,out] stream    The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in] blastRecord    The @link BlastRecord @endlink whose header you want to print.
 * @param[in] fieldList Sequence of @link BlastMatchField::Enum @endlink
 * @param[in] blastFormatTag The @link BlastFormat @endlink specifier.
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastRecord#writeRecord
 */


// Record parameter
template <typename TStream,
          typename TRecord,
          typename TSpec,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & stream,
            TRecord const & record,
            BlastDbSpecs<TSpec> const & dbSpecs,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g>  TFormat;
    writeHeader(stream,
                record.qId,
                dbSpecs.dbName,
                record.matches.size(),
                TFormat());
}

// Record parameter, custom fields
template <typename TStream,
          typename TRecord,
          typename TSpec,
          typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & stream,
            TRecord const & record,
            BlastDbSpecs<TSpec> const & dbSpecs,
            TFieldList const & fieldList,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /**/)
{
    static_assert(g == BlastFormatGeneration::BLAST_PLUS,
                  "writeMatch() with custom fields only supported with "
                  "BlastFormatGeneration::BLAST_PLUS; use default fields or "
                  "the second signature (arbitrary columns) instead.");
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;

    writeHeader(stream, record.qId, dbSpecs.dbName,
                record.matches.size(), fieldList, TFormat());
}

// NO-OPS for TABULAR
template <typename TStream, typename TString1, typename TString2,
          BlastFormatProgram p, BlastFormatGeneration g>
inline void
writeHeader(TStream &, TString1 const &, TString2 const &,
            BlastFormat<BlastFormatFile::TABULAR, p, g> const & /*tag*/)
{
}

template <typename TStream,
          typename T,
          typename T2,
          typename... TFields,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & /**/,
            T const & /**/,
            T2 const & /**/,
            BlastFormat<BlastFormatFile::TABULAR,p,g> const & /*tag*/,
            TFields const & ... /**/)
{
}

template <typename TStream,
          typename T,
          typename T2,
          typename T3,
          typename... TFields,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & /**/,
            T const & /**/,
            T2 const & /**/,
            T3 const & /**/,
            BlastFormat<BlastFormatFile::TABULAR,p,g> const & /*tag*/,
            TFields const & ... /**/)
{
}

template <typename TStream,
          typename T,
          typename T2,
          typename T3,
          typename T4,
          typename... TFields,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & /**/,
            T const & /**/,
            T2 const & /**/,
            T3 const & /**/,
            T4 const & /**/,
            BlastFormat<BlastFormatFile::TABULAR,p,g> const & /*tag*/,
            TFields const & ... /**/)
{
}

// ----------------------------------------------------------------------------
// Function writeMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastMatch#writeMatch
 * @headerfile seqan/blast.h
 * @brief write a @link BlastMatch @endlink to file
 * @signature writeMatch(stream, blastMatch, [fields,] blastFormatTag)
 *
 * This function writes a single match if @link BlastFormatFile @endlink
 * is @link BlastFormatFile::TABULAR_WITH_HEADER @endlink or
 * @link BlastFormatFile::TABULAR @endlink. Not specifying <tt>fields</tt> causes default
 * columns to be printed, which is described <a href="https://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/blastall/blastall_node93.html">
 * here</a>. Please note that BLAST is 1-indexed and considers the last position
 * to be the back, not the end, i.e. last one included in a match/sequence/...,
 * not the one behind it (as SeqAn does); this functions corrects for both of
 * these bahaviours, so you don't have to. Additionally, based on your
 * @link BlastFormatProgram @endlink, positions are transformed back to DNA space, if
 * translation has taken place.
 * Please note also that query and subject IDs are truncated at the first space
 * character in NCBI BLAST, this is also done by default here.
 *
 * By passing the <tt>fields</tt> variable you manually specify which columns you want to
 * print; the same conversions mentioned above will me made. See
 * @link BlastMatchField::Enum @endlink for a list of fields available. <tt>fields</tt>
 * is only available with @link BlastFormatGeneration @endlink ==
 * @link BlastFormatGeneration::BLAST_PLUS @endlink.
 *
 * Many guides recommend always printing the default 12 columns and using only
 * additional columns with additional (custom) data.
 * 
 * Please see @link BlastFormat#writeMatch @endlink for an implementation that
 * does not require a @link BlastMatch @endlink object.
 *
 * @param[in,out] stream    The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in] fields    A @link Sequence @endlink of @link BlastMatchField::Enum @endlink
 * @param[in] blastMatch    The @link BlastMatch @endlink you wish to print.
 * @param[in] blastFormatTag The @link BlastFormat @endlink specifier.
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastRecord#writeRecord
 * @see BlastRecord#writeHeader
 */

// Function for match object and manually specified fields
template <typename TStream,
          typename TBlastMatch,
          typename TBlastMatchFields,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeMatch(TStream & stream,
           TBlastMatch const & match,
           TBlastMatchFields const & fields,
           BlastFormat<BlastFormatFile::TABULAR, p, g> const & /*tag*/)
{
    static_assert(g == BlastFormatGeneration::BLAST_PLUS,
                  "writeMatch() with custom fields only supported with "
                  "BlastFormatGeneration::BLAST_PLUS; use default fields or "
                  "the second signature (arbitrary columns) instead.");
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    _writeFields(stream, match, fields, TFormat());
}

// Function for match object and default colummns
template <typename TStream, typename TBlastMatch,
          BlastFormatProgram p, BlastFormatGeneration g>
inline void
writeMatch(TStream & stream,
           TBlastMatch const & match,
           BlastFormat<BlastFormatFile::TABULAR, p, g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    _writeFields(stream, match, BlastMatchField<g>::defaults, TFormat());
}

// writeMatch TABULAR_WITH_HEADER equal to TABULAR
template <typename TStream,
          typename TBlastMatch,
          typename TBlastMatchFields,
          BlastFormatProgram p>
inline void
writeMatch(TStream & stream,
           TBlastMatch const & match,
           TBlastMatchFields const & fields,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       BlastFormatGeneration::BLAST_PLUS> const &/**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    writeMatch(stream, match, fields, TFormat());
}

template <typename TStream, typename TBlastMatch,
          BlastFormatProgram p, BlastFormatGeneration g>
inline void
writeMatch(TStream & stream,
           TBlastMatch const & match,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &/**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    writeMatch(stream, match, TFormat());
}

/*!
 * @fn BlastFormat#writeMatch
 * @headerfile seqan/blast.h
 * @brief write blast tabular output without a @link BlastMatch @endlink object
 * @signature writeMatch(stream, blastFormatTag, columns...)
 *
 * This function writes a single match if @link BlastFormatFile @endlink
 * is TABULAR_WITH_HEADER or TABULAR. In contrast to 
 * @link BlastMatch#writeMatch @endlink (which works on @link BlastMatch @endlink)
 * this signature allows an arbitrary amount of and
 * arbitrary typed columns to be printed. If you do this and you use the
 * TABULAR_WITH_HEADER format, you should also print custom column labels
 * with @link BlastFormat#writeHeader @endlink. Please note that this function
 * does none of the adjustments mentioned in the other signature,
 * so they have to be done by yourself (if you use
 * the mentioned fields and you want compatability).
 *
 * Many guides recommend always printing the default 12 columns and using only
 * additional columns with additional (custom) data.
 *
 * @param[in,out] stream    The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in] blastFormatTag The @link BlastFormat @endlink specifier.
 * @param[in] columns...   Custom columns
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastRecord#writeRecord
 * @see BlastRecord#writeHeader
 * @see BlastFormat#writeHeader
 */

// Function for arbitrary number and typed fields
template <typename TStream, typename TField, typename... TFields,
          BlastFormatProgram p, BlastFormatGeneration g>
inline void
writeMatch(TStream & stream,
           BlastFormat<BlastFormatFile::TABULAR,
                       p,
                       g> const & /*tag*/,
            TField const & field1,
            TFields const & ... fields)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    write(stream, field1);

    _writeFields(stream, TFormat(), fields...);
    write(stream, '\n');
}

// writeMatch TABULAR_WITH_HEADER equal to TABULAR
template <typename TStream,
          typename TField,
          typename... TFields,
          BlastFormatProgram p, BlastFormatGeneration g>
inline void
writeMatch(TStream & stream,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const & /**/,
           TField const & field1,
           TFields const & ... fields)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    writeMatch(stream, TFormat(), field1, fields... );
}

// ----------------------------------------------------------------------------
// Function writeTop()
// ----------------------------------------------------------------------------

// TABULAR formats have no special "top"

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeRecordImplTab(TStream                    & stream,
                    TRecord              const & record,
                    TDbSpecs             const & dbSpecs,
                    BlastFormat<f, p, g> const & /*tag*/)
{
    typedef BlastFormat<f, p, g> TFormat;

    //TODO if debug, do lots of sanity checks on record

    //NOOP for TABULAR
    writeHeader(stream, record, dbSpecs, TFormat());
    for (auto it = record.matches.begin(); it != record.matches.end(); ++it)
    {
        //SOME SANITY CHECKS
        SEQAN_ASSERT_EQ(CharString(record.qId), CharString(it->qId));

        writeMatch(stream, *it, TFormat());
    }
}

// main dox in blast_base.h; this only contains extra sig for fieldList
/*!
 * @fn BlastRecord#writeRecord
 * @headerfile seqan/blast.h
 * @brief write a @link BlastRecord @endlink including it's
 *  @link BlastMatch @endlink es to a file, manually specifying the fields
 * @signature writeRecord(stream, blastRecord, blastDbSpecs, fieldList, blastFormatTag)
 *
 * @param[in,out] stream        The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in] blastRecord   The @link BlastRecord @endlink you wish to print.
 * @param[in] blastDbSpecs  The @link BlastDbSpecs @endlink .
 * @param[in] blastField    A Sequence of @link BlastMatchField::Enum @endlink to print.
 * @param[in] blastFormatTag The @link BlastFormat @endlink specifier of which
 *  @link BlastFormatFile @endlink must be @link BlastFormatFile::TABULAR @endlink
 *  or @link BlastFormatFile::TABULAR_WITH_HEADER @endlink
 *
 * @see BlastFormat
 * @see BlastRecord
 */

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          typename TFieldList,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeRecordImplTab(TStream                    & stream,
                    TRecord              const & record,
                    TDbSpecs             const & dbSpecs,
                    TFieldList           const & fields,
                    BlastFormat<f, p, g> const & /*tag*/)
{
    typedef BlastFormat<f, p, g> TFormat;

    //TODO if debug, do lots of sanity checks on record

    //NOOP for TABULAR
    writeHeader(stream, record, dbSpecs, fields, TFormat());
    for (auto it = record.matches.begin(); it != record.matches.end(); ++it)
    {
        //SOME SANITY CHECKS
        SEQAN_ASSERT_EQ(CharString(record.qId), CharString(it->qId));

        writeMatch(stream, *it, fields, TFormat());
    }
}

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeRecord(TStream             & stream,
            TRecord       const & record,
            TDbSpecs      const & dbSpecs,
            BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    _writeRecordImplTab(stream, record, dbSpecs, TFormat());
}

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeRecord(TStream             & stream,
            TRecord       const & record,
            TDbSpecs      const & dbSpecs,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    _writeRecordImplTab(stream, record, dbSpecs, TFormat());
}


template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeRecord(TStream             & stream,
            TRecord       const & record,
            TDbSpecs      const & dbSpecs,
            TFieldList    const & fields,
            BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    _writeRecordImplTab(stream, record, dbSpecs, fields, TFormat());
}

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeRecord(TStream             & stream,
            TRecord       const & record,
            TDbSpecs      const & dbSpecs,
            TFieldList    const & fields,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    _writeRecordImplTab(stream, record, dbSpecs, fields, TFormat());
}

// ----------------------------------------------------------------------------
// Function writeBottom()
// ----------------------------------------------------------------------------


// template <typename TStream,
//           typename TDbSpecs,
//           typename TBlastScoringAdapater,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// writeBottom(TStream                           & /**/,
//             TDbSpecs                    const & /**/,
//             TBlastScoringAdapater       const & /**/,
//             BlastFormat<BlastFormatFile::TABULAR, p,g> const & /*tag*/)
// {
//     //TODO check if this is really NO-OP forr this format
// // }
// 
// template <typename TStream,
//           typename TDbSpecs,
//           typename TBlastScoringAdapater,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// writeBottom(TStream                           & /**/,
//             TDbSpecs                    const & /**/,
//             TBlastScoringAdapater       const & /**/,
//             BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p,g> const & /*tag*/)
// {
//     //TODO check if this is really NO-OP forr this format
// // }

} // namespace seqan
#endif // header guard
