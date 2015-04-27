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
// This file contains routines to generate BLAST tab-seperated output
// ==========================================================================

#ifndef SEQAN_BLAST_BLAST_TABULAR_H_
#define SEQAN_BLAST_BLAST_TABULAR_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class BlastTabular
 * @signature typedef Tag<BlastTabular_> BlastTabular;
 * @headerfile <seqan/blast.h>
 * @brief Support for Blast Tabular file formats (with and without headers)
 *
 * This tag is part of the
 *
 * There are three blast format related tags in SeqAn:
 *
 * <li> @link BlastReport @endlink with the FormattedFile output specialization @link BlastReportOut @endlink</li>
 * <li> @link BlastTabular @endlink with the FormattedFile output and input specializations
 * @link BlastTabularOut @endlink and @link BlastTabularIn @endlink</li>
 * <li> @link BlastTabularLL @endlink which provides light-weight, but very basic tabular IO </li>
 *
 * This is the second tag, it offers <b>high-level</b> support for reading and writing NCBI Blast compatible
 * <b>tabular</b> files, both with and without headers. These are the formats that are available in legacy Blast
 * (<tt>blastall</tt> executable) with the parameters <tt>-m 8</tt> and <tt>-m 9</tt> (with headers)
 * and in BLAST+ (<tt>blastx</tt>, <tt>blastn</tt>...) with
 * the parameters <tt>-outfmt 6</tt> and <tt>-outfmt 7</tt> respectively.
 *
 * Please consult the documentation for @link BlastIOContext @endlink to understand
 * the different options you have with this format.
 *
 * For very basic tabular IO there is the third tag, @link BlastTabularLL @endlink. *
 *
 * The reference Blast implementation used for developing the SeqAn support is NCBI Blast+ 2.2.26 and
 * NCBI Blast 2.2.26 for the legacy support.
 *
 * @section Input
 *
 * If you are interested in iteratively reading all records from a blast tabular file, you just need to call
 * readRecord:
 *
 * <ul>
 * <li>while not @link StreamConcept#atEnd @endlink </li>
 * <ul>
 * <li> BlastTabular#@link BlastTabular#readRecord @endlink</li>
 * </ul>
 * </ul>
 *
 * Alternatively to this spec's interface you may also call BlastTabular<b>In</b>#@link BlastTabularIn#readRecord @endlink
 * which uses the @link FormattedFile @endlink's interface but is the same otherwise. See @link BlastTabularIn @endlink
 * for a full code example.
 *
 * While @link BlastTabular#readHeader @endlink and @link BlastTabular#readFooter @endlink are defined, they are no-ops
 * and all information can be extracted by calling readRecord.
 *

 * @section Output
 *
 * High-level file writing consists of three different calls:
 * <ul>
 * <li> BlastTabular#@link BlastTabular#writeHeader @endlink </li>
 * <li> BlastTabular#@link BlastTabular#writeRecord @endlink repeated up to n times</li>
 * <li> BlastTabular#@link BlastTabular#writeFooter @endlink </li>
 * </ul>
 *
 * You can also use the FormattedFile interface which has the same functions, but only requires one parameter:
 * <ul>
 * <li> BlastTabularOut#@link BlastTabularOut#writeHeader @endlink </li>
 * <li> BlastTabularOut#@link BlastTabularOut#writeRecord @endlink repeated up to n times</li>
 * <li> BlastTabularOut#@link BlastTabularOut#writeFooter @endlink </li>
 * </ul>
 *
 * See @link BlastTabularOut @endlink for a full code example.
 *
 * Strictly speaking the writeHeader call is not required for BlastTabular, but for consistency with other (blast)
 * formats it is recommended. The example in @link BlastTabularOut @endlink illustrates that.
 *
 */
struct BlastTabular_;
typedef Tag<BlastTabular_> BlastTabular;

/*!
 * @enum BlastTabularSpec
 * @headerfile <seqan/blast.h>
 * @signature enum class BlastTabularSpec : uint8_t { ... };
 * @brief Spec for @link BlastIOContext @endlink
 *
 * @val BlastTabularSpec BlastTabularSpec::NO_HEADER
 * @brief Tabular format without header
 *
 * @val BlastTabularSpec BlastTabularSpec::HEADER
 * @brief Tabular format with header
 *
 * @val BlastTabularSpec BlastTabularSpec::UNKNOWN
 * @brief not defined or not known
 */
enum class BlastTabularSpec : uint8_t
{
    NO_HEADER = 0,
    HEADER = 1,
    UNKNOWN = 255,
};

/*!
 * @defgroup BlastTabularSpecTagGroup BlastTabularSpec Integral Constants
 * @brief Integral constants for @link BlastTabularSpec @endlink
 * TODO(h4nn3s): does this need to be publicly visible?
 *
 * @tag BlastTabularSpecTagGroup#BlastTabularSpecTag
 * @brief Templated typedef for @link BlastTabularSpec @endlink
 * @signature template <BlastTabularSpec h>
 * using BlastTabularSpecTag = std::integral_constant<BlastTabularSpec, h>;
 *
 */
template <BlastTabularSpec h>
using BlastTabularSpecTag = std::integral_constant<BlastTabularSpec, h>;

/*!
 * @tag BlastTabularSpecTagGroup#BlastTabularSpecTagHeader
 * @brief Integral constant for @link BlastTabularSpec::HEADER @endlink
 *
 * @tag BlastTabularSpecTagGroup#BlastTabularSpecTagNoHeader
 * @brief Integral constant for @link BlastTabularSpec::NO_HEADER @endlink
 *
 * @tag BlastTabularSpecTagGroup#BlastTabularSpecTagUnkown
 * @brief Integral constant for @link BlastTabularSpec::UNKNOWN @endlink
 */
typedef BlastTabularSpecTag<BlastTabularSpec::HEADER>     BlastTabularSpecTagHeader;
typedef BlastTabularSpecTag<BlastTabularSpec::NO_HEADER>  BlastTabularSpecTagNoHeader;
typedef BlastTabularSpecTag<BlastTabularSpec::UNKNOWN>    BlastTabularSpecTagUnkown;

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<BlastTabular, T> :
    public MagicHeader<Nothing, T> {};

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<BlastTabular, T>
{
    static constexpr char const * VALUE[5] =
    {
        ".blast",
        ".m8",
        ".bm8",
        ".m9",
        ".bm9"
    };
};

template <typename T>
constexpr char const * FileExtensions<BlastTabular, T>::VALUE[5];

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TContext, typename TDirection, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<BlastTabular, TDirection, TContext>, TStorageSpec>
{
    typedef TContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<BlastTabular, TDirection, TSpec> >
{
    typedef BlastTabular Type;
};

// ----------------------------------------------------------------------------
// Class BlastMatchField
// ----------------------------------------------------------------------------

/*!
 * @class BlastMatchField BlastMatchField
 * @brief A "meta" datastructure that contains information about members of @link BlastMatch @endlinkes
 * @headerfile seqan/blast.h
 *
 * @signature template <typename TVoidSpec = void> struct BlastMatchField;
 * @tparam TVoidSpec An extra spec to prevent global inclusion of statics members (you can safely ignore this)
 *
 * This data structure conveniently gives access to all possible fields used in
 * BLAST-compatabile tabular output formats. @link BlastMatchField::Enum @endlink is needed to
 * specify a custom field composition for a @link BlastIOContext @endlink.
 *
 * The member variables offer the correct labels for the tabular formats' I/O and strings
 * to interact with the user on the command line.
 *
 * Please note that for the legacyFormat (@link BlastIOContext::legacyFormat @endlink) specifying or reading custom
 * fields is not supported and the columnsLabels will always be printed as
 * @link BlastMatchField::legacyColumnLabels @endlink.
 *
 * <h3>Table overview</h3>
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
 * More fields will likely be implemented in the future.
 */


template <typename TVoidSpec = void>
struct BlastMatchField
{
    /*!
     * @enum BlastMatchField::Enum
     * @headerfile seqan/blast.h
     * @signature enum class BlastMatchField<TVoidSpec>::Enum : uint8_t { ... };
     * @brief A strongly typed enum mapping all fields supported by NCBI Blast to an integer
     *
     * The available values are visible in detailed description of @link BlastMatchField @endlink.
     */
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

    /*!
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

    /*!
     * @var static_constexpr_char_const_*_const BlastMatchField::optionLabels[]
     * @brief An array of CStrings representing the command line parameter name of each field
     */
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

    /*!
     * @var static_constexpr_char_const_*_const BlastMatchField::legacyColumnLabels
     * @brief A single CString representing the <b>column labels</b> of the @link BlastIOContext::legacyFormat @endlink.
     */
    static constexpr char const * const legacyColumnLabels =
    {
        "Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s."
         " end, e-value, bit score"
    };

    /*!
     * @var static_constexpr_char_const_*_const BlastMatchField::columnLabels[]
     * @brief An array of CStrings representing the <b>column label</b> of each possible field; for the
     * @link BlastIOContext::legacyFormat @endlink, use @link BlastMatchField::legacyColumnLabels @endlink instead.
     */
    static constexpr char const * const columnLabels [] =
    {
        "query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. "
        "end, evalue, bit score",
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

    /*!
     * @var static_constexpr_char_const_*_const BlastMatchField::descriptions[]
     * @brief An array of CStrings representing the human-readable descriptions of each field
     */
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

    /*!
     * @var static_constexpr_bool_const BlastMatchField::implemented[]
     * @brief An array of bools revealing whether the Blast I/O module supports printing this field
     */
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
constexpr char const * const BlastMatchField<TVoidSpec>::optionLabels[45];

template <typename TVoidSpec>
constexpr char const * const BlastMatchField<TVoidSpec>::legacyColumnLabels;

template <typename TVoidSpec>
constexpr char const * const BlastMatchField<TVoidSpec>::columnLabels[45];

template <typename TVoidSpec>
constexpr char const * const BlastMatchField<TVoidSpec>::descriptions[45];

template <typename TVoidSpec>
constexpr bool const BlastMatchField<TVoidSpec>::implemented[45];

template <typename TVoidSpec>
constexpr std::array<typename BlastMatchField<TVoidSpec>::Enum const, 12> BlastMatchField<TVoidSpec>::defaults;

// ============================================================================
// Metafunctions and global const-expressions
// ============================================================================

}

#endif
