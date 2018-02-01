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
// This file contains the BlastIOContext's code
// ==========================================================================

#ifndef SEQAN_BLAST_BLAST_IO_CONTEXT_H__
#define SEQAN_BLAST_BLAST_IO_CONTEXT_H__

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

template <typename TScore>
struct BlastScoringScheme;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Mfn BlastIOContextStringType_
// ----------------------------------------------------------------------------

template <typename TContext>
struct BlastIOContextStringType_
{
    typedef std::string Type;
};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class BlastIOContext
// ----------------------------------------------------------------------------

/*!
 * @class BlastIOContext
 * @headerfile <seqan/blast.h>
 * @signature template <typename TScore_ = Blosum62,
 * BlastProgram p = BlastProgram::UNKNOWN, BlastTabularSpec h = BlastTabularSpec::UNKNOWN>
 * struct BlastIOContext { ... };
 *
 * @brief An object that holds file global information and buffers for BlastIO
 *
 * @tparam TScore   Type of the @link Score @endlink object used.
 * @tparam p        @link BlastProgram @endlink as compile-time parameter.
 * @tparam h        @link BlastTabularSpec @endlink as compile-time parameter.
 *
 * This is a part of the Blast formatted files. Before writing, some of the context's members should be set; after
 * reading it will contain
 * all information from the file that did not belong to a @link BlastRecord @endlink, e.g. the name of the database.
 * It also contains buffers for internal use.
 *
 * You should re-use this object (i.e. only create it once for
 * every file that you read/write). And you don't need to and should not clear()
 * this, except when restarting IO on a different file.
 *
 * To speed-up file writing slightly you can set the value template parameters <tt>p</tt> and/or <tt>h</tt> to something
 * other than ::DYNAMIC at compile-time (e.g. if you know that you will be printing only BLASTX), but then you won't
 * be able to modify these values at run-time. For file reading this is also possible, but usually the
 * added flexibility of automatically detecting these values is prefferable.
 *
 * If not explicitly stated otherwise, the member variables are <i>out-parameters</i> of <tt>readHeader()</tt>,
 * <tt>readRecord()</tt> and <tt>readFooter()</tt>, i.e. they are set by these functions; and they are
 * <i>in-parameters</i> to  <tt>writeHeader()</tt>, <tt>writeRecord()</tt> and <tt>writeFooter()</tt>, i.e. they
 * influence these functions' output.
 *
 * See @link BlastTabularFileOut @endlink and @link BlastReportFileOut @endlink for more complete examples of usage.
 */

template <typename TScore_ = Blosum62,
          BlastProgram p = BlastProgram::DYNAMIC,
          BlastTabularSpec h = BlastTabularSpec::DYNAMIC>
struct BlastIOContext
{
    typedef TScore_ TScore;
    typedef typename BlastIOContextStringType_<BlastIOContext>::Type TString;

    /*!
     * @var BlastProgramSelector BlastIOContext::blastProgram;
     * @brief The @link BlastProgram @endlink.
     *
     * @section Remarks
     *
     * Behaves exactly like an enum of type @link BlastProgram @endlink, unless the second template parameter was
     * specified to make this a compile-time constant. See @link BlastProgramSelector @endlink for more information.
     */
    BlastProgramSelector<p> blastProgram;

    /*!
     * @var BlastTabularSpecSelector BlastIOContext::tabularSpec;
     * @brief The @link BlastTabularSpec @endlink.
     *
     * @section Remarks
     *
     * Behaves exactly like an enum of type @link BlastTabularSpec @endlink, unless the third template parameter was
     * specified to make this a compile-time constant. See @link BlastTabularSpecSelector @endlink for more information.
     */
    BlastTabularSpecSelector<h> tabularSpec;

    /*!
     * @var BlastScoringScheme<TScore> BlastIOContext::scoringScheme;
     * @brief The @link BlastScoringScheme @endlink.
     */
    BlastScoringScheme<TScore> scoringScheme;

    /*!
     * @var TString BlastIOContext::versionString;
     * @brief The blast version string.
     *
     * @section Remarks
     *
     * Used when writing @link BlastReportFileOut @endlink and @link BlastTabularFileOut @endlink if the context's tabularSpec
     * is set to BlastTabularSpec::COMMENTS. Defaults to a version string based on the emulated
     * blast version and the current SeqAn version.
     * When reading from @link BlastTabularFileOut @endlink the corresponding line is extracted from the comment lines
     * (if present).
     */
    TString versionString;
    void _setDefaultVersionString()
    {
        clear(versionString);
        append(versionString, _programTagToString(blastProgram));
        append(versionString, " 2.2.26");
        if (!legacyFormat)
            append(versionString, "+");
        append(versionString, " [I/O Module of SeqAn-");
        append(versionString, std::to_string(SEQAN_VERSION_MAJOR));
        append(versionString, '.');
        append(versionString, std::to_string(SEQAN_VERSION_MINOR));
        append(versionString, '.');
        append(versionString, std::to_string(SEQAN_VERSION_PATCH));
        append(versionString, ", http://www.seqan.de]");
    }

    /*!
     * @var bool BlastIOContext::legacyFormat;
     * @brief Whether to use the legacy format (only @link BlastTabular @endlink).
     *
     * @section Remarks
     *
     * Setting this flag when writing to a @link BlastTabularFileOut @endlink (that has BlastTabularSpec::COMMENTS set)
     * will result in the legacy version of the comments being written. In the legacy format the mismatches column
     * also includes all gaps in addition to mismatches.
     * Note that many other features like custom fields are not supported in this format.
     *
     * When reading @link BlastTabularFileOut @endlink this flag will automatically be set based on the comments (if a
     * they exist).
     */
    bool legacyFormat = false;

    /*!
     * @var TString BlastIOContext::dbName;
     * @brief Name of the dabase or path to the file.
     */
    TString         dbName;

    /*!
     * @var uint64_t BlastIOContext::dbTotalLength;
     * @brief Summed up sequence length of the database.
     */
    uint64_t        dbTotalLength = 0u;

    /*!
     * @var uint64_t BlastIOContext::dbNumberOfSeqs;
     * @brief Number of sequences in the database.
     */
    uint64_t        dbNumberOfSeqs = 0u;

    /*!
     * @var StringSet<TString> BlastIOContext::otherLines;
     * @brief A StringSet that will contain all comment lines that
     * could not be interpreted in another way (only @link BlastTabularFileIn @endlink).
     */
    StringSet<TString, Owner<ConcatDirect<>>> otherLines;

    /*!
     * @var std::vector<BlastMatchField::Enum> BlastIOContext::fields;
     * @brief The fields (types of columns) in @link BlastTabular @endlink-formats.
     *
     * @section Remarks
     *
     * This is an <i>out-parameter</i> for:
     * <li> @link BlastTabularFileIn#readRecord @endlink iff tabularSpec == COMMENTS (otherwise it can't be deduced)</li>
     *
     * This is an <i>in-parameter</i> for:
     * <li> @link BlastTabularFileIn#readRecord @endlink if tabularSpec != COMMENTS (specified fields will be expected)</li>
     * <li> @link BlastTabularFileOut#writeRecord @endlink (specified fields will written)
     *
     * Setting @link BlastIOContext::ignoreFieldsInComments @endlink will make this variable be an <i>in-parameter</i> for
     * the first case, as well. This variable is ignored in the legacy formats and for non-tabular formats.
     */
    std::vector<typename BlastMatchField<>::Enum> fields = { BlastMatchField<>::Enum::STD };

    /*!
     * @var StringSet<TString> BlastIOContext::fieldsAsStrings;
     * @brief The fields (types of columns) in @link BlastTabular @endlink-formats, but as uninterpreted strings.
     *
     * @section Remarks
     *
     * Useful when the comment lines do not conform to standards and you want to extract the verbatim column labels or
     * if you wish to print non-standard column labels (which you shouldn't!).
     */
    StringSet<TString, Owner<ConcatDirect<>>> fieldsAsStrings;

    /*!
     * @var bool BlastIOContext::ignoreFieldsInComments;
     * @brief Use fields as in-parameter for readRecord as well (only @link BlastTabularFileIn @endlink).
     *
     * @section Remarks
     *
     * See @link BlastTabularFileIn#readRecord @endlink. Use this when the comment lines do not
     * conform to standards (and the fields can't be read), but you know that
     * the matches are in the given, e.g. default format.
     */
    bool ignoreFieldsInComments = false;

    /*!
     * @var StringSet<TString> BlastIOContext::conformancyErrors;
     * @brief Holds non fatal error messages when reading from @link BlastTabularFileIn @endlink.
     *
     * @section Remarks
     *
     * After doing a @link BlastTabularFileIn#readRecord @endlink this will indicate whether the
     * comment lines contained non-fatal parse errors, usually the result
     * of a file written by a sloppy blast implementation or possibly a bug in SeqAn.
     * An empty StringSet indicates that all is good.
     */
    StringSet<TString, Owner<ConcatDirect<>>> conformancyErrors;

    // ------- CACHES, BUFFERS and INTERNALS --------- //

    // counted internally for TabularFooter
    uint64_t _numberOfRecords = 0;

    // cache for length adjustments in blast statistics
    std::unordered_map<uint64_t, uint64_t> _cachedLengthAdjustments;

    // io-buffers
    TString _lineBuffer; // holds the current line
    TString _stringBuffer;
    StringSet<TString, Owner<ConcatDirect<>>> _setBuffer1;
    StringSet<TString, Owner<ConcatDirect<>>> _setBuffer2;
    BlastMatch<> bufMatch;
    BlastRecord<> bufRecord;
};

}

#endif
