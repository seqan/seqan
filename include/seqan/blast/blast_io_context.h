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
// This file contains the BlastIOContext's code
// ==========================================================================

#ifndef SEQAN_BLAST_BLAST_IO_CONTEXT_H__
#define SEQAN_BLAST_BLAST_IO_CONTEXT_H__

namespace seqan
{

//forward declare
template <typename TScore, typename TSpec = void>
struct BlastScoringAdapter;

/*!
 * @class BlastIOContext
 * @headerfile <seqan/blast.h>
 * @signature template <typename TScore_ = Blosum62, typename TString = std::string,
 *           BlastProgram p = BlastProgram::UNKNOWN, BlastTabularSpec h = BlastTabularSpec::UNKNOWN>
 * struct BlastIOContext { ... };
 *
 * @brief An object that holds file global information and buffers for BlastIO
 *
 * @tparam TScore   Type of the @link Score @endlink object used.
 * @tparam TString  Type of the @link StringConcept @endlink used for members and buffers.
 * @tparam p        @link BlastProgram @endlink as compile-time parameter.
 * @tparam h        @link BlastTabularSpec @endlink as compile-time parameter.
 *
 * This needs to be passed to most read*(), skip*() and write*() functions as
 * a parameter. Before writing, some of the context's members should be set; after reading it will contain
 * all information from the file that did not belong to a @link BlastRecord @endlink, e.g. the name of the database.
 * It also contains buffers for internal use.
 *
 * You should re-use this object (i.e. only create it once for
 * every file that you read/write). And you don't need to and should not clear()
 * this, except when restarting IO on a different file.
 *
 * To speed-up file writing slightly you can set the value template parameters <tt>p</tt> and/or <tt>h</tt> to something
 * other than ::UNKNOWN at compile-time (e.g. if you know that you will be printing only BLASTX), but then you won't
 * be able to modify these values with @link BlastIOContext#setBlastProgram @endlink and
 * @link BlastIOContext#setBlastTabularSpec @endlink at run-time. For file reading this also possible, but usually the
 * added flexibility of automatically detecting these values is prefferable.
 *
 * @section Example
 *
 * Here as an example of which members to set on a context, before using it for Output:
 * @code{.cpp}
 * typedef BlastIOContext<Blosum62> TContext;
 * BlastReportOut<TContext> outfile("/tmp/output.blast");
 *
 * Blosum62 scheme;
 * setScoreGapOpen(scheme, -11);
 * setScoreGapExtend(scheme, -1);
 *
 * // upon assigning, this is converted to SeqAn's scoring behaviour
 * setBlastScoringScheme(context(outfile), scheme);
 *
 * // protein vs protein search is BLASTP
 * setBlastProgram(context(outfile), BlastProgram::BLASTP);
 *
 * // set the database properties in the context
 * context(outfile).dbName = "The Foo Database";
 * context(outfile).dbTotalLength = length(concat(subjects));
 * context(outfile).dbNumberOfSeqs = length(subjects);
 * @endcode
 *
 * See @link BlastTabularOut @endlink and @link BlastReportOut @endlink for more complete examples of usage.
 */

template <typename TScore_ = Blosum62,
          typename TString_ = std::string,
          BlastProgram _p = BlastProgram::UNKNOWN,
          BlastTabularSpec _h = BlastTabularSpec::UNKNOWN>
struct BlastIOContext
{
    typedef TScore_ TScore;
    typedef TString_ TString;
    static constexpr BlastProgram p = _p;
    static constexpr BlastTabularSpec h = _h;

    // the blastProgram as compile time parameter
    typedef BlastProgramTag<p> TBlastProgram; // compile-time parameter
    // if the upper is set to UNKNOWN, than this run-time variable is consulted instead:
    BlastProgram blastProgram = BlastProgram::UNKNOWN;

    // the BlastTabularSpec as compile time parameter
    typedef BlastTabularSpecTag<h> TBlastTabularSpec; // compile-time parameter
    // if the upper is set to UNKNOWN, than this run-time variable is consulted instead:
    BlastTabularSpec tabularSpec = BlastTabularSpec::HEADER;

    /*!
     * @var BlastScoringAdapter<TScore> BlastIOContext::scoringAdapter;
     * @brief A @link BlastScoringAdapter @endlink, you should never have to access this directly.
     */
    BlastScoringAdapter<TScore> scoringAdapter;

    /*!
     * @var TString BlastIOContext::versionString;
     * @brief the blast version string
     *
     * Used when writing @link BlastReportOut @endlink and @link BlastTabularOut @endlink if the context's tabularSpec
     * is set to BlastTabularSpec::HEADER. Defaults to a version string based on the emulated
     * blast version and the current SeqAn version.
     * When reading from @link BlastTabularOut @endlink the corresponding line is extracted from the header
     * (if present).
     */
    TString versionString;
    void _setDefaultVersionString()
    {
        clear(versionString);
        append(versionString, _programTagToString(blastProgram, TBlastProgram()));
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
     * @brief whether to use the legacy format (only @link BlastTabular @endlink).
     *
     * Setting this flag when writing to a @link BlastTabularOut @endlink (that has BlastTabularSpec::HEADER set) will
     * result in the legacy header being written. This is the slightly different header used by C-only versions of blast
     * (<tt>blastall</tt>-binary). In the legacy format the mismatches column also includes all gaps in addition to
     * mismatches. Note that many other features like custom fields are not supported in this format.
     *
     * When reading @link BlastTabularOut @endlink this flag will automatically be set based on the header (if a
     * header exists).
     */
    bool legacyFormat = false;

    /*!
     * @var TDbName BlastIOContext::dbName;
     * @brief Name of the dabase or path to the file
     */
    TString         dbName;

    /*!
     * @var uint64_t BlastIOContext::dbTotalLength;
     * @brief summed up sequence length of the database
     */
    uint64_t        dbTotalLength = 0u;

    /*!
     * @var uint64_t BlastIOContext::dbNumberOfSeqs;
     * @brief number of sequences in the database
     */
    uint64_t        dbNumberOfSeqs = 0u;

    /*!
     * @var StringSet<TString> BlastIOContext::otherLines;
     * @brief a StringSet that will contain all comment or header lines that
     * could not be interpreted in another way (only @link BlastTabularIn @endlink).
     */
    StringSet<TString, Owner<ConcatDirect<>>> otherLines;

    /*!
     * @var std::vector<BlastMatchField::Enum> BlastIOContext::fields;
     * @brief the fields (types of columns) in @link BlastTabular @endlink-formats.
     *
     * Is an in-parameter to writeRecord, writeMatch, writeRecordHeader and readMatch. In
     * the latter case it signifies the expected fields.
     * Is an out-parameter
     * of readRecord and readRecordHeader where it returns the fields specified in
     * the header.
     */
    std::vector<typename BlastMatchField<>::Enum> fields { { BlastMatchField<>::Enum::STD } };

    /*!
     * @var StringSet<TString> BlastIOContext::fieldsAsStrings;
     * @brief the fields (types of columns) in @link BlastTabular @endlink-formats, but as uninterpreted strings.
     *
     * Useful when the header does not conform to standards and you want to extract the verbatim column labels or if
     * you wish to print non-standard column labels (which you don't!).
     */
     StringSet<TString, Owner<ConcatDirect<>>> fieldsAsStrings;

    /*!
     * @var bool BlastIOContext::ignoreFieldsInHeader;
     * @brief used fields as in-parameter for readRecord as well (only @link BlastTabularIn @endlink).
     *
     * When doing @link BlastTabular#readRecord @endlink, the
     * @link BlastIOContext::fields @endlink member is used as in-parameter to
     * readRecordHeader() and as out-parameter to readMatch(); setting this bool
     * deactivates the first behaviour. Use this when the header does not
     * conform to standards (and the fields can't be read), but you know that
     * the matches are in the given, e.g. default format.
     */
    bool ignoreFieldsInHeader = false;

    /*!
     * @var StringSet<TString> BlastIOContext::conformancyErrors;
     * @brief holds non fatal error messages when reading from @link BlastTabularIn @endlink.
     *
     * After doing a @link BlastTabular#readRecord @endlink this will indicate whether the
     * record header contained non-fatal parse errors, usually the result
     * of a file written by a sloppy blast implementation or possibly a bug in SeqAn.
     * An empty StringSet indicates that all is good.
     */
    StringSet<TString, Owner<ConcatDirect<>>> conformancyErrors;

    // ------- CACHES, BUFFERS and INTERNALS --------- //

    // TODO prefix the below with _ ?
    // counted internally for TabularFooter
    uint64_t numberOfRecords = 0;

    // cache for length adjustments in blast statistics
    std::unordered_map<uint64_t, uint64_t> cachedLengthAdjustments;

    // needed for readRecord of TABULAR
    TString lastId;

    // io-buffers
    TString buffer1;
    TString buffer2;
    StringSet<TString, Owner<ConcatDirect<>>> buffers1;
    StringSet<TString, Owner<ConcatDirect<>>> buffers2;
    BlastMatch<> bufMatch;
    BlastRecord<> bufRecord;
};

/*!
* @fn BlastIOContext#getBlastProgram
* @signature BlastProgram getBlastProgram(context);
* @param[in] context  The @link BlastIOContext @endlink.
* @brief return the @link BlastProgram @endlink set in the context.
*
* Developer note: this function is constexpr when the
* context's compile-time parameter p != BlastProgram::UNKNOWN. Otherwise it is plain inline.
*/

template <typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
constexpr BlastProgram
getBlastProgram(BlastIOContext<TScore, TString, p, h> const &)
{
    return p;
}

template <typename TScore,
          typename TString,
          BlastTabularSpec h>
inline BlastProgram
getBlastProgram(BlastIOContext<TScore, TString, BlastProgram::UNKNOWN, h> const & context)
{
    return context.blastProgram;
}

/*!
* @fn BlastIOContext#setBlastProgram
* @signature void setBlastProgram(context, blastProgram);
* @param[in,out] context        The @link BlastIOContext @endlink.
* @param[in]     blastProgram   The @link BlastProgram @endlink you wish to set.
* @brief set the program type at run-time
*
* Note that this function will bail out, if the blastProgram was already set as a compile-time argument to the
* context (i.e. it is not euqal to BlastProgram::UNKNOWN) and the value you wish to set is not the same as the one set
* at compile time.
*/

template <typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
setBlastProgram(BlastIOContext<TScore, TString, p, h> &, BlastProgram const _p)
{
    if (p != _p)
        SEQAN_FAIL("ERROR: Tried to set blastProgram on context, but was already defined at compile time (and set to a "
        "different value!");
}

template <typename TScore,
          typename TString,
          BlastTabularSpec h>
inline void
setBlastProgram(BlastIOContext<TScore, TString, BlastProgram::UNKNOWN, h> & context, BlastProgram const p)
{
    context.blastProgram = p;
}

/*!
* @fn BlastIOContext#getBlastTabularSpec
* @signature BlastTabularSpec getBlastTabularSpec(context);
* @param[in] context  The @link BlastIOContext @endlink.
* @brief return the @link BlastTabularSpec @endlink set in the context.
*
* Developer note: this function is constexpr when the
* context's compile-time parameter p != BlastTabularSpec::UNKNOWN. Otherwise it is plain inline.
*/

template <typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
constexpr BlastTabularSpec
getBlastTabularSpec(BlastIOContext<TScore, TString, p, h> const &)
{
    return h;
}

template <typename TScore,
          typename TString,
          BlastProgram p>
inline BlastTabularSpec
getBlastTabularSpec(BlastIOContext<TScore, TString, p, BlastTabularSpec::UNKNOWN> const & context)
{
    return context.tabularSpec;
}

/*!
* @fn BlastIOContext#setBlastTabularSpec
* @signature void setBlastTabularSpec(context, tabularSpec);
* @param[in,out] context        The @link BlastIOContext @endlink.
* @param[in]     tabularSpec    The @link BlastTabularSpec @endlink you wish to set.
* @brief set the tabular spec at run-time
*
* Note that this function will bail out, if the tabularSpec was already set as a compile-time argument to the
* context (i.e. it is not euqal to BlastTabularSpec::UNKNOWN) and the value you wish to set is not the same as the one set
* at compile time.
*/

template <typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
setBlastTabularSpec(BlastIOContext<TScore, TString, p, h> &, BlastTabularSpec const _h)
{
    if (h != _h)
        SEQAN_FAIL("ERROR: Tried to set tabularSpec on context, but was already defined at compile time (and set to a "
        "different value!");
}

template <typename TScore,
          typename TString,
          BlastProgram p>
inline void
setBlastTabularSpec(BlastIOContext<TScore, TString, p, BlastTabularSpec::UNKNOWN> & context,
                    BlastTabularSpec const h)
{
    context.tabularSpec = h;
}

/*!
* @fn BlastIOContext#getScoringScheme
* @signature ScoringScheme const & getScoringScheme(context);
* @param[in] context  The @link BlastIOContext @endlink.
* @brief get the @link Score @endlink object of the context.
* @return ScoringScheme reference to the scheme used in the context.
*/

template <typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline TScore const &
getScoringScheme(BlastIOContext<TScore, TString, p, h> const & context)
{
    return context.scoringAdapter.scheme;
}

/*!
 * @fn BlastIOContext#setScoringScheme
 * @signature void setScoringScheme(context, scoringScheme);
 * @param[in,out] context        The @link BlastIOContext @endlink.
 * @param[in]     scoringScheme  The new @link Score @endlink object.
 * @brief set the @link Score @endlink object of the context.
 *
 * @section Remarks
 *
 * It is important to note that gap-costs are computed differently in SeqAn and in BLAST, see
 * @link BlastIOContext#setBlastScoringScheme @endlink for more details.
 *
 * After setting the scoringScheme, the scoringAdpater will be updated. If this fails (e.g. no statistical parameters
 * for the scoring scheme are available) this function will call SEQAN_FAIL. If you don't want this behaviour and you
 * absolutely know what you are doing, you can operate directly on the
 * @link BlastIOContext::scoringAdapter @endlink-member.
*/

template <typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
setScoringScheme(BlastIOContext<TScore, TString, p, h> & context, TScore const & scoringScheme)
{
    context.scoringAdapter.scheme = scoringScheme;

    if (!_selectSet(context.scoringAdapter))
        SEQAN_FAIL("ERROR: No Karlin-Altschul parameters where available for you scoring scheme --> your scoring scheme"
                   " and/or your gap costs are not supported.");
}

/*!
* @fn BlastIOContext#getBlastScoringScheme
* @signature ScoringScheme getBlastScoringScheme(context);
* @param[in] context  The @link BlastIOContext @endlink.
* @brief get the @link Score @endlink object of the context, converting to Blast behaviour.
* @return ScoringScheme a copy of the scoringScheme of the context, converted to Blast's gap notation, see
* @link BlastIOContext#setBlastScoringScheme @endlink
*/

template <typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline TScore
getBlastScoringScheme(BlastIOContext<TScore, TString, p, h> const & context)
{
    TScore newscore(context.scoringAdapter.scheme);
    seqanScoringScheme2blastScoringScheme(newscore);
    return newscore;
}

/*!
 * @fn BlastIOContext#setBlastScoringScheme
 * @signature void setBlastScoringScheme(context, scoringScheme);
 * @param[in,out] context        The @link BlastIOContext @endlink.
 * @param[in]     scoringScheme  The new @link Score @endlink object.
 * @brief set the @link Score @endlink object of the context, converting from Blast behaviour.
 *
 * @section Remarks
 *
 * It is important to note that gap-costs are computed differently in SeqAn and in BLAST.
 * Blast (and many other tools) compute scores of a stretch of gaps as
 * <tt>s = gO + n * gE</tt>
 * where gO is the gapOpen score, gE is the gap extend score and n ist the
 * total number of gap characters.
 *
 * SeqAn, however, computes them as as <tt>s = gO + (n-1) * gE</tt>.
 *
 * The context always holds a scoring scheme in SeqAn's format, but you can use this function to pass a scoring
 * scheme in blast's notation and it will be converted automatically. If you are used to Blast's notation or provide
 * a user interface where the user enters gap-scores, you will want to use this interface.
 *
 * After setting the scoringScheme, the scoringAdpater will be updated. If this fails (e.g. no statistical parameters
 * for the scoring scheme are available) this function will call SEQAN_FAIL. If you don't want this behaviour and you
 * absolutely know what you are doing, you can operate directly on the
 * @link BlastIOContext::scoringAdapter @endlink-member.
 */

template <typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
setBlastScoringScheme(BlastIOContext<TScore, TString, p, h> & context, TScore scoringScheme)
{
    blastScoringScheme2seqanScoringScheme(scoringScheme);
    setScoringScheme(context, scoringScheme);
}

}

#endif