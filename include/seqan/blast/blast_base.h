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
// Module for handling NCBI Blast I/O and E-Value computation
// ==========================================================================

#ifndef SEQAN_EXTRAS_BLAST_BLAST_BASE_H_
#define SEQAN_EXTRAS_BLAST_BLAST_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @enum BlastFormatFile
 * @brief Enum with BLAST file format specs
 * @signature enum class BlastFormatFile : uint8_t { ... };
 *
 * @headerfile seqan/blast.h
 *
 * The following table lists the available BlastFormatFile values. The
 * corresponding command line arguments of legacy BLAST (<tt>blastall</tt>
 * binary) and BLAST+ (individual binaries for <tt>blastn</tt>, <tt>blastx</tt>
 * etc) are given and the corresponding version after which this implementation
 * was modelled. It should work with newer versions, but as the format is
 * sometimes changed without notice, this might not always be the case.
 *
 * Only the most important formats are implemented, see the respective input and
 * output columns.
 *
 *
 * @htmlonly
 * <span style="font-size:90%"><table>
 * <tr><th> BlastFormatFile::            </th><th>   #  </th><th> blast* 2.2.26+  </th><th>  blastall-2.2.26 </th><th>  output   </th><th>  input     </th><th> description                                       </th></tr>
 * <tr><td> PAIRWISE                     </td><td>   0  </td><td>-outfmt 0        </td><td>   -m 0           </td><td>   &#9745; </td><td>  &#9744;   </td><td> pairwise                                          </td></tr>
 * <tr><td> MASTER_SLAVE_IDENT           </td><td>   1  </td><td>-outfmt 1        </td><td>   -m 1           </td><td>   &#9744; </td><td>  &#9744;   </td><td> query-anchored showing identities                 </td></tr>
 * <tr><td> MASTER_SLAVE_NO_IDENT        </td><td>   2  </td><td>-outfmt 2        </td><td>   -m 2           </td><td>   &#9744; </td><td>  &#9744;   </td><td> query-anchored no identities                      </td></tr>
 * <tr><td> FLAT_MASTER_SLAVE_IDENT      </td><td>   3  </td><td>-outfmt 3        </td><td>   -m 3           </td><td>   &#9744; </td><td>  &#9744;   </td><td> flat query-anchored, show identities              </td></tr>
 * <tr><td> FLAT_MASTER_SLAVE_NO_IDENT   </td><td>   4  </td><td>-outfmt 4        </td><td>   -m 4           </td><td>   &#9744; </td><td>  &#9744;   </td><td> flat query-anchored, no identities                </td></tr>
 * <tr><td> XML                          </td><td>   5  </td><td>-outfmt 5        </td><td>   -m 7           </td><td>   &#9744; </td><td>  &#9744;   </td><td> XML Blast output                                  </td></tr>
 * <tr><td> TABULAR                      </td><td>   6  </td><td>-outfmt 6        </td><td>   -m 8           </td><td>   &#9745; </td><td>  &#9745;   </td><td> tabular                                           </td></tr>
 * <tr><td> TABULAR_WITH_HEADER          </td><td>   7  </td><td>-outfmt 7        </td><td>   -m 9           </td><td>   &#9745; </td><td>  &#9745;   </td><td> tabular with comment lines                        </td></tr>
 * <tr><td> TEXT_ASN1                    </td><td>   8  </td><td>-outfmt 8        </td><td>   -m 10          </td><td>   &#9744; </td><td>  &#9744;   </td><td> Text ASN.1                                        </td></tr>
 * <tr><td> BIN_ASN1                     </td><td>   9  </td><td>-outfmt 9        </td><td>   -m 11          </td><td>   &#9744; </td><td>  &#9744;   </td><td> Binary ASN.1                                      </td></tr>
 * <tr><td> CSV                          </td><td>  10  </td><td>-outfmt 10       </td><td>                  </td><td>   &#9744; </td><td>  &#9744;   </td><td> Comma-separated values                            </td></tr>
 * <tr><td> BLAST_AR_ASN1                </td><td>  11  </td><td>-outfmt 11       </td><td>                  </td><td>   &#9744; </td><td>  &#9744;   </td><td> BLAST archive format (ASN.1)                      </td></tr>
 * <tr><td> MASTER_SLAVE_BLUNT_ENDS      </td><td>  12  </td><td>                 </td><td>   -m 5           </td><td>   &#9744; </td><td>  &#9744;   </td><td> query-anchored no identities and blunt ends       </td></tr>
 * <tr><td> FLAT_MASTER_SLAVE_BLUNT_ENDS </td><td>  13  </td><td>                 </td><td>   -m 6           </td><td>   &#9744; </td><td>  &#9744;   </td><td> flat query-anchored, no identities and blunt ends </td></tr>
 * </table></span>
 * @endhtmlonly
 */

/*
 * SeqAn can currently write PAIRWISE, TABULAR and TABULAR_WITH_HEADER. It can
 * read TABULAR and TABULAR_WITH_HEADER.
 *
 * @val BlastFormatFile BlastFormatFile::PAIRWISE;
 * @brief default blast output (blastall -m 0 &amp; blast* -outfmt 0)
 *
 * @val BlastFormatFile BlastFormatFile::MASTER_SLAVE_IDENT;
 * @brief master-slave showing identities (blastall -m 1 &amp; blast* -outfmt 1)
 *
 * @val BlastFormatFile BlastFormatFile::MASTER_SLAVE_NO_IDENT;
 * @brief master-slave without identities (blastall -m 2 &amp; blast* -outfmt 2)
 *
 * @val BlastFormatFile BlastFormatFile::FLAT_MASTER_SLAVE_IDENT;
 * @brief flat master-slave showing identities (blastall -m 3 &amp; blast* -outfmt 3)
 *
 * @val BlastFormatFile BlastFormatFile::FLAT_MASTER_SLAVE_NO_IDENT;
 * @brief master-slave without identities, with blunt ends (blastall -m 4, not available with Blast+)
 *
 * @val BlastFormatFile BlastFormatFile::MASTER_SLAVE_BLUNT_ENDS;
 * @brief flat master-slave without identities, with blunt ends (blastall -m 5, not available with Blast+)
 *
 * @val BlastFormatFile BlastFormatFile::FLAT_MASTER_SLAVE_BLUNT_ENDS;
 * @brief flat master-slave without identities, with blunt ends (blastall -m 6, not available with Blast+)
 *
 * @val BlastFormatFile BlastFormatFile::XML;
 * @brief XML (blastall -m 7 &amp; blast* -outfmt 8)
 *
 * @val BlastFormatFile BlastFormatFile::TABULAR;
 * @brief tab-seperated (blastall -m 8 &amp; blast* -outfmt 6)
 *
 * @val BlastFormatFile BlastFormatFile::TABULAR_WITH_HEADER;
 * @brief tab-seperated with Header / comments (blastall -m 9 &amp; blast* -outfmt 7)
 *
 * @val BlastFormatFile BlastFormatFile::TEXT_ASN1;
 * @brief Abstract Syntax Notation One (blast* -outfmt 8, not available in traditional BLAST)
 *
 * @val BlastFormatFile BlastFormatFile::BIN_ASN1;
 * @brief Abstract Syntax Notation One (blast* -outfmt 9, not available in traditional BLAST)
 *
 * @val BlastFormatFile BlastFormatFile::CSV;
 * @brief comma-seperated values (blast* -outfmt 10, not available in traditional BLAST)
 *
 * @val BlastFormatFile BlastFormatFile::BLAST_AR_ASN1;
 * @brief Blast Archive Format, Abstract Syntax Notation One (blast* -outfmt 11, not available in traditional BLAST)
 */

enum class BlastFormatFile : uint8_t
{
    PAIRWISE = 0,
    MASTER_SLAVE_IDENT = 1,
    MASTER_SLAVE_NO_IDENT = 2,
    FLAT_MASTER_SLAVE_IDENT = 3,
    FLAT_MASTER_SLAVE_NO_IDENT = 4,
    XML = 5,
    TABULAR = 6,
    TABULAR_WITH_HEADER = 7,
    TEXT_ASN1 = 8,                 // only available in Generation==Blast+
    BIN_ASN1 = 9,                  // only available in Generation==Blast+
    CSV = 10,                      // only available in Generation==Blast+
    BLAST_AR_ASN1 = 11,              // only available in Generation==Blast+
    MASTER_SLAVE_BLUNT_ENDS = 12,       // only available in Generation==Blast
    FLAT_MASTER_SLAVE_BLUNT_ENDS = 13,   // only available in Generation==Blast
    INVALID_File=255
};

/*!
 * @enum BlastFormatProgram
 * @brief Enum with BLAST program spec
 * @signature enum class BlastFormatProgram : uint8_t { ... };
 *
 * @headerfile seqan/blast.h
 *
 * @val BlastFormatProgram BlastFormatProgram::BLASTN
 * @brief Nucleotide Query VS Nucleotide Subject
 *
 * @val BlastFormatProgram BlastFormatProgram::BLASTP
 * @brief Protein Query VS Protein Subject
 *
 * @val BlastFormatProgram BlastFormatProgram::BLASTX
 * @brief translated Nucleotide Query VS Protein Subject
 *
 * @val BlastFormatProgram BlastFormatProgram::TBLASTN
 * @brief Protein Query VS translated Nucleotide Subject
 *
 * @val BlastFormatProgram BlastFormatProgram::TBLASTX
 * @brief translated Nucleotide Query VS translated Nucleotide Subject
 *
 */
enum class BlastFormatProgram : uint8_t
{
    BLASTN,         //              NUCL VS             NUCL
    BLASTP,         //              PROT VS             PROT
    BLASTX,         // TRANSLATED   NUCL VS             PROT
    TBLASTN,        //              PROT VS TRANSLATED  NUCL
    TBLASTX,        // TRANSLATED   NUCL VS TRANSLATED  NUCL
    INVALID_Program=255
};

/*!
 * @enum BlastFormatGeneration
 * @brief Enum with BLAST program generation
 * @signature enum class BlastFormatGeneration : uint8_t { ... };
 *
 * @headerfile seqan/blast.h
 *
 * @val BlastFormatGeneration BlastFormatGeneration::BLAST_LEGACY
 * @brief traditional Blast, written in C ("blastall" binary); all behaviour related
 * to this Format is based on NCBI BLAST-2.2.26
 *
 * @val BlastFormatGeneration BlastFormatGeneration::BLAST_PLUS
 * @brief Blast+, written in C++; all behaviour related
 * to this Format is based on NCBI BLAST-2.2.27+
 *
 */

enum class BlastFormatGeneration : uint8_t
{
    BLAST_LEGACY,
    BLAST_PLUS,
    INVALID_Generation=255
};

/*!
 * @class BlastFormat
 * @headerfile seqan/blast.h
 * @brief Blast Format specifier
 *
 * @signature template <BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
 *            struct BlastFormat;
 *
 * @section Data structures
 *
 * The easiest way to use the BLAST module is to organize your data in the
 * structures provided, i.e. create records (@link BlastRecord @endlink)
 * that contain matches (@link BlastMatch @endlink) of one query
 * sequence. Store some general information of your database in a
 * @link BlastDbSpecs @endlink object.
 *
 * @section E-Value statistics
 *
 * This part provides e-value statistics and related functions. To use it
 * you need to convert the SeqAn scoring scheme (see
 * @link BlastScoringScheme @endlink) and create a @link BlastScoringAdapter @endlink.
 * After this you can @link BlastMatch#calcStatsAndScore @endlink and then
 * @link BlastMatch#calcBitScoreAndEValue @endlink for your
 * @link BlastMatch @endlinkes.
 *
 * @code{.cpp}
 * // the match object
 * BlastMatch<> m;
 * m.sId = "";
 * m.qId = "FOOBAR11 the best read";
 *
 * String<AminoAcid> subj = "NSPGFRKVAYAQTKPRRLCFPNGGTRRGSF";
 * String<AminoAcid> qry  =        "VAYAQ""PRKLCYP"         ;
 * m.qLength = length(qry);
 *
 * // database specs
 * BlastDbSpecs<> dbSpecs;
 * dbSpecs.dbName = "The really non-redundant database";
 * dbSpecs.dbNumberOfSeqs = 1;
 * dbSpecs.dbTotalLength = length(subj); // usually sum of subject lengths
 *
 * // the scoring scheme
 * Blosum62 scheme;
 * setScoreGapOpen(scheme, -11);
 * setScoreGapExtend(scheme, -1);
 * blastScoringScheme2seqanScoringScheme(scheme);
 * BlastScoringAdapter<Blosum62> adapter(scheme);
 *
 * // setup the align object and compute the alignment
 * resize(rows(m.align), 2);
 * assignSource(row(m.align, 0), subj);
 * assignSource(row(m.align, 1), qry);
 * localAlignment(m.align, scheme);
 *
 *
 * // count identities, mismatches a.s.o.; compute raw score
 * calcStatsAndScore(m, scheme);
 * // calculate bit score and evalue
 * calcBitScoreAndEValue(m, dbSpecs, adapter);
 * @endcode
 *
 * @section Output (high-level)
 *
 * For a subset of formats and features there is high-level output support with
 * @link BlastTabularFileOut @endlink. This is a short example of how to use it:
 *
 * @code{.cpp}
 * BlastTabularOut out("/tmp/example.blast");
 *
 * BlastRecord<> r;
 *
 * TODO
 * @endcode
 *
 *
 * When writing Blast(-like) results to a file the easiest way to proceed
 * is to organize your matches in the data structures described above.
 *
 * When writing this data, call @link BlastFormat#writeTop @endlink
 * first, than iterate over the records, calling
 * @link BlastRecord#writeRecord @endlink on each, and finishing with
 * @link BlastFormat#writeBottom @endlink .
 * Certain BlastFormats have additional more fine-grained functions, e.g.
 * @link BlastRecord#writeHeader @endlink and
 * @link BlastMatch#writeMatch @endlink .
 *
 * Should you for some reason require more fine-grained control and/or do not
 * wish to use the mentioned data structures, you may also use these alternative
 * interfaces: @link BlastFormat#writeHeader @endlink and
 * @link BlastFormat#writeMatch @endlink (only available for tabular formats).
 *
 * @tparam f    template parameter of type @link BlastFormatFile @endlink
 * @tparam p    template parameter of type @link BlastFormatProgram @endlink
 * @tparam g    template parameter of type @link BlastFormatGeneration @endlink
 */

template <BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
struct BlastFormat_
{
};

template <BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
using BlastFormat = Tag<BlastFormat_<f, p, g>>;

// ============================================================================
// Metafunctions
// ============================================================================

//TODO document or _

template <BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
constexpr BlastFormatFile
getFileType(BlastFormat<f, p, g> const & /**/)
{
    return f;
}

template <BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
constexpr BlastFormatProgram
getProgramType(BlastFormat<f, p, g> const & /**/)
{
    return p;
}

template <BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
constexpr BlastFormatGeneration
getGenerationType(BlastFormat<f, p, g> const & /**/)
{
    return g;
}

// ----------------------------------------------------------------------------
// qHasRevComp()
// ----------------------------------------------------------------------------

constexpr bool
qHasRevComp(BlastFormatProgram const p)
{
    return ((p==BlastFormatProgram::BLASTP) || (p==BlastFormatProgram::TBLASTN))
            ? false
            : true;
}

template <BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
constexpr bool
qHasRevComp(BlastFormat<f,p,g> const & /**/)
{
    return ((p==BlastFormatProgram::BLASTP) || (p==BlastFormatProgram::TBLASTN))
            ? false
            : true;
}

template <typename TFormat>
using QHasRevComp = typename std::conditional<qHasRevComp(TFormat()),
                                              True,
                                              False>::type;

// ----------------------------------------------------------------------------
// qIsTranslated()
// ----------------------------------------------------------------------------

constexpr bool
qIsTranslated(BlastFormatProgram const p)
{
    return ((p==BlastFormatProgram::BLASTX) || (p==BlastFormatProgram::TBLASTX))
            ? true
            : false;
}

template <BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
constexpr bool
qIsTranslated(BlastFormat<f,p,g> const & /**/)
{
    return ((p==BlastFormatProgram::BLASTX) || (p==BlastFormatProgram::TBLASTX))
            ? true
            : false;
}

template <typename TFormat>
using QIsTranslated = typename std::conditional<qIsTranslated(TFormat()),
                                             True,
                                             False>::type;

// ----------------------------------------------------------------------------
// qNumFrames()
// ----------------------------------------------------------------------------

template <BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
constexpr uint8_t
qNumFrames(BlastFormat<f,p,g> const & /**/)
{
    using TFormat = BlastFormat<f,p,g>;
    return (qIsTranslated(TFormat())
            ? 6
            : (qHasRevComp(TFormat())
                ? 2
                : 1));
}

// ----------------------------------------------------------------------------
// sHasRevComp()
// ----------------------------------------------------------------------------


constexpr bool
sHasRevComp(BlastFormatProgram const p)
{
    return ((p==BlastFormatProgram::TBLASTX) ||
            (p==BlastFormatProgram::TBLASTN))
            ? true
            : false;
}

template <BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
constexpr bool
sHasRevComp(BlastFormat<f,p,g> const & /**/)
{
    return ((p==BlastFormatProgram::TBLASTX) ||
            (p==BlastFormatProgram::TBLASTN))
            ? true
            : false;
}

template <typename TFormat>
using SHasRevComp = typename std::conditional<sHasRevComp(TFormat()),
                                              True,
                                              False>::type;

// ----------------------------------------------------------------------------
// sIsTranslated()
// ----------------------------------------------------------------------------

constexpr bool
sIsTranslated(BlastFormatProgram const p)
{
    return ((p==BlastFormatProgram::TBLASTX) ||
            (p==BlastFormatProgram::TBLASTN))
            ? true
            : false;
}

template <BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
constexpr bool
sIsTranslated(BlastFormat<f,p,g> const & /**/)
{
    return ((p==BlastFormatProgram::TBLASTX) ||
            (p==BlastFormatProgram::TBLASTN))
            ? true
            : false;
}

template <typename TFormat>
using SIsTranslated = typename std::conditional<sHasRevComp(TFormat()),
                                             True,
                                             False>::type;

// ----------------------------------------------------------------------------
// sNumFrames()
// ----------------------------------------------------------------------------

template <BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
constexpr uint8_t
sNumFrames(BlastFormat<f,p,g> const & /**/)
{
    using TFormat = BlastFormat<f,p,g>;
    return (sIsTranslated(TFormat())
            ? 6
            : (sHasRevComp(TFormat())
                ? 2
                : 1));
}

// ----------------------------------------------------------------------------
// getBlastProgramType()
// ----------------------------------------------------------------------------

/*!
 * @mfn getBlastProgramType
 * @headerfile seqan/blast.h
 * @brief for given query and subject alphabets, return the @link BlastFormatProgram @endlink
 * @signature  getBlastProgramType(TQueryAlph const &, TSubjAlph const &)
 *
 * @tparam          TQueryAlph  Alphabet of of the query sequences
 * @tparam          TSubjAlph   Alphabet of of the subejct sequences
 * @return          Value of @link BlastFormatProgram @endlink . For nucl VS nucl BlastFormatProgram::BLASTN is returned, although TBlastX would also be legal.
 *
 * This is a convenience function for determining the Blast-Program type
 * corresponding to input type combinations.
 *
 * This is a constexpr function.
 */

// template< typename TQueryAlph, typename TSubjAlph>
// constexpr BlastFormatProgram
// getBlastProgramType(TQueryAlph const &, TSubjAlph const &)
// {
//     return BlastFormatProgram::INVALID_Program;
// }
//
// template<typename TQueryAlph, typename TSubjAlph,
//          typename TSpec, typename TSpec2>
// constexpr BlastFormatProgram
// getBlastProgramType(String<TQueryAlph, TSpec> const &,
//                     String<TSubjAlph, TSpec2> const &)
// {
//     // needs constexpr constructors of Alphabet types
//     return getBlastProgramType(TQueryAlph(), TSubjAlph());
// }
//
// // --- DNA vs DNA ---
// // NOTE that Dna VS Dna could also be TBlastX, but BlastN is more common
// constexpr BlastFormatProgram
// getBlastProgramType(Dna const &, Dna const &)
// {
//     return BlastFormatProgram::BLASTN;
// }
//
// constexpr BlastFormatProgram
// getBlastProgramType(Dna const &, Dna5 const &)
// {
//     return BlastFormatProgram::BLASTN;
// }
//
// constexpr BlastFormatProgram
// getBlastProgramType(Dna5 const &, Dna const &)
// {
//     return BlastFormatProgram::BLASTN;
// }
//
// constexpr BlastFormatProgram
// getBlastProgramType(Dna5 const &, Dna5 const &)
// {
//     return BlastFormatProgram::BLASTN;
// }
//
// // --- Protein vs Protein ---
// constexpr BlastFormatProgram
// getBlastProgramType(AminoAcid const &, AminoAcid const &)
// {
//     return BlastFormatProgram::BLASTP;
// }
//
// // --- Dna vs Protein ---
// constexpr BlastFormatProgram
// getBlastProgramType(Dna const &, AminoAcid const &)
// {
//     return BlastFormatProgram::BLASTX;
// }
//
// constexpr BlastFormatProgram
// getBlastProgramType(Dna5 const &, AminoAcid const &)
// {
//     return BlastFormatProgram::BLASTX;
// }
//
// // --- Protein vs Dna ---
// constexpr BlastFormatProgram
// getBlastProgramType(AminoAcid const &, Dna const &)
// {
//     return BlastFormatProgram::TBLASTX;
// }
//
// constexpr BlastFormatProgram
// getBlastProgramType(AminoAcid const &, Dna5 const &)
// {
//     return BlastFormatProgram::TBLASTX;
// }

// ============================================================================
// Functions
// ============================================================================

/*!
 * @defgroup BlastScoringScheme
 * @brief functions for converting to and from Blast's scoring-scheme behaviour
 *
 * Blast (and many other tools) compute scores of a stretch of gaps as
 * <tt>s = gO + n * gE</tt>
 * where gO is the gapOpen score, gE is the gap extend score and n ist the
 * total number of gaps.
 *
 * SeqAn, however, computes them as as <tt>s = gO + (n-1) * gE</tt>. The
 * functions below convert between the behaviours by adjusting the
 * gapOpen score.
 *
 * For more information, see <a href="https://trac.seqan.de/ticket/1091">https://trac.seqan.de/ticket/1091</a>.
 *
 * Please note that independent of this issue, SeqAn always works with
 * scores, never with penalties, i.e. a penalty is represented by a negative
 * score.
 *
 * @fn BlastScoringScheme#seqanScoringScheme2blastScoringScheme
 * @signature void seqanScoringScheme2blastScoringScheme(scoringScheme);
 * @brief convert to Blast's behaviour
 * @param[in,out]      scoringScheme      The @link Score @endlink object to modify.
 *
 * @fn BlastScoringScheme#blastScoringScheme2seqanScoringScheme
 * @signature void blastScoringScheme2seqanScoringScheme(scoringScheme);
 * @brief convert from Blast's behaviour
 * @param[in,out]      scoringScheme      The @link Score @endlink object to modify.
 *
 */

template <typename TValue, typename TSpec>
inline void
seqanScoringScheme2blastScoringScheme(Score<TValue, TSpec> & scheme)
{
    setScoreGapOpen(scheme, scoreGapOpen(scheme) - scoreGapExtend(scheme));
}

template <typename TValue, typename TSpec>
inline void
blastScoringScheme2seqanScoringScheme(Score<TValue, TSpec> & scheme)
{
    setScoreGapOpen(scheme, scoreGapOpen(scheme) + scoreGapExtend(scheme));
}

// ----------------------------------------------------------------------------
// _programTagToString()
// ----------------------------------------------------------------------------

template <BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
constexpr const char *
_programTagToString(BlastFormat<f, p, g> const & /*tag*/)
{
    return "UNKOWN BLAST PROGRAM";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr const char *
_programTagToString(BlastFormat<f, BlastFormatProgram::BLASTN, g> const &)
{
    return "BlastN";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr const char *
_programTagToString(BlastFormat<f, BlastFormatProgram::BLASTP, g> const &)
{
    return "BLASTP";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr const char *
_programTagToString(BlastFormat<f, BlastFormatProgram::BLASTX, g> const &)
{
    return "BLASTX";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr const char *
_programTagToString(BlastFormat<f, BlastFormatProgram::TBLASTN, g> const &)
{
    return "TBLASTN";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr const char *
_programTagToString(BlastFormat<f, BlastFormatProgram::TBLASTX, g> const &)
{
    return "TBLASTX";
}

// ----------------------------------------------------------------------------
// Function writeTop()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#writeTop
 * @headerfile seqan/blast.h
 * @brief write the top-most section of a BLAST output file (NO-OP for tabular formats)
 * @signature writeTop(stream, blastDbSpecs, blastFormatTag)
 *
 * @param[in,out] stream     The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in] blastDbSpecs   The @link BlastDbSpecs @endlink of your database-
 * @param[in] blastFormatTag The @link BlastFormat @endlink specifier.
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastMatch
 */

template <typename TStream,
          typename TDbSpecs,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeTop(TStream                    & /**/,
         TDbSpecs             const & /**/,
         BlastFormat<f, p, g> const & /*tag*/)
{
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastRecord#writeRecord
 * @headerfile seqan/blast.h
 * @brief write a @link BlastRecord @endlink including it's
 *  @link BlastMatch @endlinkes and possible headers to a file.
 * @signature writeRecord(stream, blastRecord, blastDbSpecs, [fieldList,] blastFormatTag)
 *
 * @param[in,out] stream     The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in] blastRecord    The @link BlastRecord @endlink you wish to print.
 * @param[in] blastDbSpecs   The @link BlastDbSpecs @endlink.
 * @param[in] fieldList      A Sequence of @link BlastMatchField::Enum @endlink to print.
 * @param[in] blastFormatTag The @link BlastFormat @endlink specifier.
 *
 * The fieldList parameter may only be specified when @link BlastFormatFile
 * @endlink is TABULAR or TABULAR_WITH_HEADER and if @link
 * BlastFormatGeneration @endlink is BLAST_PLUS.
 *
 * @see BlastFormat
 * @see BlastRecord
 */

// ----------------------------------------------------------------------------
// Function writeBottom()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#writeBottom
 * @headerfile seqan/blast.h
 * @brief write the top-most section of a BLAST output file (NO-OP for tabular formats)
 * @signature writeBottom(stream, blastDbSpecs, scoringAdapter, blastFormatTag)
 *
 * @param[in,out] stream            The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in] scoringAdapter    A @link BlastScoringAdapter @endlink with relevant information.
 * @param[in] blastDbSpecs      The @link BlastDbSpecs @endlink of your database-
 * @param[in] blastFormatTag The @link BlastFormat @endlink specifier.
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastMatch
 */

template <typename TStream,
          typename TDbSpecs,
          typename TBlastScoringAdapater,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeBottom(TStream                           & /**/,
            TDbSpecs                    const & /**/,
            TBlastScoringAdapater       const & /**/,
            BlastFormat<f, p,g>         const & /*tag*/)
{
}

// ----------------------------------------------------------------------------
// Function _untranslatePositions() -- retransform positions
// ----------------------------------------------------------------------------

template <typename TPos>
inline void
_untranslatePositions(TPos & effectiveStart,
                      TPos & /**/,
                      signed char const /**/,
                      TPos const & /**/,
                      False const & /*hasReverseComplement*/,
                      False const & /*hasFrames*/)
{
    // BLAST is 1-indexed, but end positions are "on" instead of behind
    // so only the begin positions need adapting
    ++effectiveStart;
}

template <typename TPos>
inline void
_untranslatePositions(TPos & effectiveStart,
                      TPos & effectiveEnd,
                      signed char const frameShift,
                      TPos const & length,
                      True const & /*hasReverseComplement*/,
                      False const & /*hasFrames*/)
{
    if (frameShift > 0)
    {
        // BLAST is 1-indexed, but end positions are "on" instead of behind
        // so only the begin positions need adapting
        ++effectiveStart;
    } else
    {
        // reverse strand coordinates have to be transformed
        effectiveStart = length - effectiveStart;
        effectiveEnd = length - effectiveEnd + 1;
        // end is incremented instead of start
    }
}

template <typename TPos>
inline void
_untranslatePositions(TPos & effectiveStart,
                      TPos & effectiveEnd,
                      signed char const frameShift,
                      TPos const & length,
                      True const & /*hasReverseComplement*/,
                      True const & /*hasFrames*/)
{
    // correct for codon translation and frameshift
    // subtract 1 because frameshift is 1-indexed
    effectiveStart = effectiveStart * 3 + std::abs(frameShift) - 1;
    effectiveEnd = effectiveEnd * 3 + std::abs(frameShift) - 1;

    _untranslatePositions(effectiveStart, effectiveEnd, frameShift,
                          length, True(), False());
}

// ----------------------------------------------------------------------------
// Function _nextPos() -- increment/decrement position
// ----------------------------------------------------------------------------

constexpr int8_t
_step(signed char const /**/,
      False const & /*hasReverseComplement*/,
      False const & /*hasFrames*/)
{
    return 1;
}

constexpr int8_t
_step(signed char const frameShift,
      True const & /*hasReverseComplement*/,
      False const & /*hasFrames*/)
{
    // iterate backwards for reverse frames
    return (frameShift < 0) ? -1 : 1;
}

constexpr int8_t
_step(signed char const frameShift,
      True const & /*hasReverseComplement*/,
      True const & /*hasFrames*/)
{
    // iterate three nucleotides per amino acid
    return (frameShift < 0) ? -3 : 3;
}

} // namespace seqan

#endif // header guard
