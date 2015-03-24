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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Functions for the normalization of gaps in alignments.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_NORMALIZE_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_NORMALIZE_H_

#include <seqan/align.h>

namespace seqan {

// Helper type for normalizing pairwise alignments.

template <typename TGapsH, typename TGapsV>
class PairwiseAlignmentNormalizer_
{
    typedef typename Iterator<TGapsH, Standard>::Type TIteratorH;
    typedef typename Iterator<TGapsV, Standard>::Type TIteratorV;

public:
    PairwiseAlignmentNormalizer_(TGapsH & gapsH, TGapsV & gapsV) :
            gapsH(gapsH), gapsV(gapsV)
    {
        SEQAN_ASSERT_EQ(length(gapsH), length(gapsV));
    }

    bool run()
    {
        bool b1 = deleteAllGapsColumns();
        bool b2 = removeNullifyingIndels();
        bool b3 = leftAlignIndels();
        return (b1 || b2 || b3);
    }

    // Remove columns that consist exclusively of gaps.
    bool deleteAllGapsColumns()
    {
        bool anyChange = false;

        // std::cerr << "GAPS H\t" << gapsH << "\n";
        // std::cerr << "GAPS V\t" << gapsV << "\n";

        unsigned posH = 0, posV = 0;
        while (posH != length(gapsH))
        {
            if (isGap(gapsH, posH) && isGap(gapsV, posV))
            {
                anyChange = true;
                unsigned n = std::min(countGaps(gapsH, posH), countGaps(gapsV, posV));
                removeGaps(iter(gapsH, posH, Standard()), n);
                removeGaps(iter(gapsV, posV, Standard()), n);
            }
            else
            {
                ++posH;
                ++posV;
            }
        }

        return anyChange;
    }

    // Remove consecutive columns having having an insertion XOR a deletion.
    //
    // Precondition: No all-gaps columns.
    bool removeNullifyingIndels()
    {
        bool anyChange = false;
        
        // build list of operations (pseudo-CIGAR string)
        std::vector<std::pair<char, unsigned> > ops;  // (op, num), op in ['I', 'D', 'M']
        {
            TIteratorH itH = begin(gapsH, Standard());
            TIteratorV itV = begin(gapsV, Standard());

            for (; itH != end(gapsH, Standard()); ++itH, ++itV)
                if (!isGap(itH) && !isGap(itV))
                {
                    if (ops.empty() || ops.back().first != 'M')
                        ops.push_back(std::make_pair('M', 1));
                    else
                        ops.back().second += 1;
                }
                else if (isGap(itH))
                {
                    if (ops.empty() || ops.back().first != 'I')
                        ops.push_back(std::make_pair('I', 1));
                    else
                        ops.back().second += 1;
                }
                else if (isGap(itV))
                {
                    if (ops.empty() || ops.back().first != 'D')
                        ops.push_back(std::make_pair('D', 1));
                    else
                        ops.back().second += 1;
                }
        }

        // for (unsigned i = 0; i < ops.size(); ++i)
        //     std::cerr << ops[i].first << "\t" << ops[i].second << "\n";

        // Update ops to reflect the alignment we aim at.
        for (unsigned i = 0, j = 1; j < ops.size(); ++i, ++j)
        {
            if ((ops[i].first == 'I' && ops[j].first == 'D') ||
                (ops[i].first == 'D' && ops[j].first == 'I'))
            {
                if (ops[i].second == ops[j].second)
                {
                    anyChange = true;
                    ops.erase(ops.begin() + j);
                    ops[i].first = 'M';
                }
                else if (ops[i].second > ops[j].second)
                {
                    anyChange = true;
                    ops[i].second -= ops[j].second;
                    ops[j].first = 'M';
                }
                else
                {
                    anyChange = true;
                    ops[j].second -= ops[i].second;
                    ops[i].first = 'M';
                }
            }
        }

        // Apply ops to the alignment.
        clearGaps(gapsH);
        clearGaps(gapsV);
        TIteratorH itH = begin(gapsH, Standard());
        TIteratorV itV = begin(gapsV, Standard());

        for (unsigned i = 0; i < ops.size(); ++i)
        {
            switch (ops[i].first)
            {
                case 'I':
                    insertGaps(itH, ops[i].second);
                    break;
                case 'D':
                    insertGaps(itV, ops[i].second);
                    break;
                case 'M':
                    break;
                default:
                    SEQAN_FAIL("Invalid operation %c!", ops[i].first);
            }
            goFurther(itH, ops[i].second);
            goFurther(itV, ops[i].second);
        }

        return anyChange;
    }

    // Shift indels to the left.
    //
    // Precondition: Normalized gaps.
    //
    // This is following the algorithm from
    //
    //   http://genome.sph.umich.edu/wiki/Variant_Normalization
    bool leftAlignIndels()
    {
        bool anyChange = false;

        unsigned minPos = 0;  // end pos of last gap
        unsigned pos = 0;     // begin pos of next gap
        unsigned posEnd = 0;  // end pos of next gap
        while (minPos != length(gapsH))
        {
            // std::cerr << "H\t" << gapsH << "\n"
            //           << "V\t" << gapsV << "\n"
            //           << "length(gapsH) == " << length(gapsH) << "\n"
            //           << "minPos == " << minPos << "\n"
            //           << "\n";
            
            pos = minPos;

            // forward to the next gap
            while (pos < length(gapsH) && !isGap(gapsH, pos) && !isGap(gapsV, pos))
                ++pos;
            if (pos == length(gapsH))
                break;
            posEnd = pos;
            if (isGap(gapsH, pos))
            {
                while (posEnd < length(gapsH) && isGap(gapsH, posEnd))
                    ++posEnd;
            }
            else
            {
                while (posEnd < length(gapsV) && isGap(gapsV, posEnd))
                    ++posEnd;
            }
            // shift gap to the left
            unsigned sPos = pos, sPosEnd = posEnd;
            if (isGap(gapsH, pos))
            {
                while (sPos > minPos && gapsV[sPos - 1] == gapsV[sPosEnd - 1])
                {
                    --sPos;
                    --sPosEnd;
                }
                if (sPos != pos)  // any change
                {
                    anyChange = true;
                    removeGaps(gapsH, pos, posEnd - pos);
                    insertGaps(gapsH, sPos, posEnd - pos);
                }
            }
            else
            {
                while (sPos > minPos && gapsH[sPos - 1] == gapsH[sPosEnd - 1])
                {
                    --sPos;
                    --sPosEnd;
                }
                if (sPos != pos)  // any change
                {
                    anyChange = true;
                    removeGaps(gapsV, pos, posEnd - pos);
                    insertGaps(gapsV, sPos, posEnd - pos);
                }
            }
            minPos = sPos + (posEnd - pos);
        }

        // Iterate until there is no change any more, but return indicator flag from first iteration.
        if (anyChange)
            leftAlignIndels();
        return anyChange;
    }

private:
    TGapsH & gapsH;
    TGapsV & gapsV;
};

// ---------------------------------------------------------------------------
// Function removeAllGapsColumns()
// ---------------------------------------------------------------------------

/*!
 * @fn removeAllGapsColumns
 * @headerfile <seqan/align/normalize.h>
 * @brief Remove columns that only have gaps from pairwise alignments.
 *
 * @signature bool removeAllGapsColumns(gapsH, gapsV);
 * @signature bool removeAllGapsColumns(align);
 *
 * @param[in,out] gapsH @link Gaps @endlink object for the first sequence.
 * @param[in,out] gapsV @link Gaps @endlink object for the second sequence.
 * @param[in,out] align @link Align @endlink object with the pairwise alignment.
 * @return        bool  <tt>true</tt> if any columns were removed.
 */

template <typename TSeqH, typename TSpecH, typename TSeqV, typename TSpecV>
bool removeAllGapsColumns(Gaps<TSeqH, TSpecH> & gapsH, Gaps<TSeqV, TSpecV> & gapsV)
{
    PairwiseAlignmentNormalizer_<Gaps<TSeqH, TSpecH>, Gaps<TSeqV, TSpecV> > normalizer(gapsH, gapsV);
    return normalizer.deleteAllGapsColumns();
}

template <typename TSeq, typename TSpec>
bool removeAllGapsColumns(Align<TSeq, TSpec> & align)
{
    return removeAllGapsColumns(row(align, 0), row(align, 1));
}

// ---------------------------------------------------------------------------
// Function removeNullifyingIndels()
// ---------------------------------------------------------------------------

/*!
 * @fn removeNullifyingIndels
 * @headerfile <seqan/align/normalize.h>
 * @brief Compactify alignment for adjacent insertion/deletion columns.
 *
 * @signature bool removeNullifyingIndels(gapsH, gapsV);
 * @signature bool removeNullifyingIndels(align);
 *
 * @param[in,out] gapsH @link Gaps @endlink object for the first sequence.
 * @param[in,out] gapsV @link Gaps @endlink object for the second sequence.
 * @param[in,out] align @link Align @endlink object with the pairwise alignment.
 * @return        bool  <tt>true</tt> if any shifting took place.
 *
 * A precondition for this function is that there are no all-gaps columns (see @link removeAllGapsColumns @endlink.
 *
 * The alignment is considered from the left to the right and in each case, the flanking insertion/deletion (or
 * deletion/insertion) stretch is compressed as far as possible.
 */

template <typename TSeqH, typename TSpecH, typename TSeqV, typename TSpecV>
bool removeNullifyingIndels(Gaps<TSeqH, TSpecH> & gapsH, Gaps<TSeqV, TSpecV> & gapsV)
{
    PairwiseAlignmentNormalizer_<Gaps<TSeqH, TSpecH>, Gaps<TSeqV, TSpecV> > normalizer(gapsH, gapsV);
    return normalizer.removeNullifyingIndels();
}

template <typename TSeq, typename TSpec>
bool removeNullifyingIndels(Align<TSeq, TSpec> & align)
{
    return removeNullifyingIndels(row(align, 0), row(align, 1));
}

// ---------------------------------------------------------------------------
// Function leftAlignIndels()
// ---------------------------------------------------------------------------

/*!
 * @fn leftAlignIndels
 * @headerfile <seqan/align/normalize.h>
 * @brief Left-shift indels in pairwise alignments.
 *
 * @signature bool leftAlignIndels(gapsH, gapsV);
 * @signature bool leftAlignIndels(align);
 *
 * @param[in,out] gapsH @link Gaps @endlink object for the first sequence.
 * @param[in,out] gapsV @link Gaps @endlink object for the second sequence.
 * @param[in,out] align @link Align @endlink object with the pairwise alignment.
 * @return        bool  <tt>true</tt> if any left-shifting took place.
 *
 * A precondition for this function is that there are no all-gaps columns (see @link removeAllGapsColumns @endlink.
 *
 * Note that this can lead to new leading and trailing gaps, thus analyze the resulting gapsH and gapsV.
 */

template <typename TSeqH, typename TSpecH, typename TSeqV, typename TSpecV>
bool leftAlignIndels(Gaps<TSeqH, TSpecH> & gapsH, Gaps<TSeqV, TSpecV> & gapsV)
{
    PairwiseAlignmentNormalizer_<Gaps<TSeqH, TSpecH>, Gaps<TSeqV, TSpecV> > normalizer(gapsH, gapsV);
    return normalizer.leftAlignIndels();
}

template <typename TSeq, typename TSpec>
bool leftAlignIndels(Align<TSeq, TSpec> & align)
{
    return leftAlignIndels(row(align, 0), row(align, 1));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_NORMALIZE_H_
