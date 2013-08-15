// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STORE_ALIGN_H
#define SEQAN_HEADER_STORE_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Aligned Read Store
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class AlignedReadStoreElement
 * @headerfile <seqan/store.h>
 * @brief Represents an alignment of a read to a contig.
 *
 * @signature template <typename TPos, typename TGapAnchors[, typename TSpec]>
 *            struct AlignedReadStoreElement;
 *
 * @tparam TPos        The position type to use.
 * @tparam TGapAnchors Type of a read @link GapAnchor @endlink.
 * @tparam TSpec       The specializing type.  Default: <tt>void</tt>.
 *
 * Value type of the @link FragmentStore::alignedReadStore @endlink.  In contrast to all other @link FragmentStore
 * @endlink stores, the @link AlignedReadStoreElement::id @endlink of an aligned read is explicitely stored as a member
 * to allow for reordering the @link FragmentStore::alignedReadStore @endlink and still recover the id.
 */

/*!
 * @fn AlignedReadStoreElement::AlignedReadStoreElement
 * @brief Constructor.
 *
 * @signature AlignedReadStoreElement::AlignedReadStoreElement();
 * @signature AlignedReadStoreElement::AlignedReadStoreElement(id, readId, contigId, beginPos, endPos[, gaps]);
 *
 * @param[in] id       The alignment id refers to associated alignment information in @link
 *                     FragmentStore::alignQualityStore @endlink or @link FragmentStore::alignedReadTagStore @endlink.
 * @param[in] readId   Refers to the aligned read in the @link FragmentStore::readStore @endlink.
 * @param[in] contigID Refers to the contig in the @link FragmentStore::contigStore @endlink the read is aligned to.
 * @param[in] beginPos Begin position of the alignment in gap-space.
 * @param[in] endPos   End position of the alignment in gap-space.
 * @param[in] gaps     A @link String @endlink of @link GapAnchor @endlink objects.
 *
 * The default constructors sets all ids to <tt>AlignedReadStoreElement::INVALID_ID</tt> and @link
 * AlignedReadStoreElement::beginPos @endlink and @link AlignedReadStoreElement::endPos @endlink to 0.
 */

 /*!
 * @typedef AlignedReadStoreElement::TId
 * @brief Type of all store ids.
 *
 * <tt>TId</tt> is the result of <tt>Id&lg;AlignedReadStoreElement&lgt;&gt; &gt;::Type</tt>, see @link Id @endlink.
 *
 * @typedef AlignedReadStoreElement::TPos
 * @brief Type of the @link AlignedReadStoreElement::beginPos @endlink and @link AlignedReadStoreElement::endPos
 *        @endlink.
 *
 * @typedef AlignedReadStoreElement::TGapAnchors
 * @brief Type of the @link AlignedReadStoreElement::gaps @endlink member.
 *
 * @typedef AlignedReadStoreElement::TSpec
 * @brief The specializing type.
 */

/*!
 * @var TId AlignedReadStoreElement::INVALID_ID;
 * @brief Constant expression that is the value of an invalid id.
 *
 * @var TId AlignedReadStoreElement::id;
 * @brief The alignment id refers to associated alignment information in @link FragmentStore::alignQualityStore
 *        @endlink or @link FragmentStore::alignedReadTagStore @endlink.
 *
 * @var TId AlignedReadStoreElement::readId;
 * @brief Refers to the aligned read in @link FragmentStore::readStore @endlink.
 *
 * @var TId AlignedReadStoreElement::contigId;
 * @brief Refers to the contig in the @link FragmentStore::contigStore @endlink that the read is aligned with.
 *
 * @var TId AlignedReadStoreElement::pairMatchId;
 * @brief Two read alignments having the same pairMatchId form a valid pair match.  If it equals <tt>INVALID_ID</tt>
 *        then the read is either not paired or could not be aligned as part of a pair.
 *
 * @var TPos AlignedReadStoreElement::beginPos;
 * @brief Begin position in gap-space.
 *
 * If beginPos &lt; endPos then the read is aligned to the reverse strand where beginPos and endPos are the
 * corresponding positions on the forward strand.
 *
 * @var TPos AlignedReadStoreElement::endPos;
 * @brief End position in gap-space.
 *
 * If beginPos &lt; endPos then the read is aligned to the reverse strand where beginPos and endPos are the
 * corresponding positions on the forward strand.
 *
 * @var TGaps AlignedReadStoreElement::gaps;
 * @brief String of read @link GapAnchor @endlink objects.  Can be used to create a @link AnchorGaps @endlink
 *        alignment row.
 */

/**
.Class.AlignedReadStoreElement
..summary:Represents an alignment between read and contig.
..cat:Fragment Store
..signature:AlignedReadStoreElement<>
..signature:AlignedReadStoreElement<TPos[, TGapAnchor[, TSpec]]>
..param.TPos:Type to store (gap-space) positions.
..param.TGapAnchor:Type of a read gap anchor.
...type:Class.GapAnchor
..param.TSpec:The specialization type.
...default:$void$
..remarks:Value type of the @Memvar.FragmentStore#alignedReadStore@ string.
In contrast to all other @Class.FragmentStore@ stores, the @Memvar.AlignedReadStoreElement#id@ of an aligned read is explicitly stored
as a member to allow for random reordering of the @Memvar.FragmentStore#alignedReadStore@, e.g. sorting by genomic position via @Function.sortAlignedReads@.
..include:seqan/store.h

.Typedef.AlignedReadStoreElement#TId
..summary:Type of all stored ids.
..remarks:$TId$ equals the result of $Id<AlignedReadStoreElement<> >::Type$, see @Metafunction.Id@.
..class:Class.AlignedReadStoreElement
.Typedef.AlignedReadStoreElement#TPos
..summary:Type of the @Memvar.AlignedReadStoreElement#beginPos@ and @Memvar.AlignedReadStoreElement#endPos@ members.
..class:Class.AlignedReadStoreElement
.Typedef.AlignedReadStoreElement#TGapAnchors
..summary:Type of the @Memvar.AlignedReadStoreElement#gaps@ member.
..class:Class.AlignedReadStoreElement
.Typedef.AlignedReadStoreElement#TSpec
..summary:The specialization type.
..class:Class.AlignedReadStoreElement

.Memfunc.AlignedReadStoreElement#AlignedReadStoreElement
..summary:Constructor
..signature:AlignedReadStoreElement()
..signature:AlignedReadStoreElement(id, readId, contigId, beginPos, endPos[, gaps])
..param.id:The alignment id refers to associated alignment information in @Memvar.FragmentStore#alignQualityStore@ or @Memvar.FragmentStore#alignedReadTagStore@.
..param.readId:Refers to the aligned read in the @Memvar.FragmentStore#readStore@.
..param.contigId:Refers to the contig in the @Memvar.FragmentStore#contigStore@ the read is aligned with.
..param.beginPos:Begin position of the alignment in gap-space.
..param.endPos:End position of the alignment in gap-space.
..param.gaps:Read gap anchors.
..remarks:The default constructor sets all ids to @Memvar.AlignedReadStoreElement#INVALID_ID@ and
@Memvar.AlignedReadStoreElement#beginPos@ and @Memvar.AlignedReadStoreElement#endPos@ to $0$.
..class:Class.AlignedReadStoreElement

.Memvar.AlignedReadStoreElement#id
..summary:The alignment id refers to associated alignment information in @Memvar.FragmentStore#alignQualityStore@ or @Memvar.FragmentStore#alignedReadTagStore@.
..type:Metafunction.Id
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#readId
..summary:Refers to the aligned read in the @Memvar.FragmentStore#readStore@.
..type:Metafunction.Id
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#contigId
..summary:Refers to the contig in the @Memvar.FragmentStore#contigStore@ the read is aligned with.
..type:Metafunction.Id
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#pairMatchId
..summary:Two read alignments having the same @Memvar.AlignedReadStoreElement#pairMatchId@ form a valid pair match.
If it equals @Memvar.AlignedReadStoreElement#INVALID_ID@, the read is either not paired or could not be aligned as part of a pair match.
..type:Metafunction.Id
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#beginPos
..summary:Begin position of the alignment in gap-space.
..type:Typedef.AlignedReadStoreElement#TPos
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#endPos
..summary:End position of the alignment in gap-space. If @Memvar.AlignedReadStoreElement#endPos@ < @Memvar.AlignedReadStoreElement#beginPos@, the read is aligned to the reverse strand, 
where @Memvar.AlignedReadStoreElement#beginPos@ and @Memvar.AlignedReadStoreElement#endPos@ are the corresponding positions on the forward strand.
..type:Typedef.AlignedReadStoreElement#TPos
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#gaps
..summary:String of read gap anchors. Can be used to create a @Spec.AnchorGaps@ alignment row.
..type:Typedef.AlignedReadStoreElement#TGapAnchors
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#INVALID_ID
..summary:Constant to represent an invalid id.
..type:Metafunction.Id
..class:Class.AlignedReadStoreElement
*/

template <typename TPos_, typename TGapAnchor_, typename TSpec_ = void>
struct AlignedReadStoreElement
{
	typedef typename Id<AlignedReadStoreElement>::Type	TId;
	typedef TPos_										TPos;
	typedef TGapAnchor_									TGapAnchor;
	typedef TSpec_										TSpec;
	typedef String<TGapAnchor>							TGapAnchors;

	static const TId INVALID_ID;
	
	TId			id;
	TId			readId;
	TId			contigId;
	TId			pairMatchId;	// unique id. for multiple mate-pair matches (not matePairId)
	TPos		beginPos;		// begin position of the gapped sequence in gapped contig sequence
	TPos		endPos;			// end position of ..., for reverse aligned reads holds end < begin
	TGapAnchors	gaps;

	AlignedReadStoreElement() : id(INVALID_ID), readId(INVALID_ID), contigId(INVALID_ID), pairMatchId(INVALID_ID), beginPos(0), endPos(0) {}

	AlignedReadStoreElement(TId _id, TId _readId, TId _contigId, TPos _beginPos, TPos _endPos) :
		id(_id), 
		readId(_readId), 
		contigId(_contigId), 
		pairMatchId(INVALID_ID), 
		beginPos(_beginPos), 
		endPos(_endPos) {}

	AlignedReadStoreElement(TId _id, TId _readId, TId _contigId, TPos _beginPos, TPos _endPos, TGapAnchors const &_gaps) :
		id(_id), 
		readId(_readId), 
		contigId(_contigId), 
		pairMatchId(INVALID_ID), 
		beginPos(_beginPos), 
		endPos(_endPos),
		gaps(_gaps) {}

    inline bool operator==(AlignedReadStoreElement const & other) const
    {
        return id == other.id &&
                readId == other.readId &&
                contigId == other.contigId &&
                pairMatchId == other.pairMatchId &&
                beginPos == other.beginPos &&
                endPos == other.endPos &&
                gaps == other.gaps;
    }
};


// TODO(holtgrew): I find this useful for debugging purposes. Keep it?
template <typename TStream, typename TPos, typename TGapAnchor, typename TSpec>
TStream & operator<<(TStream & stream, AlignedReadStoreElement<TPos, TGapAnchor, TSpec> const & e) {
    return stream << "AlignedReadStore(id=" << e.id << ", readId=" << e.readId << ", contigId=" << e.contigId << ", pairMatchId=" << e.pairMatchId << ", beginPos=" << e.beginPos << ", endPos=" << e.endPos << ", {gaps})";
}


//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec> 
const typename Id<AlignedReadStoreElement<TPos, TGapAnchor, TSpec> >::Type 
AlignedReadStoreElement<TPos, TGapAnchor, TSpec>::INVALID_ID = MaxValue<typename Id<AlignedReadStoreElement<TPos, TGapAnchor, TSpec> >::Type>::VALUE; 

//////////////////////////////////////////////////////////////////////////////

/*!
 * @class AlignQualityStoreElement
 * @headerfile <seqan/store.h>
 * @brief Stores the alignment qualities.
 *
 * @signature template <typename TScore[, typename TSpec]>
 *            struct AlignQualityStoreElement;
 *
 * @tparam TScore Type to store align and pair score values.
 * @tparam TSpec  The specializing type.  Default: <tt>void</tt>.
 *
 * Value type of @link FragmentStore::alignQualityStore @endlink string.
 *
 *
 * @fn AlignQualityStoreElement::AlignQualityStoreElement
 * @brief Constructor, sets all members to 0.
 *
 * @signature AlignQualityStoreElement::AlignQualityStoreElement();
 * 
 *
 * @var TScore AlignQualityStoreElement::pairScore;
 * @brief Combined score of both alignmetns and pair match.
 *
 * @var TScore AlignQualityStoreElement::score;
 * @brief Score of the alignment.
 *
 * @var TCount AlignQualityStoreElement::errors;
 * @brief Absolute number of errors in the alignment (<tt>unsigned char</tt>).
 */

/**
.Class.AlignQualityStoreElement
..summary:Stores alignment qualities.
..cat:Fragment Store
..signature:AlignQualityStoreElement<TScore[, TSpec]>
..param.TScore:Type to store align and pair score values.
..param.TSpec:The specialization type.
...default:$void$
..remarks:Value type of the @Memvar.FragmentStore#alignQualityStore@ string.

.Memfunc.AlignQualityStoreElement#AlignQualityStoreElement
..summary:Constructor
..signature:AlignQualityStoreElement<TScore[, TSpec]> ()
..remarks:Sets all members to $0$.

..class:Class.AlignQualityStoreElement
.Memvar.AlignQualityStoreElement#pairScore
..summary:Combined score of both alignments of a pair match.
..class:Class.AlignQualityStoreElement
.Memvar.AlignQualityStoreElement#score
..summary:Score of the alignment.
..class:Class.AlignQualityStoreElement
.Memvar.AlignQualityStoreElement#errors
..summary:Absolute number of errors in the alignment.
..type:nolink:unsigned char
..class:Class.AlignQualityStoreElement
..include:seqan/store.h
*/

template <typename TScore, typename TSpec = void>
struct AlignQualityStoreElement
{
	TScore				pairScore;		// score of the mate-pair alignment (this read is part of)
	TScore				score;			// score of the single read alignment
	unsigned char		errors;			// absolute number of errors (Hamming or edit distance)
	
	AlignQualityStoreElement():
		pairScore(0),
		score(0),
		errors(0) {}

	AlignQualityStoreElement(TScore _pairScore, TScore _score, unsigned char _errors):
		pairScore(_pairScore),
		score(_score),
		errors(_errors) {}

    inline bool operator==(AlignQualityStoreElement const & other)
    {
        return pairScore == other.pairScore &&
                score == other.score &&
                errors == other.errors;
    }
};


//////////////////////////////////////////////////////////////////////////////
// Sorting tags
//////////////////////////////////////////////////////////////////////////////

/*!
 * @defgroup SortAlignedReadTags Tags for sortAlignedReads
 * @brief Tags to select a specific field to sort the @link FragmentStore::alignedReadStore @endlink by.
 *
 * @see sortAlignedReads
 * @see lowerBoundAlignedReads
 * @see upperBoundAlignedReads
 *
 * @tag SortAlignedReadTags#SortContigId
 * @headerfile <seqan/store.h>
 * @brief Sort aligned reads by @link AlignedReadStoreElement::contigId @endlink.
 *
 * @signature typedef Tag<SortContigId_> const SortContigId;
 *
 * @tag SortAlignedReadTags#SortId
 * @headerfile <seqan/store.h>
 * @brief Sort aligned reads by @link AlignedReadStoreElement::id @endlink.
 *
 * @signature typedef Tag<SortId_> const SortId;
 *
 * @tag SortAlignedReadTags#SortBeginPos
 * @headerfile <seqan/store.h>
 * @brief Sort aligned reads by @link AlignedReadStoreElement::beginPos @endlink.
 *
 * @signature typedef Tag<SortBeginPos_> const SortBeginPos;
 *
 * @tag SortAlignedReadTags#SortEndPos
 * @headerfile <seqan/store.h>
 * @brief Sort aligned reads by @link AlignedReadStoreElement::endPos @endlink.
 *
 * @signature typedef Tag<SortEndPos_> const SortEndPos;
 *
 * @tag SortAlignedReadTags#SortPairMatchId
 * @headerfile <seqan/store.h>
 * @brief Sort aligned reads by @link AlignedReadStoreElement::pairMatchId @endlink.
 *
 * @signature typedef Tag<SortPairMatchId_> const SortPairMatchId;
 *
 * @tag SortAlignedReadTags#SortReadId
 * @headerfile <seqan/store.h>
 * @brief Sort aligned reads by @link AlignedReadStoreElement::readId @endlink.
 *
 * @signature typedef Tag<SortReadId_> const SortReadId;
 */

/**
.Tag.sortAlignedRead Tags
..summary:Tag to select a specific field to stably sort the @Memvar.FragmentStore#alignedReadStore@ by.
..cat:Fragment Store
..see:Function.sortAlignedReads
..see:Function.lowerBoundAlignedReads
..see:Function.upperBoundAlignedReads
..tag.SortContigId:
...summary:Sort alignedReads by @Memvar.AlignedReadStoreElement#contigId@.
...signature:SortContigId
..tag.SortId:
...summary:Sort alignedReads by @Memvar.AlignedReadStoreElement#id@.
...signature:SortId
..tag.SortBeginPos:
...summary:Sort alignedReads by @Memvar.AlignedReadStoreElement#beginPos@.
...signature:SortBeginPos
..tag.SortEndPos:
...summary:Sort alignedReads by @Memvar.AlignedReadStoreElement#endPos@.
...signature:SortEndPos
..tag.SortPairMatchId:
...summary:Sort alignedReads by @Memvar.AlignedReadStoreElement#pairMatchId@.
...signature:SortPairMatchId
..tag.SortReadId:
...summary:Sort alignedReads by @Memvar.AlignedReadStoreElement#readId@.
...signature:SortReadId
..include:seqan/store.h
*/

struct SortContigId_;
typedef Tag<SortContigId_> const SortContigId;

struct SortId_;
typedef Tag<SortId_> const SortId;

struct SortBeginPos_;
typedef Tag<SortBeginPos_> const SortBeginPos;

struct SortEndPos_;
typedef Tag<SortEndPos_> const SortEndPos;

struct SortPairMatchId_;
typedef Tag<SortPairMatchId_> const SortPairMatchId;

struct SortReadId_;
typedef Tag<SortReadId_> const SortReadId;


//////////////////////////////////////////////////////////////////////////////
// Sorting functors
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec>
inline void swap(
    AlignedReadStoreElement<TPos, TGapAnchor, TSpec> &a,
    AlignedReadStoreElement<TPos, TGapAnchor, TSpec> &b)
{
    std::swap(a.id, b.id);
    std::swap(a.readId, b.readId);
    std::swap(a.contigId, b.contigId);
    std::swap(a.pairMatchId, b.pairMatchId);
    std::swap(a.beginPos, b.beginPos);
    std::swap(a.endPos, b.endPos);
    swap(a.gaps, b.gaps);
}

template <typename TAlignedRead, typename TTag>
struct _LessAlignedRead;

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return (a1.id) < (a2.id);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortContigId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.contigId < a2.contigId;
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortBeginPos> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return _min(a1.beginPos, a1.endPos) < _min(a2.beginPos, a2.endPos);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortEndPos> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return _max(a1.beginPos, a1.endPos) < _max(a2.beginPos, a2.endPos);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortPairMatchId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.pairMatchId < a2.pairMatchId;
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortReadId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.readId < a2.readId;
	}
};

//////////////////////////////////////////////////////////////////////////////
// Sorting function
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): Maybe add equalRangeAlignedReads?

/*!
 * @fn sortAlignedReads
 * @headerfile <seqan/store.h>
 * @brief Stably read alignments, e.g. in @link FragmentStore::alignedReadStore @endlink.
 *
 * @signature void sortAlignedReads(alignStore, sortTag);
 * @signature void sortAlignedReads(alignStore, lessFunctor);
 *
 * @param[in,out] alignStore  The @link SequenceConcept sequence @endlink of @link AlignedReadStoreElement @endlink
 *                            to be sorted, e.g. @link FragmentStore::alignedReadStore @endlink.
 * @param[in]     sortTag     Tag for selecting the member to sort by.  See @link SortAlignedReadTags @endlink.
 * @param[in]     lessFunctor A functor to pass to <tt>std::stable_sort</tt> for sorting the sequence.
 *
 * @see SortAlignedReadTags
 * @see lowerBoundAlignedReads
 * @see upperBoundAlignedReads
 */

/*!
 * @fn lowerBoundAlignedReads
 * @headerfile <seqan/store.h>
 * @brief Performs a binary lower bound search on read alignments.
 *
 * @signature TIter1 lowerBoundAlignedReads(alignStore, value, sortTag);
 * @signature TIter2 lowerBoundAlignedReads(itBegin, itEnd, value, sortTag);
 *
 * @param[in,out] alignStore  The @link SequenceConcept sequence @endlink of @link AlignedReadStoreElement @endlink
 *                            to be searched, e.g. @link FragmentStore::alignedReadStore @endlink.
 * @param[in]     itBegin     Iterator to the begin of the sequence to search.
 * @param[in]     itEnd       Iterator to the end of the sequence to search.
 * @param[in]     value       The value to search for.
 * @param[in]     sortTag     Tag for selecting the member to compare by.  See @link SortAlignedReadTags @endlink.
 *
 * @return TIter1 Iterator to the lower bound item.  If <tt>TAlignStore</tt> is the type of <tt>alignStore</tt> then
 *                TIter1 is the result of <tt>Iterator&lt;TAlignStore, Standard&gt;::Type</tt>.
 * @return TIter2 Iterator to the lower bound item.  Has the same type as <tt>itBegin</tt> and <tt>itEnd</tt>.
 *
 * @see SortAlignedReadTags
 * @see sortAlignedReads
 * @see upperBoundAlignedReads
 */

/*!
 * @fn upperBoundAlignedReads
 * @headerfile <seqan/store.h>
 * @brief Performs a binary upper bound search on read alignments.
 *
 * @signature TIter1 upperBoundAlignedReads(alignStore, value, sortTag);
 * @signature TIter2 upperBoundAlignedReads(itBegin, itEnd, value, sortTag);
 *
 * @param[in,out] alignStore  The @link SequenceConcept sequence @endlink of @link AlignedReadStoreElement @endlink
 *                            to be searched, e.g. @link FragmentStore::alignedReadStore @endlink.
 * @param[in]     itBegin     Iterator to the begin of the sequence to search.
 * @param[in]     itEnd       Iterator to the end of the sequence to search.
 * @param[in]     value       The value to search for.
 * @param[in]     sortTag     Tag for selecting the member to compare by.  See @link SortAlignedReadTags @endlink.
 *
 * @return TIter1 Iterator to the upper bound item.  If <tt>TAlignStore</tt> is the type of <tt>alignStore</tt> then
 *                TIter1 is the result of <tt>Iterator&lt;TAlignStore, Standard&gt;::Type</tt>.
 * @return TIter2 Iterator to the upper bound item.  Has the same type as <tt>itBegin</tt> and <tt>itEnd</tt>.
 *
 * @see SortAlignedReadTags
 * @see sortAlignedReads
 * @see lowerBoundAlignedReads
 */

/**
.Function.sortAlignedReads
..summary:Stably sort aligned reads.
..cat:Fragment Store
..signature:sortAlignedReads(alignStore, sortTag)
..signature:sortAlignedReads(alignStore, lessFunctor)
..param.alignStore:A sequence of @Class.AlignedReadStoreElement@ to be sorted, e.g. @Memvar.FragmentStore#alignedReadStore@.
..param.sortTag:Selects the field to sort by.
...type:Tag.sortAlignedRead Tags
..param.lessFunctor:STL-less functor to compare two @Class.AlignedReadStoreElement.AlignedReadStoreElements@.
..remarks:This function calls $std::stable_sort$ to sort the @Memvar.FragmentStore#alignedReadStore@.
..include:seqan/store.h
..see:Function.lowerBoundAlignedReads
..see:Function.upperBoundAlignedReads

.Function.lowerBoundAlignedReads
..summary:Performs a binary lower bound search on the aligned reads.
..cat:Fragment Store
..signature:lowerBoundAlignedReads(alignStore, value, sortTag)
..signature:lowerBoundAlignedReads(itBegin, itEnd, value, sortTag)
..param.alignStore:A sequence of @Class.AlignedReadStoreElement@ to be searched through, e.g. @Memvar.FragmentStore#alignedReadStore@.
..param.itBegin:An iterator to the first element of the sequence of @Class.AlignedReadStoreElement@ to be searched through.
..param.itEnd:An iterator behind the last element of the sequence of @Class.AlignedReadStoreElement@ to be searched through.
..param.value:The value to use for the comparison.
..param.sortTag:Selects the field for the comparison in the binary search.
...type:Tag.sortAlignedRead Tags
..remarks:This is equivalent to calling $std::lower_bound$ on @Memvar.FragmentStore#alignedReadStore@ with according parameters.
..include:seqan/store.h
..see:Function.sortAlignedReads
..see:Function.upperBoundAlignedReads

.Function.upperBoundAlignedReads
..summary:Performs a binary upper bound search on the aligned reads.
..cat:Fragment Store
..signature:upperBoundAlignedReads(alignStore, value, sortTag)
..signature:upperBoundAlignedReads(itBegin, itEnd, value, sortTag)
..param.alignStore:A sequence of @Class.AlignedReadStoreElement@ to be searched through, e.g. @Memvar.FragmentStore#alignedReadStore@.
..param.itBegin:An iterator to the first element of the sequence of @Class.AlignedReadStoreElement@ to be searched through.
..param.itEnd:An iterator behind the last element of the sequence of @Class.AlignedReadStoreElement@ to be searched through.
..param.value:The value to use for the comparison.
..param.sortTag:Selects the field for the comparison in the binary search.
...type:Tag.sortAlignedRead Tags
..remarks:This is equivalent to calling $std::upper_bound$ on @Memvar.FragmentStore#alignedReadStore@ with according parameters.
..include:seqan/store.h
..see:Function.sortAlignedReads
..see:Function.lowerBoundAlignedReads
*/

template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign& alignStore, Tag<TSortSpec> const &) 
{
	std::stable_sort(
		begin(alignStore, Standard() ), 
		end(alignStore, Standard() ), 
		_LessAlignedRead<typename Value<TAlign>::Type, Tag<TSortSpec> const>() );
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign const & alignStore, Tag<TSortSpec> const &) 
{
	std::stable_sort(
		begin(const_cast<TAlign&>(alignStore), Standard() ), 
		end(const_cast<TAlign&>(alignStore), Standard() ), 
		_LessAlignedRead<typename Value<TAlign>::Type, Tag<TSortSpec> const>() );
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign & alignStore, TFunctorLess const &less) 
{
	std::stable_sort(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		less);
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign const & alignStore, TFunctorLess const &less) 
{
	std::stable_sort(
		begin(const_cast<TAlign&>(alignStore), Standard()), 
		end(const_cast<TAlign&>(alignStore), Standard()), 
		less);
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
lowerBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortId const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
upperBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortId const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
lowerBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigId const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
upperBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigId const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
lowerBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
upperBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
lowerBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
upperBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
lowerBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchId const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
upperBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchId const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
lowerBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadId const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign const, Standard>::Type
upperBoundAlignedReads(TAlign const & alignStore, 
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadId const>() );
}

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign & alignStore, 
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortId const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortBeginPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortBeginPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortEndPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortEndPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortBeginPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortBeginPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortEndPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortEndPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
