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

#ifndef SEQAN_HEADER_GRAPH_IMPL_FRAGMENT_H
#define SEQAN_HEADER_GRAPH_IMPL_FRAGMENT_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Fragment Specs
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class ExactFragment
 * @extends Fragment
 * @headerfile <seqan/graph_types.h>
 * @brief A type for ungapped, pairwise segment matches.
 *
 * @signature template <[typename TSize[, typename TSpec]]>
 *            class Fragment<TSize, ExactFragment<TSpec> >;
 *
 * @tparam TSize The size type of the underlying sequence.  Default: <tt>Size&lt;CharString&gt;Type</tt>.
 * @tparam TSpec Specializing type.  Default: <tt>ExactFragment&lt;&gt;</tt>.
 */

/**
.Spec.ExactFragment
..cat:Alignments
..general:Class.Fragment
..summary:A type for ungapped, pairwise segment matches.
..signature:Fragment<TSize, ExactFragment<TSpec> > 
..param.TSize: The Size type of the underlying sequences.
...metafunction:Metafunction.Size
..param.TSpec:The specializing type for the graph.
...metafunction:Metafunction.Spec
..include:seqan/graph_types.h
..see:Spec.ExactReversableFragment
*/

template<typename TSpec = Default>
struct ExactFragment;	


/*!
 * @class ExactReversableFragment
 * @extends Fragment
 * @headerfile <seqan/graph_types.h>
 * @brief A type for ungapped, pairwise segment matches that maybe in reverse orientation.
 *
 * Compared to the @link ExactFragment @endlink specialzing type of @link Fragment @endlink, a @link
 * ExactReversableFragment @endlink stores an additional bool value to indicate whether a match is in reverse
 * orientation or not.
 *
 * @signature template <[typename TSize[, typename TSpec]]>
 *            class Fragment<TSize, ExactReversableFragment<TSpec> >;
 *
 * @tparam TSize The size type of the underlying sequence.  Default: <tt>Size&lt;CharString&gt;Type</tt>.
 * @tparam TSpec Specializing type.  Default: <tt>ExactFragment&lt;&gt;</tt>.
 */

/**
.Spec.ExactReversableFragment
..cat:Alignments
..general:Class.Fragment
..summary:A type for ungapped, pairwise segment matches that may be in reverse orientation.
..signature:Fragment<TSize, ExactReversableFragment<TSpec> > 
..param.TSize: The Size type of the underlying sequences.
...metafunction:Metafunction.Size
..param.TSpec:The specializing type for the graph.
...metafunction:Metafunction.Spec
..remarks:Compared to the @Spec.ExactFragment@ specialzing type of @Class.Fragment@, a @Spec.ExactReversableFragment@ stores an additional bool value to indicate whether a match is in reverse orientation or not.
..include:seqan/graph_types.h
..see:Spec.ExactFragment
*/

template<typename TSpec = Default>
struct ExactReversableFragment;	


//////////////////////////////////////////////////////////////////////////////
// Default Fragment is the exact one
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class Fragment
 * @headerfile <seqan/graph_types.h>
 * @brief A type for pairwise segment matches.
 *
 * @signature template <[typename TSize[, typename TSpec]]>
 *            class Fragment;
 *
 * @tparam TSize The size type of the underlying sequence.  Default: <tt>Size&lt;CharString&gt;Type</tt>.
 * @tparam TSpec Specializing type.  Default: <tt>ExactFragment&lt;&gt;</tt>.
 *
 * @section Examples
 *
 * @code{.cpp}
 * // Construct fragment.
 * unsigned seqId1 = 0, beg1 = 0, seqId2 = 32, beg2 = 42, len = 33;
 * Fragment<> fragment(seqId1, beg1, seqId2, beg2, len);
 * 
 * // Update fragment's properties.
 * fragmentBegin(fragment, 0) = 10;
 * fragmentBegin(fragment, 1) = 10;
 * sequenceId(fragment, 0) = 33;
 * sequenceId(fragment, 1) = 44;
 * fragmentLength(fragment) += 42;
 * @endcode
 */

/**
.Class.Fragment:
..cat:Alignments
..summary:A type for ungapped, pairwise segment matches.
..signature:Fragment<TSize, TSpec>
..param.TSize:The size type of the underlying sequences.
...metafunction:Metafunction.Size
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:@Spec.ExactFragment@
..include:seqan/graph_types.h
..example:A small example using fragments.
..example.code:
// Construct fragment.
unsigned seqId1 = 0, beg1 = 0, seqId2 = 32, beg2 = 42, len = 33;
Fragment<> fragment(seqId1, beg1, seqId2, beg2, len);

// Update fragment's properties.
fragmentBegin(fragment, 0) = 10;
fragmentBegin(fragment, 1) = 10;
sequenceId(fragment, 0) = 33;
sequenceId(fragment, 1) = 44;
fragmentLength(fragment) += 42;
*/


template<typename TSize = typename Size<String<char> >::Type, typename TSpec = ExactFragment<> >
class Fragment;


//////////////////////////////////////////////////////////////////////////////
// Size Metafunction
//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec>
struct Size<Fragment<TSize, TSpec> > {
	typedef TSize Type;
};


template<typename TSize, typename TSpec>
struct Size<Fragment<TSize, TSpec> const> {
	typedef TSize Type;
};


//////////////////////////////////////////////////////////////////////////////
// Exact Fragment
//////////////////////////////////////////////////////////////////////////////
	

template<typename TSize, typename TSpec>
class Fragment<TSize, ExactFragment<TSpec> > {
public:
    typedef typename Id<Fragment>::Type TId;

    TId seqId1;
    TSize begin1;
    TId seqId2;
    TSize begin2;
    TSize len;

/*!
 * @fn ExactFragment::Fragment
 * @brief Constructor.
 *
 * @signature Fragment::Fragment();
 * @signature Fragment::Fragment(seqID1, beg1, seqID2, beg2, l);
 *
 * @param[in] seqID1  ID of the first sequence.  Type: <tt>TId</tt>.
 * @param[in] beg1    Begin position of segment match in first sequence.  Type: <tt>TSize</tt>.
 * @param[in] seqID2  ID of the second sequence.  Type: <tt>TId</tt>.
 * @param[in] beg2    Begin position of segment match in second sequence.  Type: <tt>TSize</tt>.
 * @param[in] l       The length of the segment match.  Type: <tt>TSize</tt>.
 */
  
/**
.Memfunc.ExactFragment#Fragment:
..class:Spec.ExactFragment
..summary:Constructor.
..signature:Fragment()
..signature:Fragment(seqId1, beg1, seqId2, beg2, len)
..param.seqId1:The id of the first sequence.
...type:Metafunction.Id
..param.beg1:The TSize begin position on the first sequence.
..param.seqId2:The id of the second sequence.
...type:Metafunction.Id
..param.beg2:The TSize begin position on the second sequence.
..param.len:The TSize length of the segment match.
*/

    Fragment() : seqId1(0), begin1(0), seqId2(0), begin2(0), len(0) {}

    Fragment(TId sqId1, TSize beg1, TId sqId2, TSize beg2, TSize l) :
            seqId1(sqId1), begin1(beg1), seqId2(sqId2), begin2(beg2), len(l) 
    {}

};

template<typename TSize, typename TSpec>
inline bool
operator==(Fragment<TSize, ExactFragment<TSpec> > const & left,
           Fragment<TSize, ExactFragment<TSpec> > const & right)
{
    return (left.seqId1 == right.seqId1 &&
            left.begin1 == right.begin1 &&
            left.seqId2 == right.seqId2 &&
            left.begin2 == right.begin2 &&
            left.len == right.len);
}

template<typename TSize, typename TSpec>
inline bool
operator<(Fragment<TSize, ExactFragment<TSpec> > const & left,
          Fragment<TSize, ExactFragment<TSpec> > const & right)
{
    if (left.seqId1 < right.seqId1)
        return true;
    if (left.seqId1 > right.seqId1)
        return false;
    if (left.begin1 < right.begin1)
        return true;
    if (left.begin1 > right.begin1)
        return false;
    if (left.seqId2 < right.seqId2)
        return true;
    if (left.seqId2 > right.seqId2)
        return false;
    if (left.begin2 < right.begin2)
        return true;
    if (left.begin2 > right.begin2)
        return false;
    if (left.len < right.len)
        return true;
    // if (left.len > right.len)
    //     return false;
    return false;
}

//////////////////////////////////////////////////////////////////////////////
// Exact Fragment that is a forward or reverse match
//////////////////////////////////////////////////////////////////////////////
	

template<typename TSize, typename TSpec>
class Fragment<TSize, ExactReversableFragment<TSpec> > {
public:
    typedef typename Id<Fragment>::Type TId_;

    TId_ seqId1;
    TSize begin1;
    TId_ seqId2;
    TSize begin2;
    TSize len;
    bool reversed;

/*!
 * @fn ExactReversableFragment::Fragment
 * @brief Constructor.
 *
 * @signature Fragment::Fragment();
 * @signature Fragment::Fragment(seqID1, beg1, seqID2, beg2, l[, reversed]);
 *
 * @param[in] seqID1   ID of the first sequence.  Type: <tt>TId</tt>.
 * @param[in] beg1     Begin position of segment match in first sequence.  Type: <tt>TSize</tt>.
 * @param[in] seqID2   ID of the second sequence.  Type: <tt>TId</tt>.
 * @param[in] beg2     Begin position of segment match in second sequence.  Type: <tt>TSize</tt>.
 * @param[in] l        The length of the segment match.  Type: <tt>TSize</tt>.
 * @param[in] reversed A bool; <tt>true</tt> if the segments match in reverse orientation, <tt>false</tt> otherwise.
 */

/**
.Memfunc.ExactReversableFragment#Fragment:
..class:Spec.ExactReversableFragment
..summary:Constructor.
..signature:Fragment()
..signature:Fragment(seqId1, beg1, seqId2, beg2, len[, reversed])
..param.seqId1:The id of the first sequence.
...type:Metafunction.Id
..param.beg1:The TSize begin position on the first sequence.
..param.seqId2:The id of the second sequence.
...type:Metafunction.Id
..param.beg2:The TSize begin position on the second sequence.
..param.len:The TSize length of the segment match.
..param.reversed:$true$ if the segments match in reverse orientation, $false$ otherwise.
...default:$false$
...type:nolink:$bool$
*/
    
    Fragment() : seqId1(0), begin1(0), seqId2(0), begin2(0), len(0), reversed(false) {}
    
    Fragment(TId_ sqId1, TSize beg1, TId_ sqId2, TSize beg2, TSize l) :
            seqId1(sqId1), begin1(beg1), seqId2(sqId2), begin2(beg2), len(l), reversed(false) 
    {}
	
    Fragment(TId_ sqId1, TSize beg1, TId_ sqId2, TSize beg2, TSize l, bool rev) :
            seqId1(sqId1), begin1(beg1), seqId2(sqId2), begin2(beg2), len(l), reversed(rev) 
    {}
};

template<typename TSize, typename TSpec>
inline bool
operator==(Fragment<TSize, ExactReversableFragment<TSpec> > const & left,
           Fragment<TSize, ExactReversableFragment<TSpec> > const & right)
{
    return (left.seqId1 == right.seqId1 &&
            left.begin1 == right.begin1 &&
            left.seqId2 == right.seqId2 &&
            left.begin2 == right.begin2 &&
            left.len == right.len &&
            left.reversed == right.reversed);
}

template<typename TSize, typename TSpec>
inline bool
operator<(Fragment<TSize, ExactReversableFragment<TSpec> > const & left,
          Fragment<TSize, ExactReversableFragment<TSpec> > const & right)
{
    if (left.seqId1 < right.seqId1)
        return true;
    if (left.seqId1 > right.seqId1)
        return false;
    if (left.begin1 < right.begin1)
        return true;
    if (left.begin1 > right.begin1)
        return false;
    if (left.seqId2 < right.seqId2)
        return true;
    if (left.seqId2 > right.seqId2)
        return false;
    if (left.begin2 < right.begin2)
        return true;
    if (left.begin2 > right.begin2)
        return false;
    if (left.len < right.len)
        return true;
    if (left.len > right.len)
        return false;
    if (left.reversed < right.reversed)
        return true;
    // if (left.reversed > right.reversed)
    //     return false;
    return false;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Fragment#label
 * @brief Access to the Fragment's label.
 *
 * @signature TInfix label(frag, stringSet, seqID);
 *
 * @param[in] frag      The Fragment to query.
 * @param[in] stringSet The @link StringSet @endlink with the sequences.
 * @param[in] seqID     The id of the sequence for which the label should be retrieved.
 */

/**
.Function.label
..class:Class.Fragment
..signature:label(f,str,seqId)
..param.f:A fragment.
...type:Class.Fragment
..param.str:The string set underlying the fragment.
..param.seqId:The id of the sequence for which the label should be retrieved.
...remarks:
*/
template<typename TSize, typename TSpec, typename TStringSet, typename TVal>
inline typename Infix<typename Value<TStringSet>::Type>::Type
label(Fragment<TSize, TSpec> const& f,
      TStringSet& str,
      TVal const seqId)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Fragment<TSize, TSpec> >::Type TId;
	return ((TId) seqId == (f.seqId1)) ? infix(getValueById(str, (TId) seqId), f.begin1, f.begin1 + f.len) : infix(getValueById(str, (TId) seqId), f.begin2, f.begin2 + f.len);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Fragment#sequenceId
 * @brief Access to the sequence ID of a fragment.
 *
 * @signature TId sequenceId(frag, seqNum);
 *
 * @param[in] frag   A Fragment.
 * @param[in] seqNum The sequence number fo rwhich the id should be retrieved.  Note that @link Fragment @endlink
                     stores information about exactly two sequences which can be accessed with seqNum 0 or 1 but whose
                     ids may differ from their seqNum.
 *
 * @return TId Reference to the sequence fragment id member.
 */

/**
.Function.sequenceId
..class:Class.Fragment
..signature:sequenceId(f,seqNum)
..param.f:A fragment.
...type:Class.Fragment
..param.seqNum:The sequence number for which the id should be retrieved.
...remarks:Note that @Class.Fragment@ stores information about exactly two sequences which can be accessed with seqNum 0 or 1, but whose ids may differ from their seqNum.
*/
template<typename TSize, typename TSpec, typename TVal>
inline typename Id<Fragment<TSize, TSpec> >::Type &
sequenceId(Fragment<TSize, TSpec> const& f,
		   TVal const seqId)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Fragment<TSize, TSpec> >::Type TId;
	return ((TId) seqId == 0) ? const_cast<TId &>(f.seqId1) : const_cast<TId &>(f.seqId2);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Fragment#fragmentBegin
 * @label Return fragment begin.
 *
 * @signature TSize fragmentBegin(frag, seqId);
 *
 * @param[in] frag  The Fragment to query.
 * @param[in] seqId The id of the sequence to get the begin for.
 *
 * @return TSize Reference to the fragment begin position member.
 */

/**
.Function.fragmentBegin
..class:Class.Fragment
..signature:fragmentBegin(f, seqId)
..param.f:A fragment.
...type:Class.Fragment
..param.seqId:The sequence id for which the begin position should be retrieved.
...remarks:Retrieve with @Function.sequenceId@.
*/
template<typename TSize, typename TSpec, typename TVal>
inline TSize&
fragmentBegin(Fragment<TSize, TSpec> const& f,
			  TVal const seqId)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Fragment<TSize, TSpec> >::Type TId;
	return ((TId) seqId == f.seqId1) ? const_cast<TSize&>(f.begin1) : const_cast<TSize&>(f.begin2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TVal>
inline TSize&
fragmentLength(Fragment<TSize, TSpec> const& f,
			   TVal const)
{
	SEQAN_CHECKPOINT
	return const_cast<TSize&>(f.len);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Fragment#fragmentLength
 * @brief Return fragment begin.
 *
 * @signature TSize fragmentLength(frag);
 *
 * @param[in] frag The Fragment to query for its length.
 *
 * @return TSize Reference to the Fragment's length.
 */

/**
.Function.fragmentLength
..class:Class.Fragment
..signature:fragmentBegin(f)
..param.f:A fragment.
...type:Class.Fragment
*/
template<typename TSize, typename TSpec>
inline TSize&
fragmentLength(Fragment<TSize, TSpec> const& f)
{
	SEQAN_CHECKPOINT
	return const_cast<TSize&>(f.len);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getProjectedPosition
..cat:Alignments
..class:Class.Fragment
..signature:getProjectedPosition(f,seqId,pos,seqId2,pos2)
..summary:Projects a position of one sequence taking part in a pairwise match onto the other sequence.
..signature:getProjectedPosition(f,seqId1,pos1,seqId2,pos2)
..param.f:A fragment.
...type:Class.Fragment
..param.seqId:The id of the sequence to project from.
...type:Metafunction.Id
..param.pos:The position to project.
...type:Metafunction.Size
..param.seqId2:The resulting id of the sequence that pos was projected onto.
...type:Metafunction.Id
..param.pos2:The resulting projected position.
...type:Metafunction.Size
*/
template<typename TSize, typename TSpec, typename TId1, typename TPosition1, typename TId2, typename TPosition2>
inline void
getProjectedPosition(Fragment<TSize, ExactFragment<TSpec> > const& f,
					 TId1 const seqId,
					 TPosition1 const pos,
					 TId2& seqId2,
					 TPosition2& pos2)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Fragment<TSize, TSpec> >::Type TId;
	
	if ((TId) seqId == f.seqId1) {
		SEQAN_ASSERT((TPosition1)f.begin1<=pos);
		SEQAN_ASSERT(pos - f.begin1 < f.len)	;
		pos2 = f.begin2 + (pos - f.begin1);
		seqId2 = f.seqId2;
		return;
	} else {
		SEQAN_ASSERT((TPosition1)f.begin2<=pos);
		SEQAN_ASSERT(pos - f.begin2 < f.len);
		pos2 = f.begin1 + (pos - f.begin2);
		seqId2 = f.seqId1;
		return;
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TValue, typename TId1, typename TPosition1, typename TId2, typename TPosition2>
inline void
getProjectedPosition(Fragment<TSize, ExactFragment<TSpec> > const& f,
					TValue seg_num,
					 TId1 const seqId,
					 TPosition1 const pos,
					 TId2& seqId2,
					 TPosition2& pos2)
{
	(void) seqId;  // When compiled without assertions.
	SEQAN_ASSERT((seg_num == 0 && seqId == f.seqId1) || (seg_num == 1 && seqId == f.seqId2));

	if (seg_num == 0) {
		SEQAN_ASSERT((TPosition1)f.begin1<=pos);
		SEQAN_ASSERT(pos - f.begin1 < f.len)	;
		pos2 = f.begin2 + (pos - f.begin1);
		seqId2 = f.seqId2;
		return;
	} else {
		SEQAN_ASSERT((TPosition1)f.begin2<=pos);
		SEQAN_ASSERT(pos - f.begin2 < f.len);
		pos2 = f.begin1 + (pos - f.begin2);
		seqId2 = f.seqId1;
		return;
	}
}



//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TId1, typename TPosition1, typename TId2, typename TPosition2>
inline void
getProjectedPosition(Fragment<TSize, ExactReversableFragment<TSpec> > const& f,
					 TId1 const seqId,
					 TPosition1 const pos,
					 TId2& seqId2,
					 TPosition2& pos2)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Fragment<TSize, TSpec> >::Type TId;
	
	if ((TId) seqId == f.seqId1) {
		SEQAN_ASSERT((TPosition1)f.begin1<=pos);
		SEQAN_ASSERT(pos - f.begin1 < f.len)	;
		if (f.reversed) pos2 = (f.begin2 + f.len - 1) - (pos - f.begin1);
		else pos2 = f.begin2 + (pos - f.begin1);
		seqId2 = f.seqId2;
		return;
	} else {
		SEQAN_ASSERT((TPosition1)f.begin2<=pos);
		SEQAN_ASSERT(pos - f.begin2 < f.len);
		if (f.reversed) pos2 = (f.begin1 + f.len - 1) - (pos - f.begin2);
		else pos2 = f.begin1 + (pos - f.begin2);
		seqId2 = f.seqId1;
		return;
	}
}


/////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TValue, typename TId1, typename TPosition1, typename TId2, typename TPosition2>
inline void
getProjectedPosition(Fragment<TSize, ExactReversableFragment<TSpec> > const& f,
					 TValue seg_num,
					 TId1 const seqId,
					 TPosition1 const pos,
					 TId2& seqId2,
					 TPosition2& pos2)
{
	SEQAN_CHECKPOINT
	(void) seqId;  // When compiled without assertions.
	SEQAN_ASSERT((seg_num == 0 && seqId==f.seqId1) || (seg_num == 1 && seqId==f.seqId2));

	if (seg_num == 0) {
		SEQAN_ASSERT((TPosition1)f.begin1<=pos);
		SEQAN_ASSERT(pos - f.begin1 < f.len)	;
		if (f.reversed) pos2 = (f.begin2 + f.len - 1) - (pos - f.begin1);
		else pos2 = f.begin2 + (pos - f.begin1);
		seqId2 = f.seqId2;
		return;
	} else {
		SEQAN_ASSERT((TPosition1)f.begin2<=pos);
		SEQAN_ASSERT(pos - f.begin2 < f.len);
		if (f.reversed) pos2 = (f.begin1 + f.len - 1) - (pos - f.begin2);
		else pos2 = f.begin1 + (pos - f.begin2);
		seqId2 = f.seqId1;
		return;
	}
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn ExactReversableFragment#isReversed
 * @brief Return whether a segment match is in reverse orientation.
 *
 * @signature bool isReversed(frag);
 *
 * @param[in] frag The Fragment to query for reverseness.
 *
 * @return bool <tt>true</tt> if the fragment is reversed and <tt>false</tt> otherwise.
 */

/**
.Function.isReversed
..cat:Alignments
..class:Class.Fragment
..signature:isReversed<Fragment<TSize,ExactReversableFragment<TSpec> >(f)
..summary:Returns true if the segment match is in reverse orientation.
..signature:isReversed(f)
..param.f:A fragment.
...type:Class.Fragment
*/
template<typename TSize, typename TSpec>
inline bool
isReversed(Fragment<TSize, ExactReversableFragment<TSpec> > const& f)
{
	SEQAN_CHECKPOINT
	return f.reversed;
}

// Compare lexicographically as tuple.

template<typename TSize, typename TSpec>
inline bool operator>(Fragment<TSize, ExactFragment<TSpec> > const & lhs,
                      Fragment<TSize, ExactFragment<TSpec> > const & rhs)
{
    if (lhs.seqId1 > rhs.seqId1)
        return true;
    if (lhs.seqId1 == rhs.seqId1 && lhs.begin1 > rhs.begin1)
        return true;
    if (lhs.seqId1 == rhs.seqId1 && lhs.begin1 == rhs.begin1 && lhs.seqId2 > rhs.seqId2)
        return true;
    if (lhs.seqId1 == rhs.seqId1 && lhs.begin1 == rhs.begin1 && lhs.seqId2 == rhs.seqId2 && lhs.begin2 > rhs.begin2)
        return true;
    if (lhs.seqId1 == rhs.seqId1 && lhs.begin1 == rhs.begin1 && lhs.seqId2 == rhs.seqId2 && lhs.begin2 == rhs.begin2 && lhs.len > rhs.len)
        return true;
    return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
