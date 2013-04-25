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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// TODO(holtgrew): Switch to Host interface.

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GAPS_BASE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GAPS_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TSpec>
struct GapsIterator;

struct ArrayGaps_;
typedef Tag<ArrayGaps_> ArrayGaps;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Gaps
// ----------------------------------------------------------------------------

/**
.Class.Gaps
..cat:Alignments
..implements:Concept.SequenceConcept
..summary:Efficient storage of gaps for a sequence.
..signature:Gaps<TSequence, TSpec>
..description:
Gaps wrap a @Concept.SequenceConcept@ and allows to (1) insert gaps into the sequence and (2) select an infix of the gapped sequence (clipping).
The gaps are not inserted into the underlying sequence (source) but stored separately.
Using the clipping is optional and meant for selecting parts of the alignment as a part of the result of a local alignment algorithm.
..description.image:gaps_illustration|Illustration of Gaps object and positions with clipping.
..description:
In the figure above, the source sequence has seven characters, the gapped sequence has four gaps and thus consists of eleven characters.
The gapped sequence is clipped to start at position 0 in the gapped sequence and to end at position 8 in the gapped sequence (the positions given as half-open intervals $[begin, end)$).
..description:
The figure shows the three coordinate systems that are used with Gaps objects.
The source position is the position in the underlying sequence.
The unclipped view position is the position in the gapped sequence without gaps.
The view position is the position in the gapped sequence but including the clipping:
All (clipped) view positions have the clipping begin position subtracted from them.
..description:

..param.TSequence:The type of the underlying sequence.
...type:Concept.SequenceConcept
..param.TSpec:Specialization tag.
..include:seqan/align.h
 */

template <typename TSequence, typename TSpec = ArrayGaps>
class Gaps;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

///.Metafunction.Value.param.T.type:Class.Gaps
///.Metafunction.Value.class:Class.Gaps

template <typename TSequence, typename TSpec>
struct Value<Gaps<TSequence, TSpec> >
{
    typedef typename Value<TSequence>::Type           TAlphabet;
    typedef typename GappedValueType<TAlphabet>::Type Type;
};

template <typename TSequence, typename TSpec>
struct Value<Gaps<TSequence, TSpec> const> : Value<Gaps<TSequence, TSpec> >
{};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

///.Metafunction.Iterator.param.T.type:Class.Gaps
///.Metafunction.Iterator.class:Class.Gaps

template <typename TSequence, typename TSpec, typename TIteratorSpec>
struct Iterator<Gaps<TSequence, TSpec>, TIteratorSpec>
{
	typedef Iter<Gaps<TSequence, TSpec>, GapsIterator<TSpec> > Type;
};

template <typename TSequence, typename TSpec, typename TIteratorSpec>
struct Iterator<Gaps<TSequence, TSpec> const, TIteratorSpec>
{
	typedef Iter<Gaps<TSequence, TSpec> const, GapsIterator<TSpec> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

///.Metafunction.GetValue.param.T.type:Class.Gaps
///.Metafunction.GetValue.class:Class.Gaps

template <typename TSequence, typename TSpec>
struct GetValue<Gaps<TSequence, TSpec> > : Value<Gaps<TSequence, TSpec> >
{};

template <typename TSequence, typename TSpec>
struct GetValue<Gaps<TSequence, TSpec> const> : GetValue<Gaps<TSequence, TSpec> >
{};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

///.Metafunction.Position.param.T.type:Class.Gaps
///.Metafunction.Position.class:Class.Gaps

template <typename TSequence, typename TSpec>
struct Position<Gaps<TSequence, TSpec> >
{
    typedef typename Position<TSequence>::Type TSeqPos_;
    typedef typename MakeSigned<TSeqPos_>::Type Type;
};

template <typename TSequence, typename TSpec>
struct Position<Gaps<TSequence, TSpec> const> : Position<Gaps<TSequence, TSpec> >
{};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

///.Metafunction.Reference.param.T.type:Class.Gaps
///.Metafunction.Reference.class:Class.Gaps

template <typename TSequence, typename TSpec>
struct Reference<Gaps<TSequence, TSpec> >
{
	typedef typename Iterator<Gaps<TSequence, TSpec>, Standard>::Type TIterator_;
	typedef Proxy<IteratorProxy<TIterator_> > Type;
};

template <typename TSequence, typename TSpec>
struct Reference<Gaps<TSequence, TSpec> const>
{
	typedef typename Iterator<Gaps<TSequence, TSpec> const, Standard>::Type TIterator_;
	typedef Proxy<IteratorProxy<TIterator_> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

///.Metafunction.Size.param.T.type:Class.Gaps
///.Metafunction.Size.class:Class.Gaps

template <typename TSequence, typename TSpec>
struct Size<Gaps<TSequence, TSpec> >
{
    typedef typename Size<TSequence>::Type Type;
};

template <typename TSequence, typename TSpec>
struct Size<Gaps<TSequence, TSpec> const> : Size<Gaps<TSequence, TSpec> >
{};

// ----------------------------------------------------------------------------
// Metafunction Source
// ----------------------------------------------------------------------------

// TODO(holtgrew): Switch to Hosted Type interface

/**
.Metafunction.Source
..cat:Alignments
..summary:Return underlying sequence of Gaps/Alignments.
..signature:Source<T>::Type
..param.T:The type to query for underlying sequence.
..include:seqan/align.h
*/

///.Metafunction.Source.param.T.type:Class.Gaps
///.Metafunction.Source.class:Class.Gaps

template <typename TSequence, typename TSpec>
struct Source<Gaps<TSequence, TSpec> >
{
    typedef TSequence Type;
};

template <typename TSequence, typename TSpec>
struct Source<Gaps<TSequence, TSpec> const> : Source<Gaps<TSequence, TSpec> >
{};

// TODO(holtgrew): Also prefix/suffix/infix? Should work!

// ----------------------------------------------------------------------------
// Metafunction IsSequence
// ----------------------------------------------------------------------------

///.Metafunction.IsSequence.param.T.type:Class.Gaps
///.Metafunction.IsSequence.class:Class.Gaps

template <typename TSequence, typename TSpec>
struct IsSequence<Gaps<TSequence, TSpec> >
{
    typedef True Type;
    static const bool VALUE = true;
};

template <typename TSequence, typename TSpec>
struct IsSequence<Gaps<TSequence, TSpec> const> : IsSequence<Gaps<TSequence, TSpec> >
{};


// ============================================================================
// Functions
// ============================================================================

// TODO(holtgrew): Adding links to implemented sequence. This should be cleaned up once we have better documentation with concepts.
///.Function.begin.class:Class.Gaps
///.Function.end.class:Class.Gaps
///.Function.iter.class:Class.Gaps

// ----------------------------------------------------------------------------
// Function setSource()
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function createSource()
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function clearClipping()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#clearClipping
..class:Class.Gaps
..summary:Clear clipping from @Class.Gaps@ object.
..cat:Alignments
..signature:void clearClipping(gaps)
..param.gaps:The @Class.Gaps@ object to clear.
...type:Class.Gaps
..returns:$void$
..see:Function.Gaps#clearGaps
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function clearGaps()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#clearGaps
..class:Class.Gaps
..summary:Clear gaps and clipping from @Class.Gaps@ object.
..cat:Alignments
..signature:void clearGaps(gaps)
..param.gaps:The @Class.Gaps@ object to clear.
...type:Class.Gaps
..returns:$void$
..see:Function.Gaps#clearClipping
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#length
..class:Class.Gaps
..summary:Return length of the gapped sequence.
..cat:Alignments
..signature:TSize length(gaps)
..param.gaps:The @Class.Gaps@ object to query.
...type:Class.Gaps
..returns:Length of the gapped sequence.
...type:Metafunction.Size
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function unclippedLength()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#unclippedLength
..class:Class.Gaps
..summary:Return length of the gapped sequence without clipping.
..cat:Alignments
..signature:TSize unclippedLength(gaps)
..param.gaps:The @Class.Gaps@ object to query.
...type:Class.Gaps
..returns:Length of the gapped sequence, ignoring any clipping.
...type:Metafunction.Size
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function toViewPosition()
// ----------------------------------------------------------------------------

/**
.Function.toViewPosition
..class:Class.Gaps
..summary:Transforms source to view position.
..cat:Alignments
..signature:toViewPosition(gaps, pos)
..param.gap:A Gaps object, e.g. a row in the alignment.
...type:Class.Gaps
..param.pos:Position in the original sequence to get the view position for.
..returns:The position in the view/gaps position.
..remarks:If $gap$ is a clipped alignment row, gaps in the clipped part will influence the result. The position $pos$ is counted from the unclipped begin position and must be greater or equal the clipped begin position of $gap$.
..see:Function.toSourcePosition
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function toSourcePosition()
// ----------------------------------------------------------------------------

/**
.Function.toSourcePosition
..class:Class.Gaps
..summary:Transforms view to source position, if the view position is a gap, the original position of the next non-gap entry is returned.
..cat:Alignments
..signature:toSourcePosition(gaps, pos)
..param.gap:A Gaps object, e.g. a row in the alignment.
...type:Class.Gaps
..param.pos:Position in the view sequence (this includes gaps) to get the original position for.
..returns:The position in the source sequence.
..remarks:If $gap$ is a clipped alignment row, gaps in the clipped part will influence the result. The position $pos$ is counted from the unclipped begin position.
..see:Function.toViewPosition
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function isGap()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#isGap
..class:Class.Gaps
..summary:Query whether a given clipped view position is a gap.
..cat:Alignments
..signature:bool isGap(gaps, clippedViewPos)
..param.gaps:The @Class.Gaps@ object to query.
...type:Class.Gaps
..param.clippedViewPos:The position in the view to query.
...type:Metafunction.Position
..returns:Whether or not there is a gap at the given clipped view position.
...type:nolink:$bool$
..see:Function.insertGap
..see:Function.removeGap
..see:Function.removeGaps
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function insertGaps()
// ----------------------------------------------------------------------------

/**
.Function.insertGaps
..class:Class.Gaps
..summary:Insert multiple gaps into a gapped sequence.
..cat:Alignments
..signature:void insertGaps(gaps, clippedViewPos, count)
..param.gaps:The @Class.Gaps@ object to insert gaps into.
...type:Class.Gaps
..param.clippedViewPos:The position in the view to insert gaps at.
...type:Metafunction.Position
..param.count:The number of gaps to insert
...type:nolink:$unsigned$
..returns:$void$
..see:Function.insertGap
..see:Function.removeGap
..see:Function.removeGaps
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function insertGap()
// ----------------------------------------------------------------------------

/**
.Function.insertGap
..class:Class.Gaps
..summary:Insert one gap into a gapped sequence.
..cat:Alignments
..signature:void insertGap(gaps, clippedViewPos)
..param.gaps:The @Class.Gaps@ object to insert gap into.
...type:Class.Gaps
..param.clippedViewPos:The position in the view to insert gap at.
...type:Metafunction.Position
..returns:$void$
..see:Function.insertGaps
..see:Function.removeGap
..see:Function.removeGaps
..include:seqan/align.h
*/

// Forward to removeGaps() which has to be implemented in each subclass.

template <typename TSequence, typename TSpec, typename TPosition>
inline void
insertGap(Gaps<TSequence, TSpec> & gaps, TPosition clippedViewPos)
{
    insertGaps(gaps, clippedViewPos, 1u);
}

// ----------------------------------------------------------------------------
// Function removeGaps()
// ----------------------------------------------------------------------------

/**
.Function.removeGaps
..class:Class.Gaps
..summary:Remove multiple gaps from a gapped sequence.
..cat:Alignments
..signature:TSize removeGaps(gaps, clippedViewPos, count)
..param.gaps:The @Class.Gaps@ object to remove gaps into.
...type:Class.Gaps
..param.clippedViewPos:The position in the view to remove gaps from.
...type:Metafunction.Position
..param.count:The number of gaps to remove
...type:nolink:$unsigned$
..returns:The number of removed gaps.
...type:Metafunction.Size
..see:Function.insertGap
..see:Function.insertGaps
..see:Function.removeGap
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function removeGap()
// ----------------------------------------------------------------------------

/**
.Function.removeGap
..class:Class.Gaps
..summary:Remove one gap from a gapped sequence.
..cat:Alignments
..signature:TSize removeGap(gaps, clippedViewPos)
..param.gaps:The @Class.Gaps@ object to remove gap into.
...type:Class.Gaps
..param.clippedViewPos:The position in the view to remove gap from.
...type:Metafunction.Position
..returns:The number of removed gaps.
...type:Metafunction.Size
..see:Function.insertGap
..see:Function.insertGaps
..see:Function.removeGaps
..include:seqan/align.h
*/

// Forward to removeGaps() which has to be implemented in each subclass.

template <typename TSequence, typename TSpec, typename TPosition>
inline typename Size<Gaps<TSequence, TSpec> >::Type
removeGap(Gaps<TSequence, TSpec> & gaps, TPosition clippedViewPos)
{
    return removeGaps(gaps, clippedViewPos, 1u);
}
    
// ----------------------------------------------------------------------------
// Function countGaps()
// ----------------------------------------------------------------------------
    
/**
.Function.Gaps#countGaps
..class:Class.Gaps
..summary:Reports number of continues gaps right of current iterator position.
..cat:Alignments
..signature:TSize countGaps(iter)
..param.iter:Iterator of the  @Class.Gaps@ object to count gaps for.
...type:Metafunction.Iterator
..returns:The number of gaps right of the current iterator position, including the current position, or $0$ if there is no gap.
...type:Metafunction.Size
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function setClippedBeginPosition()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#setClippedBeginPosition
..class:Class.Gaps
..summary:Sets the begin position of the clipping.
..signature:void setClippedBeginPosition(gaps, unclippedViewPosition)
..param.gaps:The @Class.Gaps@ object to query.
...type:Class.Gaps
..param.unclippedViewPosition:The position in the unclipped view to set as the clipping begin position.
...type:Metafunction.Position
..returns:$void$
..remarks:Note that the position is *not* a clipped view position but an uncliped one!
..see:Function.Gaps#beginPosition
..see:Function.Gaps#endPosition
..see:Function.Gaps#setBeginPosition
..see:Function.Gaps#setEndPosition
..see:Function.Gaps#clippedBeginPosition
..see:Function.Gaps#clippedEndPosition
..see:Function.Gaps#setClippedEndPosition
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function setClippedEndPosition()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#setClippedEndPosition
..class:Class.Gaps
..summary:Sets the end position of the clipping.
..signature:void setClippedEndPosition(gaps, unclippedViewPosition)
..param.gaps:The @Class.Gaps@ object to query.
...type:Class.Gaps
..param.unclippedViewPosition:The position in the unclipped view to set as the clipping end position.
...type:Metafunction.Position
..returns:$void$
..remarks:Note that the position is *not* a clipped view position but an uncliped one!
..see:Function.Gaps#beginPosition
..see:Function.Gaps#endPosition
..see:Function.Gaps#setBeginPosition
..see:Function.Gaps#setEndPosition
..see:Function.Gaps#clippedBeginPosition
..see:Function.Gaps#clippedEndPosition
..see:Function.Gaps#setClippedBeginPosition
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function clippedBeginPosition()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#clippedBeginPosition
..class:Class.Gaps
..summary:Return the begin position of the clipping in the unclipped gapped sequence.
..signature:TPosition clippedBeginPosition(gaps)
..param.gaps:The @Class.Gaps@ object to query.
...type:Class.Gaps
..returns:The begin position of the current clipped view in the unclipped gapped sequence.
...type:Metafunction.Position
..see:Function.Gaps#beginPosition
..see:Function.Gaps#endPosition
..see:Function.Gaps#setBeginPosition
..see:Function.Gaps#setEndPosition
..see:Function.Gaps#clippedEndPosition
..see:Function.Gaps#setClippedBeginPosition
..see:Function.Gaps#setClippedEndPosition
..example:
In the following gaps configuration, the result of $clippedBeginPosition(gaps)$ is $1$.
..example.code:
clipping                   [     )
  (half-open interval)           

gapped sequence:          X--XXX-XX-

source position:          0111234456
unclipped view position:  0123456789
clipped view position:     0123456
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function clippedEndPosition()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#clippedEndPosition
..class:Class.Gaps
..summary:Return the end position of the clipping in the unclipped gapped sequence.
..signature:TPosition clippedEndPosition(gaps)
..param.gaps:The @Class.Gaps@ object to query.
...type:Class.Gaps
..returns:The end position of the current clipped view in the unclipped gapped sequence.
...type:Metafunction.Position
..see:Function.Gaps#beginPosition
..see:Function.Gaps#endPosition
..see:Function.Gaps#setBeginPosition
..see:Function.Gaps#setEndPosition
..see:Function.Gaps#clippedBeginPosition
..see:Function.Gaps#setClippedBeginPosition
..see:Function.Gaps#setClippedEndPosition
..example:
In the following gaps configuration, the result of $clippedEndPosition(gaps)$ is $7$.
..example.code:
clipping                   [     )
  (half-open interval)           

gapped sequence:          X--XXX-XX-

source position:          0111234456
unclipped view position:  0123456789
clipped view position:     0123456
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function setBeginPosition()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#setBeginPosition
..class:Class.Gaps
..summary:Set the begin position of the clipped gapped sequence, given a source position.
..signature:void setBeginPosition(gaps, sourcePosition)
..param.gaps:The @Class.Gaps@ object to query.
...type:Class.Gaps
..param.sourcePosition:The source position to set the clipping begin to.
...type:Metafunction.Position
..returns:$void$
..see:Function.Gaps#beginPosition
..see:Function.Gaps#endPosition
..see:Function.Gaps#setEndPosition
..see:Function.Gaps#clippedBeginPosition
..see:Function.Gaps#clippedEndPosition
..see:Function.Gaps#setClippedBeginPosition
..see:Function.Gaps#setClippedEndPosition
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function setEndPosition()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#setEndPosition
..class:Class.Gaps
..summary:Set the end position of the clipped gapped sequence, given a source position.
..signature:void setEndPosition(gaps, sourcePosition)
..param.gaps:The @Class.Gaps@ object to query.
...type:Class.Gaps
..param.sourcePosition:The source position to set the clipping end to.
...type:Metafunction.Position
..returns:$void$
..see:Function.Gaps#beginPosition
..see:Function.Gaps#endPosition
..see:Function.Gaps#setBeginPosition
..see:Function.Gaps#clippedBeginPosition
..see:Function.Gaps#clippedEndPosition
..see:Function.Gaps#setClippedBeginPosition
..see:Function.Gaps#setClippedEndPosition
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function beginPosition()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#beginPosition
..class:Class.Gaps
..summary:Return the clipping begin position as a source position.
..signature:TPosition beginPosition(gaps)
..param.gaps:The @Class.Gaps@ object to query.
...type:Class.Gaps
..returns:The begin position of the current clipped view in the source.
...type:Metafunction.Position
..see:Function.Gaps#endPosition
..see:Function.Gaps#setBeginPosition
..see:Function.Gaps#setEndPosition
..see:Function.Gaps#clippedBeginPosition
..see:Function.Gaps#clippedEndPosition
..see:Function.Gaps#setClippedBeginPosition
..see:Function.Gaps#setClippedEndPosition
..example:
In the following gaps configuration, the result of $beginPosition(gaps)$ is $1$.
The clipping starts in a gap and the source position of the first non-gap character right of the clipping begin has source position $1$.
..example.code:
clipping                   [     )
  (half-open interval)           

gapped sequence:          X--XXX-XX-

source position:          0111234456
unclipped view position:  0123456789
clipped view position:     0123456
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function endPosition()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#endPosition
..class:Class.Gaps
..summary:Return the clipping end position as a source position.
..signature:TPosition endPosition(gaps)
..param.gaps:The @Class.Gaps@ object to query.
...type:Class.Gaps
..returns:The end position of the current clipped view in the source.
...type:Metafunction.Position
..see:Function.Gaps#beginPosition
..see:Function.Gaps#setBeginPosition
..see:Function.Gaps#setEndPosition
..see:Function.Gaps#clippedBeginPosition
..see:Function.Gaps#clippedEndPosition
..see:Function.Gaps#setClippedBeginPosition
..see:Function.Gaps#setClippedEndPosition
..example:
In the following gaps configuration, the result of $endPosition(gaps)$ is $4$.
..example.code:
clipping                   [     )
  (half-open interval)           

gapped sequence:          X--XXX-XX-

source position:          0111234456
unclipped view position:  0123456789
clipped view position:     0123456
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function write()
// ----------------------------------------------------------------------------

template <typename TFile, typename TSource, typename TIDString, typename TSpec>
inline void
write(TFile & target,
	  Gaps<TSource, TSpec> const & source, 
	  TIDString const &,
	  Raw)
{
//IOREV _nodoc_ specialization not documented

	// Print gaps row
	typedef typename Iterator<Gaps<TSource, TSpec> const>::Type TIter;
	TIter begin_ = begin(source);
	TIter end_ = end(source);
	for (; begin_ != end_; ++begin_) {
		if (isGap(begin_))
			streamPut(target, gapValue<char>());
		else 
			streamPut(target, convert<char>(*begin_));
	}
}

// ----------------------------------------------------------------------------
// Function operator<<()                                      [stream operator]
// ----------------------------------------------------------------------------

// TODO(holtgrew): Document appropriately.

template <typename TStream, typename TSource, typename TSpec>
inline TStream &
operator<<(TStream & stream, Gaps<TSource, TSpec> const & gaps)
{
    typedef Gaps<TSource, TSpec> const             TGaps;
    typedef typename Iterator<TGaps, Rooted>::Type TIter;

    for (TIter it = begin(gaps, Rooted()); !atEnd(it); goNext(it))
    {
        // TODO(holtgrew): Ideally, we could simply print the expanded alphabet char but that is broken.
        if (isGap(it))
            stream << gapValue<char>();
        else
            stream << convert<char>(*it);
    }

    return stream;
}

// ----------------------------------------------------------------------------
// Function _pumpTraceToGaps()
// ----------------------------------------------------------------------------

// Internal function for converting AlignTrace<> objects into alignments in two Gaps objects.  Note that the traceback
// in the trace is stored in reverse, from back to front.  We insert the gaps in descending order of their position.
// The reason is that Gaps<> objects store the gaps in ascending order of coordinates in String<> objects and inserting
// at the end is in O(1) while inserting in the front is O(n).

template <typename TSequenceH, typename TGapsSpecH, typename TSequenceV, typename TGapsSpecV, typename TSize>
void _pumpTraceToGaps(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                      Gaps<TSequenceV, TGapsSpecV> & gapsV,
                      AlignTraceback<TSize> const & trace)
{
    typedef Gaps<TSequenceH, TGapsSpecH> TGapsH;
    typedef typename Iterator<TGapsH, Standard>::Type TGapsHIter;

    typedef Gaps<TSequenceV, TGapsSpecV> TGapsV;
    typedef typename Iterator<TGapsV, Standard>::Type TGapsVIter;

    // TODO(holtgrew): I don't understand the following.  Originally, this function used Align objects, but I did not understand it there either.
	// TODO(rausch): Pump trace into align_ (note: this is relatively slow code here. it could be improved if specialized to the Align Specs).
    clearGaps(gapsH);
    clearClipping(gapsH);
    clearGaps(gapsV);
    clearClipping(gapsV);

	TSize i = length(trace.sizes);  // Scan trace backwards.
	TGapsHIter itH = begin(gapsH);
	TGapsVIter itV = begin(gapsV);
	while (i > 0)
	{
		--i;
		TSize size = trace.sizes[i];
		switch ((int) trace.tvs[i])
		{
		case 1:  // Go horizontal.
			insertGaps(itV, size);
			break;

		case 2:  // Go vertical.
			insertGaps(itH, size);
			break;
		}
		goFurther(itH, size);
		goFurther(itV, size);
	}
}

// ----------------------------------------------------------------------------
// Function source()
// ----------------------------------------------------------------------------

/**
.Function.source
..summary:Return underlying object.
..cat:Basic
..signature:source(obj)
..param.obj:The object to get underlying sequence of.
...type:Class.Gaps
..returns:The underlying object.
...type:Metafunction.Source
..include:seqan/align.h
*/

/*
.Function.Gaps#source
..summary:Return underlying sequence.
..cat:Alignments
..signature:TSource source(gaps)
..param.obj:The object to get underlying sequence of.
...type:Class.Gaps
..returns:Reference to the underlying sequence.
...type:Metafunction.Source
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function sourceSegment()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Rename/remove?

// We need some forwards for function sourceSegment().

template <typename TSequence, typename TSpec>
inline typename Position<Gaps<TSequence, TSpec> >::Type clippedBeginPosition(Gaps<TSequence, TSpec> const & gaps);
template <typename TSequence, typename TSpec>
inline typename Position<Gaps<TSequence, TSpec> >::Type clippedEndPosition(Gaps<TSequence, TSpec> const & gaps);
template <typename TSequence, typename TSpec, typename TPosition>
inline typename Position<TSequence>::Type
toSourcePosition(Gaps<TSequence, TSpec> const & gaps, TPosition clippedViewPos);

template <typename TSequence, typename TSpec>
inline typename Infix<TSequence>::Type
sourceSegment(Gaps<TSequence, TSpec> const & gaps)
{
    return infix(source(gaps), toSourcePosition(gaps, clippedBeginPosition(gaps)), toSourcePosition(gaps, clippedEndPosition(gaps)));
}

template <typename TSequence, typename TSpec>
inline typename Infix<TSequence>::Type
sourceSegment(Gaps<TSequence, TSpec> & gaps)
{
    return infix(source(gaps), toSourcePosition(gaps, clippedBeginPosition(gaps)), toSourcePosition(gaps, clippedEndPosition(gaps)));
}

// ----------------------------------------------------------------------------
// Function assignSource()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#assignSource
..class:Class.Gaps
..summary:Assign the source of a Gaps object, copying data.
..cat:Alignments
..signature:void assignSource(gaps, sequence)
..param.gaps:The @Class.Gaps@ object to assign the source of.
...type:Class.Gaps
..param.sequence:The @Concept.SequenceConcept@ to assign as the source.
...type:Metafunction.Source
..remarks:This will copy $sequence$ into the source of $gaps$.
..returns:$void$
..see:Function.Gaps#setSource
..see:Function.source
..include:seqan/align.h
*/

// TOOD(holtgrew): Switch to Hosted Type?

template <typename TSequence, typename TSpec, typename TValue>
inline void
assignSource(Gaps<TSequence, TSpec> & gaps, TValue const & value)
{
    assign(source(gaps), value);
}

// ----------------------------------------------------------------------------
// Function setSource()
// ----------------------------------------------------------------------------

/**
.Function.Gaps#setSource
..class:Class.Gaps
..summary:Set the source of a Gaps object, do not copy if possible.
..cat:Alignments
..signature:void setSource(gaps, sequence)
..param.gaps:The @Class.Gaps@ object to set the source of.
...type:Class.Gaps
..param.sequence:The @Concept.SequenceConcept@ to set as the source.
...type:Metafunction.Source
..remarks:This will avoid copying if possible.
..returns:$void$
..see:Function.Gaps#assignSource
..see:Function.source
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

///.Function.clear.class:Class.Gaps

template <typename TSequence, typename TSpec>
inline void clearGaps(Gaps<TSequence, TSpec> & gaps);
template <typename TSequence, typename TSpec>
inline void clearClipping(Gaps<TSequence, TSpec> & gaps);

template <typename TSequence, typename TSpec>
inline void clear(Gaps<TSequence, TSpec> & gaps)
{
    clearGaps(gaps);
    clearClipping(gaps);
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
inline bool operator==(Gaps<TSequence, TSpec> const & lhs,
                       Gaps<TSequence, TSpec> const & rhs)
{
    typename Comparator<Gaps<TSequence, TSpec> >::Type lex(lhs, rhs);
    return isEqual(lex);
}

template <typename TSequence, typename TSpec, typename TRightHandSide>
inline bool operator==(Gaps<TSequence, TSpec> const & lhs,
                       TRightHandSide const & rhs)
{
    typename Comparator<Gaps<TSequence, TSpec> >::Type lex(lhs, rhs);
    return isEqual(lex);
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
inline bool operator!=(Gaps<TSequence, TSpec> const & lhs,
                       Gaps<TSequence, TSpec> const & rhs)
{
    return !(lhs == rhs);
}


}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GAPS_BASE_H_
