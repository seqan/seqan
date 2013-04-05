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

#ifndef SEQAN_HEADER_SHAPE_BASE_H
#define SEQAN_HEADER_SHAPE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	template <unsigned q>
	struct UngappedShape {};
	typedef UngappedShape<0> SimpleShape;

	template <typename TSpec>
	struct GappedShape {};
	typedef GappedShape<Default> GenericShape;


/**
.Class.Shape:
..cat:Index
..summary:Stores hash value and shape for an ungapped or gapped q-gram.
..signature:Shape<TValue, TSpec>
..param.TValue:The @Metafunction.Value@ type of the string the shape is applied to (e.g. $Dna$).
..param.TSpec:The specializing type.
...default:@Spec.SimpleShape@, for ungapped q-grams.
..remarks:The @Metafunction.ValueSize@ of Shape is the ValueSize of TValue which is the alphabet size.
..remarks:To get the span or the weight of a shape call @Function.length@ or @Function.weight@.
.Memfunc.Shape#Shape:
..class:Class.Shape
..summary:Constructor
..signature:Shape<TValue, TSpec> ()
..signature:Shape<TValue, TSpec> (shape)
..param.shape:Other Shape object. (copy constructor)
..include:seqan/index.h
*/
	template <typename TValue = Dna, typename TSpec = SimpleShape>
	class Shape;

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Shape
///.Metafunction.Value.class:Class.Shape
	template <typename TValue, typename TSpec>
	struct Value<Shape<TValue,TSpec> >
	{
		typedef __uint64 Type;
	};

///.Metafunction.Size.param.T.type:Class.Shape
///.Metafunction.Size.class:Class.Shape
	template <typename TValue, typename TSpec>
	struct Size<Shape<TValue,TSpec> >
	{
		typedef unsigned long Type;
	};

///.Metafunction.LENGTH.param.T.type:Class.Shape
///.Metafunction.LENGTH.class:Class.Shape
    template <typename TValue, unsigned q>
	struct LENGTH< Shape<TValue, UngappedShape<q> > >
	{
		enum { VALUE = q };
	};

///.Metafunction.WEIGHT.param.T.type:Class.Shape
///.Metafunction.WEIGHT.class:Class.Shape
    template <typename TValue, unsigned q>
	struct WEIGHT< Shape<TValue, UngappedShape<q> > >
	{
		enum { VALUE = q };
	};

///.Metafunction.ValueSize.param.T.type:Class.Shape
///.Metafunction.ValueSize.class:Class.Shape
	template <typename TValue, typename TSpec>
	struct ValueSize< Shape<TValue, TSpec> > 
	{
		typedef typename Value<Shape<TValue, TSpec> >::Type THashValue;
		static const THashValue VALUE = Power<
						ValueSize<TValue>::VALUE, 
						WEIGHT< Shape<TValue, TSpec> >::VALUE >::VALUE;
	};

///.Metafunction.Host.param.T.type:Class.Shape
///.Metafunction.Host.class:Class.Shape
	template <typename TValue, typename TSpec>
	struct Host<Shape<TValue,TSpec> >
	{
		typedef TValue Type;
	};


//////////////////////////////////////////////////////////////////////////////

/**
.Spec.SimpleShape:
..cat:Index
..summary:A variable length ungapped shape (also called q-gram or k-mer).
..general:Class.Shape
..signature:Shape<TValue, SimpleShape>
..param.TValue:The @Metafunction.Value@ type of the string the shape is applied to (e.g. $Dna$).
..remarks:A SimpleShape must be resized first to a valid length. To do so, call @Function.resize@.
..see:Spec.UngappedShape
..include:seqan/index.h
*/

	//////////////////////////////////////////////////////////////////////////////
	// ungapped shape with variable length
	//////////////////////////////////////////////////////////////////////////////

	template <typename TValue>
	class Shape<TValue, SimpleShape>
	{
	public:
//____________________________________________________________________________

		unsigned					span;
		typename Value<Shape>::Type	hValue;
		typename Value<Shape>::Type	XValue;
		typename Value<Shape>::Type	leftFactor;
		typename Value<Shape>::Type	leftFactor2;
		TValue						leftChar;
//____________________________________________________________________________
		
/**
.Memfunc.SimpleShape#Shape:
..class:Spec.SimpleShape
..summary:Constructor
..signature:Shape<TValue, SimpleShape> ()
..signature:Shape<TValue, SimpleShape> (shape)
..signature:Shape<TValue, SimpleShape> (q)
..param.shape:Other Shape object. (copy constructor)
..param.q:Length of the ungapped q-gram.
*/
		Shape():
			span(0),
			hValue(0),
			XValue(0),
			leftFactor(0),
			leftFactor2(0),
            leftChar(0) {}
		
		Shape(unsigned _span):
			hValue(0),
			XValue(0),
			leftFactor(0),
			leftFactor2(0),
			leftChar(0)
		{
			resize(*this, _span);
		}

		template <unsigned q>
		Shape(Shape<TValue, UngappedShape<q> > const &other)
		{
			*this = other;
		}	

//____________________________________________________________________________

		template <unsigned q>
		inline Shape &
		operator=(Shape<TValue, UngappedShape<q> > const &other)
		{
			span = other.span;
			hValue = other.hValue;
			XValue = other.XValue;
			leftFactor = other.leftFactor;
			leftFactor2 = other.leftFactor2;
			leftChar = other.leftChar;
			return *this;
		}
	};

	//////////////////////////////////////////////////////////////////////////////
	// ungapped shape with fixed length q
	//////////////////////////////////////////////////////////////////////////////

/**
.Spec.UngappedShape:
..cat:Index
..summary:A fixed length ungapped shape (also called q-gram or k-mer).
..general:Class.Shape
..signature:Shape<TValue, UngappedShape<q> >
..param.TValue:The @Metafunction.Value@ type of the sequence the shape is applied to (e.g. $Dna$).
..param.q:The length of the shape.
..include:seqan/index.h
*/

	template <typename TValue, unsigned q>
	class Shape<TValue, UngappedShape<q> >
	{
	public:
		typedef typename Value<Shape>::Type THashValue;
//____________________________________________________________________________

		static const unsigned span = q;
		static const THashValue leftFactor = Power<ValueSize<TValue>::VALUE, q - 1>::VALUE;
		static const THashValue leftFactor2 = (Power<ValueSize<TValue>::VALUE, q>::VALUE - 1) / (ValueSize<TValue>::VALUE - 1);
		// Sigma^(q-1) + Sigma^(q-2) + ... + Sigma + 1

		THashValue	hValue;		// current hash value
		THashValue	XValue;		// Sum_{i=0..q-1} (x_i + 1)
		TValue		leftChar;	// leftmost character
//____________________________________________________________________________
		Shape():
			hValue(0),
			XValue(0),
            leftChar(0) {}
	};



//////////////////////////////////////////////////////////////////////////////

///.Function.value.param.object.type:Class.Shape
///.Function.value.class:Class.Shape
	template <typename TValue, typename TSpec>
	inline typename Value< Shape<TValue, TSpec> >::Type
	value(Shape<TValue, TSpec> &me)
	{
		return me.hValue;
	}

	template <typename TValue, typename TSpec>
	inline typename Value< Shape<TValue, TSpec> >::Type
	value(Shape<TValue, TSpec> const &me)
	{
		return me.hValue;
	}


//____________________________________________________________________________

///.Function.length.param.object.type:Class.Shape
///.Function.length.class:Class.Shape
	template <typename TValue, typename TSpec>
	inline typename Size< Shape<TValue, TSpec> >::Type
	length(Shape<TValue, TSpec> const &me)
	{
	SEQAN_CHECKPOINT
		return me.span;
	}

//____________________________________________________________________________

/**.Function.weight:
..cat:Index
..summary:Number of relevant positions in a shape.
..signature:weight(shape)
..class:Class.Shape
..param.shape:Shape object for which the number of relevant positions is determined.
...type:Class.Shape
..returns:Number of relevant positions.
..remarks.text:For ungapped shapes the return value is the result of the @Function.length@ function.
For gapped shapes this is the number of '1's.
*/
	template <typename TValue, typename TSpec>
	inline typename Size< Shape<TValue, TSpec> >::Type
	weight(Shape<TValue, TSpec> const &me)
	{
	SEQAN_CHECKPOINT
		return length(me);
	}

//____________________________________________________________________________

///.Function.resize.param.object.type:Spec.SimpleShape
///.Function.resize.class:Spec.SimpleShape
	template <typename TValue, typename TSize>
	inline typename Size< Shape<TValue, SimpleShape> >::Type
	resize(Shape<TValue, SimpleShape> & me, TSize new_length)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, SimpleShape> >::Type	THValue;
		me.leftFactor = _intPow((THValue)ValueSize<TValue>::VALUE, new_length - 1);
		me.leftFactor2 = (_intPow((THValue)ValueSize<TValue>::VALUE, new_length) - 1) / (ValueSize<TValue>::VALUE - 1);
		return me.span = new_length;
	}

//____________________________________________________________________________

/**.Function.hash:
..cat:Index
..summary:Computes a (lower) hash value for a shape applied to a sequence.
..signature:hash(shape, it)
..signature:hash(shape, it, charsLeft)
..class:Class.Shape
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the shape.
..param.charsLeft:The distance of $it$ to the string end. 
If $charsLeft$ is smaller than the shape's span, the hash value corresponds to the smallest shape beginning with $charsLeft$ characters.
..returns:Hash value of the shape.
..see:Function.hashNext
..see:Function.hashUpper
..see:Function.hash2
*/

	template <typename TValue, typename TIter>
	inline typename Value< Shape<TValue, SimpleShape> >::Type
	hash(Shape<TValue, SimpleShape> &me, TIter it)
	{
		//typedef typename Value< Shape<TValue, SimpleShape> >::Type	THValue;
		typedef typename Size< Shape<TValue, SimpleShape> >::Type	TSize;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		me.hValue = ordValue(me.leftChar = *it);
		for(TSize i = 1; i < me.span; ++i) {
			++it;
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
		}
		return me.hValue;
	}

	template <typename TValue, typename TIter>
	inline void
	hashInit(Shape<TValue, SimpleShape> &me, TIter it)
	{
		//typedef typename Value< Shape<TValue, SimpleShape> >::Type	THValue;
		typedef typename Size< Shape<TValue, SimpleShape> >::Type	TSize;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        me.leftChar = 0;
		me.hValue = ordValue(*it);
		for(TSize i = 2; i < me.span; ++i) {
			++it;
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
		}
	}

//____________________________________________________________________________
// fixed ungapped shapes

	// loop unrolling ...
	template <typename THValue, typename TValue, typename TIter>
	inline THValue
	_hashFixedShape(THValue hash, TIter &, TValue const, UngappedShape<1> const) {
		return hash;
	}

	template <typename THValue, typename TValue, typename TIter, unsigned q>
	inline THValue
	_hashFixedShape(THValue hash, TIter &it, TValue const, UngappedShape<q> const) {
		++it;
		return _hashFixedShape(
			hash * ValueSize<TValue>::VALUE + ordValue((TValue)*it),
			it, TValue(), UngappedShape<q - 1>());
	}

	// ... for fixed ungapped shapes
	template <typename TValue, unsigned q, typename TIter>
	inline typename Value< Shape<TValue, UngappedShape<q> > >::Type
	hash(Shape<TValue, UngappedShape<q> > &me, TIter it)
	{
		//typedef typename Value< Shape<TValue, UngappedShape<q> > >::Type	THValue;
		//typedef typename Size< Shape<TValue, UngappedShape<q> > >::Type     TSize;

		me.hValue = ordValue(me.leftChar = *it);
		return me.hValue = _hashFixedShape(me.hValue, it, TValue(), UngappedShape<q>());
	}

	template <typename TValue, unsigned q, typename TIter>
	inline typename Value< Shape<TValue, UngappedShape<q> > >::Type
	hashInit(Shape<TValue, UngappedShape<q> > &me, TIter it)
	{
		//typedef typename Value< Shape<TValue, UngappedShape<q> > >::Type	THValue;
		//typedef typename Size< Shape<TValue, UngappedShape<q> > >::Type	TSize;

        me.leftChar = 0;
		me.hValue = ordValue(*it);

        if (q > 1)
            me.hValue = _hashFixedShape(me.hValue, it, TValue(), UngappedShape<q-1>());

		return me.hValue;
    }

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
		//typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = ordValue(me.leftChar = *it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
			}
		} else
			return me.hValue = 0;

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			me.hValue *= ValueSize<TValue>::VALUE;
		return me.hValue;
	}

//____________________________________________________________________________
// Tuple -> fixed ungapped shapes

	template <typename THValue, typename TValue, typename TTValue, unsigned SIZE, typename TPack>
	inline THValue
	_hashTuple2FixedShape(
		THValue const, 
		Tuple<TTValue, SIZE, TPack> const &tuple,
		TValue const,
		UngappedShape<1> const) 
	{
		return ordValue(tuple[0]);
	}

	template <typename THValue, typename TValue, typename TTValue, unsigned SIZE, typename TPack, unsigned q>
	inline THValue
	_hashTuple2FixedShape(
		THValue const, 
		Tuple<TTValue, SIZE, TPack> const &tuple,
		TValue const,
		UngappedShape<q> const) 
	{
		return _hashTuple2FixedShape(THValue(), tuple, TValue(), UngappedShape<q - 1>()) 
			* ValueSize<TValue>::VALUE + ordValue(tuple[q-1]);
	}

	// ... for fixed ungapped shapes
	template <
		typename TValue,
		typename TTValue, 
		unsigned SIZE, 
		unsigned q>
	typename Value< Shape<TValue, UngappedShape<q> > >::Type
	hash(
		Shape<TValue, UngappedShape<q> > &me, 
		Tuple<TTValue, SIZE, BitPacked<> > /*const &*/tuple)
	{
	SEQAN_CHECKPOINT
		if (ValueSize<TValue>::VALUE == (1 << BitsPerValue<TTValue>::VALUE))
			if (q == SIZE)
				return tuple.i;
			else
				return tuple >> (q - SIZE);
		else
			return me.hValue = _hashTuple2FixedShape(me.hValue, tuple, TValue(), UngappedShape<q>());
	}

	template <
		typename TValue,
		typename TTValue, 
		unsigned SIZE, 
		typename TPack, 
		unsigned q>
	typename Value< Shape<TValue, UngappedShape<q> > >::Type
	hash(
		Shape<TValue, UngappedShape<q> > &me, 
		Tuple<TTValue, SIZE, TPack> /*const &*/tuple)
	{
	SEQAN_CHECKPOINT
		return me.hValue = _hashTuple2FixedShape(me.hValue, tuple, TValue(), UngappedShape<q>());
	}

//____________________________________________________________________________

/**.Function.hashUpper:
..cat:Index
..summary:Computes an upper hash value for a shape applied to a sequence.
..signature:hashUpper(shape, it, charsLeft)
..class:Class.Shape
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the shape.
..param.charsLeft:The distance of $it$ to the string end. 
..returns:Upper hash value of the shape.
The hash value corresponds to the maximal @Function.hash@ value of a shape beginning with $min(charsLeft,length(shape))$ characters + 1.
..remarks:This function in conjunction with @Function.hash@ is useful to search a q-gram index for p-grams with p<q.
*/

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hashUpper(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
		//typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = ordValue(me.leftChar = *it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
			}
			++me.hValue;
		} else
			me.hValue = 1;

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			me.hValue *= ValueSize<TValue>::VALUE;
		return me.hValue;
	}

//____________________________________________________________________________

/**
.Function.hashNext:
..cat:Index
..summary:Computes the hash value for the adjacent shape.
..signature:hashNext(shape, it)
..class:Class.Shape
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the adjacent shape.
..returns:Hash value of the q-gram.
..remarks:@Function.hash@ has to be called before.
..include:seqan/index.h
*/

	template <typename TValue, typename TSpec, typename TIter>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hashNext(Shape<TValue, TSpec> &me, TIter const &it)
	{
	SEQAN_CHECKPOINT
		// remove first, shift left, and add next character
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;
		typedef typename Size< Shape<TValue, TSpec> >::Type		TSize;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		me.hValue = 
			(me.hValue - ordValue(me.leftChar) * (THValue)me.leftFactor) * ValueSize<TValue>::VALUE
			+ ordValue((TValue)*(it + ((TSize)me.span - 1)));
		me.leftChar = *it;
		return me.hValue;
	}

//____________________________________________________________________________

/**.Function.hash2:
..cat:Index
..summary:Computes an unique hash value of a shape applied to a sequence, even if the sequence is shorter than the shape span
..signature:hash2(shape, it, charsLeft)
..class:Class.Shape
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the shape.
..param.charsLeft:The distance of $it$ to the string end. 
..returns:Hash value of the shape.
..see:Function.hash2Next
..see:Function.hash2Upper
*/

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash2(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
		//typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = me.XValue = ordValue(me.leftChar = *it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				// update sum of x_i
				me.XValue += ordValue((TValue)*it);
				// shift hash
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + me.XValue;
			}
		} else
			return me.hValue = me.XValue = 0;

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + me.XValue;
		return me.hValue += iEnd;
	}

/**.Function.hash2Upper:
..cat:Index
..summary:Computes an upper unique hash value of a shape applied to a sequence, even if the sequence is shorter than the shape span.
..signature:hash2Upper(shape, it, charsLeft)
..class:Class.Shape
..param.shape:Shape to be used for hashing.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the shape.
..param.charsLeft:The distance of $it$ to the string end. 
..returns:Upper hash value of the shape.
The hash value corresponds to the maximal @Function.hash2@ value of a shape beginning with the $min(charsLeft,length(shape))$ characters + 1
..remarks:This function in conjunction with @Function.hash2@ is useful to search a q-gram index for p-grams with p<q.
*/

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash2Upper(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		THValue hValue, XValue;
		TSize i = 0;
		if (iEnd > 0) {
			hValue = XValue = ordValue((TValue)*it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				// update sum of x_i
				XValue += ordValue((TValue)*it);
				// shift hash
				hValue = hValue * ValueSize<TValue>::VALUE + XValue;
			}
		} else
			hValue = XValue = 0;

		if (charsLeft <= me.span) {
			++XValue;
			++hValue;
		}

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			hValue = hValue * ValueSize<TValue>::VALUE + XValue;
		return hValue += iEnd;
	}

//____________________________________________________________________________

/**
.Function.hash2Next:
..cat:Index
..summary:Computes a unique hash value for the adjacent shape, even if it is shorter than q.
..signature:hash2Next(shape, it)
..class:Class.Shape
..param.shape:Shape to be used for hashing the q-gram.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the adjacent shape.
..returns:Hash value of the shape.
..remarks:@Function.hash@ has to be called before with $shape$ on the left adjacent q-gram.
..include:seqan/index.h
*/

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash2Next(Shape<TValue, TSpec> &me, TIter &it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		// remove first, shift left, and add next character
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		SEQAN_ASSERT_GT((unsigned)me.span, 0u);

		if (charsLeft >= me.span) {
			// update sum of x_i
			me.XValue = me.XValue + ordValue((TValue)*(it + me.span - 1)) - ordValue(me.leftChar);
			// shift hash
			me.hValue = (me.hValue - ordValue(me.leftChar) * (THValue)me.leftFactor2) * ValueSize<TValue>::VALUE + me.XValue
						- me.span * (ValueSize<TValue>::VALUE - 1);
		} else {
			// update sum of x_i
			me.XValue -= ordValue(me.leftChar);
			// shift hash
			me.hValue = (me.hValue - ordValue(me.leftChar) * (THValue)me.leftFactor2) * ValueSize<TValue>::VALUE + me.XValue
				        - charsLeft * (ValueSize<TValue>::VALUE - 1) - ValueSize<TValue>::VALUE;
		}

		me.leftChar = *it;
		return me.hValue;
	}

/**
.Function.unhash:
..cat:Index
..summary:Inverse of the @Function.hash@ function; for ungapped shapes.
..signature:unhash(result, hash, q)
..class:Class.Shape
..param.result:String to write the result to.
...type:Class.String
..param.hash:The hash value previously computed with @Function.hash@.
...type:nolink:A number.
..param.q:The $q$-gram length.
...type:nolink:$unsigned$
..remarks:
..see:Function.hash
..see:Function.hash2
..include:seqan/index.h
*/

	template <typename TString, typename THash>
	inline void unhash(TString &result, THash hash, unsigned q)
	{
	SEQAN_CHECKPOINT
		typedef typename Value<TString>::Type	TValue;

		resize(result, q);
		for (unsigned i = q; i > 0; ) 
		{
			result[--i] = (TValue)(hash % ValueSize<TValue>::VALUE);
			hash /= ValueSize<TValue>::VALUE;
		}
	}

//____________________________________________________________________________

/**.Function.stringToShape:
..cat:Index
..summary:Takes a shape given as a string of '1' (relevant position) and '0' 
(irrelevant position) and converts it into a Shape object.
..signature:stringToShape(shape, bitmap)
..class:Class.Shape
..param.shape:Shape object that is manipulated.
...type:Spec.SimpleShape
..param.bitmap:A character string of '1' and '0' representing relevant and irrelevant positions (blanks) respectively.
...remarks:This string must begin with a '1'. Trailing '0's are ignored.
...remarks:If $shape$ is a @Spec.SimpleShape@ at most one contiguous sequences of $1$s is allowed.
...type:Class.String
..see:Function.shapeToString
..see:Function.reverse
*/

	template <typename TValue, typename TShapeString>
	inline bool
	stringToShape(
		Shape<TValue, SimpleShape> &me, 
		TShapeString const &bitmap)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TShapeString const>::Type		TIter;
		typedef typename Size<TShapeString const>::Type			TSize;

		TIter it = begin(bitmap, Standard());
		TIter itEnd = end(bitmap, Standard());

		TSize ones = 0;
		for(; it != itEnd && *it == '0' ; ++it) ;
		for(; it != itEnd && *it == '1' ; ++it)	++ones;
		for(; it != itEnd && *it == '0' ; ++it) ;

		resize(me, ones);

		return it == itEnd;
	}

//____________________________________________________________________________

/**.Function.shapeToString:
..cat:Index
..class:Class.Shape
..summary:Converts a given shape into a sequence of '1' (relevant position) and '0' 
(irrelevant position).
..signature:shapeToString(bitmap, shape)
..class:Class.Shape
..param.bitmap:The resulting sequence object.
...type:Class.String
..param.shape:Shape object.
...type:Class.Shape
..see:Function.stringToShape
*/

	template <typename TShapeString, typename TValue, unsigned q>
	inline void
	shapeToString(
		TShapeString &bitmap,
		Shape<TValue, UngappedShape<q> > const &me)
	{
	SEQAN_CHECKPOINT

		clear(bitmap);
		resize(bitmap, length(me), '1');
	}

//____________________________________________________________________________
	
///.Function.reverse.param.object.type:Spec.SimpleShape
///.Function.reverse.class:Spec.SimpleShape

	template <typename TValue, typename TSpec>
	inline void
	reverse(Shape<TValue, TSpec> &)
	{
	}
	
}	// namespace seqan

#endif
