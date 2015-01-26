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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_SHAPE_MINIMIZER_H
#define SEQAN_HEADER_SHAPE_MINIMIZER_H

namespace seqan
{

	template <typename TValue, unsigned SPAN, unsigned WEIGHT>
	class Shape<TValue, MinimizerShape<SPAN, WEIGHT> >
	{
	public:
        typedef Shape<TValue, UngappedShape<WEIGHT> >   TMinimizerShape;

        TMinimizerShape minShape;
	};

    template <typename TValue, unsigned q>
	struct LENGTH< Shape<TValue, UngappedShape<q> > >
	{
		enum { VALUE = q };
	};

    template <typename TValue, unsigned q>
	struct WEIGHT< Shape<TValue, UngappedShape<q> > >
	{
		enum { VALUE = q };
	};


	template <typename TValue, unsigned SPAN, unsigned WEIGHT, typename TIter>
    inline typename Value< Shape<TValue, UngappedShape<q> > >::Type
	hash(Shape<TValue, MinimizerShape<SPAN, WEIGHT> > &me, TIter it)
	{
        TMinimizerIt minIt = getMinimizer(it);

        return hash(me.minShape, minIt);
	}

    template <typename TValue, typename TSpec>
	inline SEQAN_HOST_DEVICE
	typename Size< Shape<TValue, TSpec> >::Type
	length(Shape<TValue, TSpec> const &me)
	{
	SEQAN_CHECKPOINT
		return me.span;
	}

	template <typename TValue, typename TSpec>
	inline SEQAN_HOST_DEVICE
    typename Size< Shape<TValue, TSpec> >::Type
	weight(Shape<TValue, TSpec> const &me)
	{
	SEQAN_CHECKPOINT
		return length(me);
	}

}	// namespace seqan

#endif
