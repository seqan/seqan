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

#ifndef SEQAN_HEADER_SHAPE_MINIMIZER_H
#define SEQAN_HEADER_SHAPE_MINIMIZER_H

namespace seqan
{

struct ReverseComplement_;
typedef Tag<ReverseComplement_> const   ReverseComplement;

// ----------------------------------------------------------------------------
// Struct MinimizerShape
// ----------------------------------------------------------------------------

template <unsigned TSPAN, unsigned TWEIGHT, typename TSpec = void>
struct MinimizerShape;

// ----------------------------------------------------------------------------
// Class Shape<MinimizerShape>
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
class Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> >
{
public:
    typedef typename Value<Shape>::Type THashValue;

    unsigned span;
    unsigned weight;
    THashValue hValue;      //current minimizer hash

    Shape():
        span(TSPAN),
        weight(TWEIGHT)
    {}
};

// ----------------------------------------------------------------------------
// Metafunction LENGTH
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
struct LENGTH<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >
{
    enum { VALUE = TSPAN };
};

// ----------------------------------------------------------------------------
// Metafunction WEIGHT
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
struct WEIGHT<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >
{
    enum { VALUE = TWEIGHT };
};

// ----------------------------------------------------------------------------
// Function weight()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
inline SEQAN_HOST_DEVICE
typename Size< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
weight(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > const &me)
{
    return me.weight;
}

//return  lexicographically smaller of S as the minimizer
template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter> 
inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type 
_getMinimizer(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    typedef typename Value<Shape<TValue, UngappedShape<TWEIGHT> > >::Type THValue;
  
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    SEQAN_ASSERT_GT((unsigned)me.span, (unsigned)me.weight); 

    Shape<TValue, UngappedShape<TWEIGHT> > tmpShape;
    TIter leftIt = it;
    THValue miniTmp = hash(tmpShape, leftIt);
    for (unsigned i = 1; i < me.span - me.weight + 1; i++)
    {
        if(miniTmp > hashNext(tmpShape, leftIt + i))
            miniTmp = tmpShape.hValue;        
    }
    me.hValue = miniTmp;
    return me.hValue;

}

//return  lexicographically smaller of S and the reverse complement of S as the minimizer 
template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TIter>
inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > >::Type
_getMinimizer(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, ReverseComplement> > &me, TIter const &it)
{
    typedef typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, ReverseComplement> > >::Type THValue;

    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    SEQAN_ASSERT_GT((unsigned)me.span, (unsigned)me.weight); 
    TIter leftIt = it;
    Shape<TValue, UngappedShape<TWEIGHT> > tmpShape;
    hashInit(tmpShape, leftIt);
    THValue miniTmp = hash(tmpShape, leftIt);

    for (int i = 1 ; i < tmpShape.span; i++)
    {
        if(miniTmp > hashNext(tmpShape, leftIt + i))
            miniTmp = tmpShape.hValue;        
    }

    String<TValue> rc;
    resize(rc, leftIt -it + me.span - 1 );
    TIter rcIt = begin(rc);
    arrayCopy(it, leftIt + length(me), rcIt);
    reverseComplement(rc);
    leftIt = begin(rc);

    THValue miniTmpRC = hash(tmpShape, leftIt);
    for (int i = 1 ; i < tmpShape.span; i++)
    {
        if(miniTmpRC > hashNext(tmpShape, leftIt + i))
            miniTmpRC = me.hValue;
    }
    me.hValue = (miniTmp < miniTmpRC ? miniTmp : miniTmpRC);
    return me.hValue;
}

// ----------------------------------------------------------------------------
// Function hash()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
hash(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    return _getMinimizer(me, it);
}

// ----------------------------------------------------------------------------
// Function hashNext()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
hashNext(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    return hash(me, it);
}


}	// namespace seqan

#endif
