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

const float _boundAlpha = 0.8;
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
    int first;
    int bound;
    THashValue hValue;      //minimizer hash
    THashValue mhValue;  
    THashValue last_hValue;
    THashValue u_hValue;
    static const THashValue leftFactor = Power<ValueSize<TValue>::VALUE, TSPAN - 1>::VALUE;
    static const THashValue m_leftFactor = Power<ValueSize<TValue>::VALUE, TWEIGHT - 1>::VALUE;
    TValue  m_leftChar; 
    TValue  leftChar;

    Shape():
        span(TSPAN),
        weight(TWEIGHT),
        first(-1),
        bound((unsigned)(TWEIGHT * _boundAlpha))
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

// ----------------------------------------------------------------------------
// Function _minHash()
// ----------------------------------------------------------------------------
// return lexicographically smaller hash as the minimizer

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TString>
inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type 
_minHash(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TString const & str)
{
    typedef typename Iterator<TString const, Standard>::Type                TIter;
    typedef typename Value<Shape<TValue, UngappedShape<TWEIGHT> > >::Type   THValue;
  
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    SEQAN_ASSERT_GT((unsigned)me.span, (unsigned)me.weight); 

    Shape<TValue, UngappedShape<TWEIGHT> > tmpShape;

    TIter strIt = begin(str, Standard());
    TIter strEnd = end(str, Standard()) - weight(me) + 1;

    THValue miniTmp = hash(tmpShape, strIt);
    THValue preMiniTmp; 
    me.first = 0;
    for (strIt++; strIt != strEnd; strIt++)
    {
        preMiniTmp = miniTmp; 
        miniTmp = _min(miniTmp, hashNext(tmpShape, strIt));
        if(preMiniTmp != miniTmp)
            me.first = strIt - begin(str, Standard()); 
    } 
    
    me.mhValue = miniTmp;

    if (me.first < me.bound) 
    {
        me.hValue = _max(me.mhValue, (Power<ValueSize<TValue>::VALUE, TWEIGHT>::VALUE - me.mhValue - 1));
    }
    else 
        me.hValue = me.mhValue;

    return me.hValue;
}


template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TString>
inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type 
_minXorMaxHash(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TString const & str)
{
    typedef typename Iterator<TString const, Standard>::Type                TIter;
    typedef typename Value<Shape<TValue, UngappedShape<TWEIGHT> > >::Type   THValue;
  
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    SEQAN_ASSERT_GT((unsigned)me.span, (unsigned)me.weight); 

    Shape<TValue, UngappedShape<TSPAN> > u_tmpShape;
    Shape<TValue, UngappedShape<TWEIGHT> > tmpShape;
   
    TIter strBegin = begin(str, Standard());

    THValue miniTmp;
    //THValue maxTmp = tmpShape.hValue;
    //THValue miniTmp2 = tmpShape.hValue;
    //THValue maxTmp2 = tmpShape.hValue;

    //me.u_hValue = hashInit(u_tmpShape, strBegin);    

    me.first = 0;
    hashInit(tmpShape, strBegin);
    miniTmp = hashNext(tmpShape, strBegin);
    for (unsigned k = 1; k < TSPAN - TWEIGHT + 1; k++)
    {
       // miniTmp2 = _min(miniTmp2, _max(miniTmp,hashNext(tmpShape, strIt)));
        //miniTmp = _min(miniTmp, tmpShape.hValue);
        //maxTmp2 = _max(maxTmp2, _min(maxTmp, tmpShape.hValue));
        //maxTmp = _max(maxTmp, tmpShape.hValue);
        hashNext(tmpShape, strBegin + k);
        //std::cout << miniTmp << " " << tmpShape.hValue << std::endl;
        if (miniTmp >= tmpShape.hValue)
        {
            miniTmp = tmpShape.hValue;
            me.first = k; 
        }
    } 
    //std::cout << "done" << std::endl;
    me.hValue = miniTmp;
    //me.hValue = atomicXor(maxTmp, miniTmp);
    //me.hValue = atomicXor(maxTmp2,  me.hValue);
    //me.hValue = atomicXor(miniTmp2, me.hValue);
    me.last_hValue = tmpShape.hValue; 
    return me.hValue;
}

/*
template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TString>
inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type 
_minXorMaxHash(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TString const & str)
{
    typedef typename Value<Shape<TValue, UngappedShape<TWEIGHT> > >::Type   THValue;
  
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    SEQAN_ASSERT_GT((unsigned)me.span, (unsigned)me.weight); 
    Shape<TValue, UngappedShape<TWEIGHT> > tmpShape;

    //TIter miniIt = it;
    //THValue maxTmp = tmpShape.hValue;
    //THValue miniTmp2 = tmpShape.hValue;
    //THValue maxTmp2 = tmpShape.hValue;
    me.first = 0;
    unsigned level = 0;
    unsigned chr[TSPAN];
    unsigned miniIt = 0;
  
    for (unsigned k = 0; k < TSPAN; k++)
        chr[k] = (unsigned)ordValue((TValue)str[k]);
    
    for (unsigned k = 1; k < TSPAN - TWEIGHT + 1 ; k++)
    {
        unsigned thisLevel = level;
        for (unsigned j = 0; j <= thisLevel; j++)
        {
            if (chr[k] < chr[miniIt + j])
            {
                miniIt = k - j;
                level = 0;
                break;
            } 
            if(chr[k] == chr[miniIt + j])
            {
                if (j == thisLevel)
                    if ( ++level == TWEIGHT)
                    {
                        miniIt = k - TWEIGHT + 1;
                        level = TWEIGHT - 1;
                    }
                continue;
            }
            else
                level = 0;
        }
       // miniTmp2 = _min(miniTmp2, _max(miniTmp,hashNext(tmpShape, strIt)));
        //miniTmp = _min(miniTmp, tmpShape.hValue);
        //maxTmp2 = _max(maxTmp2, _min(maxTmp, tmpShape.hValue));
        //maxTmp = _max(maxTmp, tmpShape.hValue);
        //std::cout << miniTmp << " " << tmpShape.hValue << std::endl;
    } 
    if (level != 0)
    for (unsigned k = TSPAN - TWEIGHT + 1; k <= TSPAN - 1; k++)
    {
        if (chr[k] < chr[miniIt + level])
        {
            miniIt = k - level;
            break;         
        }
        if (chr[k] == chr[miniIt + level])
        {
            if(++level == TWEIGHT)
            {
                miniIt = k - TWEIGHT + 1;
                break;
            }
        }
        else
            break; 
    }
    //std::cout << "done" << std::endl;
    //me.hValue = atomicXor(maxTmp, miniTmp);
    //me.hValue = atomicXor(maxTmp2,  me.hValue);
    //me.hValue = atomicXor(miniTmp2, me.hValue);
    me.first = miniIt; 
    me.last_hValue = hash(tmpShape, begin(str) + TSPAN - TWEIGHT);
    return hash(tmpShape, begin(str) + miniIt);
}
*/

/*
template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TString>
inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type 
_minXorMaxHash(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TString const & str)
{
    typedef typename Value<Shape<TValue, UngappedShape<TWEIGHT> > >::Type   THValue;
  
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    SEQAN_ASSERT_GT((unsigned)me.span, (unsigned)me.weight); 
    Shape<TValue, UngappedShape<TWEIGHT> > tmpShape;

    //TIter miniIt = it;
    //THValue maxTmp = tmpShape.hValue;
    //THValue miniTmp2 = tmpShape.hValue;
    //THValue maxTmp2 = tmpShape.hValue;
    me.first = 0;
    unsigned level = 0;
    unsigned chr[TSPAN];
    int loc[TSPAN];
    unsigned mini = 10;//ValueSize<TValue>::VALUE;
    unsigned size = TSPAN - TWEIGHT + 1;
            
    for (unsigned k = 0; k < TSPAN; k++)
    {
        chr[k] = (unsigned)ordValue((TValue)str[k]);
        loc[k] = k;
    }
    chr[TSPAN] = ValueSize<TValue>::VALUE;
  
    unsigned j;
    for (j = 0; j < TWEIGHT; j++)
    {
        unsigned this_size = size;
        for (unsigned k = 0; k < this_size; k++)
        {
            if(chr[loc[k] + j] < mini)
            {
                mini = chr[loc[k] + j];
                loc[0] = loc[k];
                size = 1;
                continue;
            }
            if(chr[loc[k] + j] == mini)
            {
                loc[size++] = loc[k];
            }
        }        
        if(size == 1)
            break;
        mini = 10;//ValueSize<TValue>::VALUE;
    } 
            
    me.first = loc[0]; 
    me.last_hValue = hash(tmpShape, begin(str) + TSPAN - TWEIGHT);
    //return hash(tmpShape, begin(str) + me.first);
    //std::cout << ValueSize<TValue>::VALUE << std::endl;
    
    return hash(tmpShape, begin(str) + me.first);
}
*/
// ----------------------------------------------------------------------------
// Function hash()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
hash(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    Range<TIter> range(it, it + length(me));
    Shape<TValue, UngappedShape<TSPAN> > u_tmpShape;

    me.hValue = _minXorMaxHash(me, range);
    me.m_leftChar = *(it + TSPAN -TWEIGHT);
    return me.hValue;
}

// ----------------------------------------------------------------------------
// Function hash(); ReverseComplement
// ----------------------------------------------------------------------------
// Uses the lexicographically smaller of S and the reverse complement of S as the minimizer

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TIter>
inline typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, ReverseComplement> > >::Type
hash(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, ReverseComplement> > &me, TIter const &it)
{
    typedef Range<TIter>                        TContainer;
    typedef typename Value<TContainer>::Type    TAlphabet;
    typedef typename RemoveConst<TAlphabet>::Type TNCAlphabet;
    typedef ModifiedString<TContainer, ModView<FunctorComplement<TNCAlphabet> > > TComplement;
    typedef ModifiedString<TComplement, ModReverse>                             TRC;

    Range<TIter> range(it, it + length(me));
    TRC revComplRange(range);

    me.hValue = _min(_minHash(me, range), _minHash(me, revComplRange));
    return me.hValue;
}

// ----------------------------------------------------------------------------
// Function hashNext()
// ----------------------------------------------------------------------------

/*template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
hashNext(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    Shape<TValue, UngappedShape<TWEIGHT> > tmpShape;
    if(me.first > 0) 
    {
        hash(tmpShape, it + me.span - me.weight);
        if(me.mhValue < tmpShape.hValue)
            me.first--;
        else 
        {
            me.mhValue = tmpShape.hValue;
            me.first = TSPAN - TWEIGHT;
        }

        if (me.first < me.bound) 
        {
            me.hValue = _max(me.mhValue, (Power<ValueSize<TValue>::VALUE, TWEIGHT>::VALUE - me.mhValue - 1));
        }
        else 
            me.hValue = me.mhValue;
        me.hValue = hash(me, it);
    return me.hValue;
}
*/

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline void hashInit(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    hash(me, it);
    Shape<TValue, UngappedShape<TSPAN> > tmpShape;
    Shape<TValue, UngappedShape<TWEIGHT> > u_shape;
    me.last_hValue = hash(u_shape, it + TSPAN - TWEIGHT -1);
    me.m_leftChar = *(it + TSPAN - TWEIGHT - 1);
    me.u_hValue = hashInit(tmpShape, it);
    me.leftChar = 0;
}


template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
hashNext(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    Shape<TValue, UngappedShape<TWEIGHT> > tmpShape;
    Shape<TValue, UngappedShape<TSPAN> > u_tmpShape;  
    typedef typename Value<Shape<TValue, UngappedShape<TWEIGHT> > >::Type THValue;
    typedef typename Size<Shape<TValue, UngappedShape<TWEIGHT> > >::Type TSize;
    if(me.first > 0) 
    {
        me.last_hValue =
        (me.last_hValue - ordValue(me.m_leftChar) * (THValue)me.m_leftFactor) * ValueSize<TValue>::VALUE
        + ordValue((TValue)*(it + ((TSize)TSPAN - 1)));
        me.m_leftChar = *(it + TSPAN - TWEIGHT);
        if (me.hValue < me.last_hValue)
            me.first--; 
        else 
        {
            me.hValue = me.last_hValue;
            me.first = TSPAN - TWEIGHT;
        }

    }
    else 
        me.hValue = hash(me, it);

    //std::cout << "me.u_hValue = " << me.u_hValue << "me.leftChar = " << me.leftChar << std::endl;
    
    me.u_hValue = 
            (me.u_hValue - ordValue(me.leftChar) * (THValue)me.leftFactor) * ValueSize<TValue>::VALUE
            + ordValue((TValue)*(it + ((TSize)TSPAN - 1)));
    me.leftChar = *it;
    return me.hValue;
}

}	// namespace seqan


#endif
