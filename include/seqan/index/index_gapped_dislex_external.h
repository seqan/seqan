// ==========================================================================
//                       index_gapped_dislex_external.h
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
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_DISLEX_EXTERNAL_H_
#define CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_DISLEX_EXTERNAL_H_

//#define DISLEX_EXTERNAL_RUNNING_TIMES

namespace SEQAN_NAMESPACE_MAIN
{
    template <typename TShape, typename TSACA = Skew7>
    struct DislexExternal {};

// --------------------------------------------------------------------------
// Filter to transform positions to lengths                          [String]
// --------------------------------------------------------------------------

template <typename TSize>
struct PositionToLengthTransform_
{
    TSize N;
    PositionToLengthTransform_(TSize len): N(len)
    {}

    inline TSize operator()(TSize pos) const
    {
        return N - pos;
    }
};

// --------------------------------------------------------------------------
// Filter to transform positions to lengths                       [StringSet]
// --------------------------------------------------------------------------

template <typename TLimitString, typename TIdAndPosPair>
struct PositionToLengthTransformMulti_
{
    typedef typename Value<TIdAndPosPair, 2>::Type TSize;
    TLimitString const & _limits;

    PositionToLengthTransformMulti_(TLimitString const & limitStr) : _limits(limitStr)
    {}

    inline TIdAndPosPair operator() (TIdAndPosPair x) const
    {
        TSize N = _limits[x.i1 + 1] - _limits[x.i1];
        x.i2 = N - x.i2;
        return x;
    }
};

// --------------------------------------------------------------------------
// Comparator for naming tuples                                      [String]
// --------------------------------------------------------------------------

/*
 * @signature DislexTupleComp_<TValue, TShape, TResult = int>
 *
 * @tparam TValue expects a Pair<TSize, Tuple> where the 1st parameter contains the
 *                <b>length</b> of the underlying suffix and the 2nd parameter the
 *                fixed-length sequence tuple (possibly bitpacked)
 * @tparam TShape expects a fixed CyclicShape (CyclicShape<FixedShape<...> >)
 *
 * Only for hardwired cyclic shapes! There is a overload for bitpacked tuples
 *
 * @see DislexTupleCompMulti_
 */
template <typename TValue, typename TShape, typename TResult=int>
struct DislexTupleComp_ :
    public std::binary_function<TValue, TValue, TResult>
{
    typedef typename Value<TValue, 1>::Type                 TSize;
    typedef typename Value<TValue, 2>::Type                 TTuple;
    typedef typename Value<TTuple>::Type                    TTupleValue;
    typedef typename StoredTupleValue_<TTupleValue>::Type TStoredValue;

    enum
    {
        _span = TShape::span,
        _weight = WEIGHT<TShape>::VALUE
    };
    TSize realLengths[_span];
    PositionToLengthTransform_<TSize> posToLen;

    DislexTupleComp_(TSize strLen) : posToLen(strLen)
    {
        #ifdef DISLEX_EXTERNAL_RUNNING_TIMES
        std::cout << "   | String, Tuple version"  << std::endl;
        #endif
        cyclicShapeToSuffixLengths(realLengths, TShape());
    }

    inline TResult operator() (const TValue &a, const TValue &b) const
    {
        const TStoredValue * sa = a.i2.i;
        const TStoredValue * sb = b.i2.i;

        TSize la = posToLen(a.i1);
        TSize lb = posToLen(b.i1);

        // find out the real lengths of the gapped strings
        TSize rla = (la < static_cast<TSize>(_span) ? realLengths[la] : static_cast<TSize>(_weight));
        TSize rlb = (lb < static_cast<TSize>(_span) ? realLengths[lb] : static_cast<TSize>(_weight));

        // lexicographical comparison
        TSize n = std::min(static_cast<TSize>(_weight), std::min(rla, rlb) );
        for (TSize i = 0; i < n; i++, ++sa, ++sb)
        {
            if (*sa == *sb)
                continue;
            return (*sa < *sb)? -1 : 1;
        }

        // both cyclic shapes are "full"
        if (la >= static_cast<TSize>(_span) && lb >= static_cast<TSize>(_span))
            return 0;

        // if they are NOT equally long
        if (rla != rlb)
            return (rla < rlb ? -1 : 1);

        // if they are equally long
        if (la > lb) return 1;
        if (la < lb) return -1;

        // only occurs when q-grams from the exact same position are passed
        SEQAN_ASSERT_EQ(a.i1, b.i1);
        return 0;
    }
};

// BitPacked version of DislexTupleComp_

template <typename TSize,
          typename TTupleValue,
          typename TShape,
          typename TResult>
struct DislexTupleComp_<Pair<TSize, Tuple<TTupleValue, WEIGHT<TShape>::VALUE, BitPacked<> >, Pack>,
                        TShape,
                        TResult> :
    public std::binary_function<Pair<TSize, Tuple<TTupleValue, WEIGHT<TShape>::VALUE, BitPacked<> >, Pack>,
                                Pair<TSize, Tuple<TTupleValue, WEIGHT<TShape>::VALUE, BitPacked<> >, Pack>,
                                TResult>
{
    typedef Tuple<TTupleValue, WEIGHT<TShape>::VALUE, BitPacked<> > TTuple;
    typedef Pair<TSize, TTuple, Pack>                               TValue;

    enum
    {
        _span = TShape::span,
        _weight = WEIGHT<TShape>::VALUE
    };
    TSize realLengths[_span];
    PositionToLengthTransform_<TSize> posToLen;

    DislexTupleComp_(TSize strLen) : posToLen(strLen)
    {
        #ifdef DISLEX_EXTERNAL_RUNNING_TIMES
        std::cout << "   | String, Bitpacked version"  << std::endl;
        #endif
        cyclicShapeToSuffixLengths(realLengths, TShape());
    }

    inline TResult operator() (const TValue &a, const TValue &b) const
    {
        // compare Tuples right away (filled with 0s in the rear)
        if (a.i2 < b.i2) return -1;
        if (a.i2 > b.i2) return 1;

        TSize la = posToLen(a.i1);
        TSize lb = posToLen(b.i1);

        // find out the real lengths of the gapped strings
        TSize rla = (la < static_cast<TSize>(_span) ? realLengths[la] : static_cast<TSize>(_weight));
        TSize rlb = (lb < static_cast<TSize>(_span) ? realLengths[lb] : static_cast<TSize>(_weight));

        // both cyclic shapes are "full"
        if (la >= static_cast<TSize>(_span) && lb >= static_cast<TSize>(_span))
            return 0;

        // if they are NOT equally long
        if (rla != rlb)
            return (rla < rlb ? -1 : 1);

        // if they are equally long
        if (la > lb) return 1;
        if (la < lb) return -1;

        // only occurs when q-grams from the exact same position are passed
        SEQAN_ASSERT_EQ(a.i1, b.i1);
        return 0;
    }
};


// --------------------------------------------------------------------------
// Comparator for naming tuples                                   [StringSet]
// --------------------------------------------------------------------------

/*
 * @signature DislexTupleCompMulti_<TValue, TShape, TResult = int>
 *
 * @tparam TValue expects a Pair<Pair<TSize, TSize>, Tuple> where the 1st parameter
 *                is a Pair of sequence ID and suffix <b>length</b> and the 2nd parameter
 *                the fixed-length sequence tuple (possibly bitpacked)
 * @tparam TShape expects a fixed CyclicShape (CyclicShape<FixedShape<...> >)
 *
 * @see DislexTupleComp_
 */
template <typename TValue,
          typename TShape,
          typename TLimitString,
          typename TResult=int>
struct DislexTupleCompMulti_  :
    public std::binary_function<TValue, TValue, TResult>
{
    typedef typename Value<TValue, 1>::Type                 TSetPos;
    typedef typename Value<TSetPos, 2>::Type                TSize;
    typedef typename Value<TValue, 2>::Type                 TTuple;
    typedef typename Value<TTuple>::Type                    TTupleValue;
    typedef typename StoredTupleValue_<TTupleValue>::Type TStoredValue;

    enum {
        _span = TShape::span,
        _weight = WEIGHT<TShape>::VALUE
    };
    TSize realLengths[_span];
    PositionToLengthTransformMulti_<TLimitString, TSetPos> posToLen;

    DislexTupleCompMulti_(TLimitString const & strSetLimits) : posToLen(strSetLimits)
    {
        #ifdef DISLEX_EXTERNAL_RUNNING_TIMES
        std::cout << "   | String Set, Tuple version"  << std::endl;
        #endif
        cyclicShapeToSuffixLengths(realLengths, TShape());
    }

    inline TResult operator() (const TValue &a, const TValue &b) const
    {
        const TStoredValue * sa = a.i2.i;
        const TStoredValue * sb = b.i2.i;

        TSetPos la = posToLen(a.i1);
        TSetPos lb = posToLen(b.i1);

        // find out the real lengths of the gapped strings
        TSize rla = (la.i2 < static_cast<TSize>(_span) ? realLengths[la.i2] : static_cast<TSize>(_weight));
        TSize rlb = (lb.i2 < static_cast<TSize>(_span) ? realLengths[lb.i2] : static_cast<TSize>(_weight));
        
        // lexicographical comparison
        TSize n = std::min(static_cast<TSize>(_weight), std::min(rla, rlb) );
        for (TSize i = 0; i < n; i++, ++sa, ++sb)
        {
            if (*sa == *sb)
                continue;
            return (*sa < *sb)? -1 : 1;
        }

        // both cyclic shapes are more than "full"
        if (la.i2 > static_cast<TSize>(_span) && lb.i2 > static_cast<TSize>(_span))
            return 0;

        // if they are NOT equally long
        if (rla != rlb)
            return (rla < rlb ? -1 : 1);

        // if they are equally long
        if (la.i2 == lb.i2)
        {
            if (la.i1 < lb.i1) return 1;
            if (la.i1 > lb.i1) return -1;
        } else {
            if (la.i2 > lb.i2) return 1;
            if (la.i2 < lb.i2) return -1;
        }

        // only occurs when q-grams from the exact same position are passed
        SEQAN_ASSERT_EQ(a.i1, b.i1);
        return 0;
    }
};

// BitPacked version of DislexTupleCompMulti_

template <typename TSetPos,
          typename TTupleValue,
          typename TShape,
          typename TLimitString,
          typename TResult>
struct DislexTupleCompMulti_<Pair<TSetPos, Tuple<TTupleValue, WEIGHT<TShape>::VALUE, BitPacked<> >, Pack>,
                             TShape,
                             TLimitString,
                             TResult> :
    public std::binary_function<Pair<TSetPos, Tuple<TTupleValue, WEIGHT<TShape>::VALUE, BitPacked<> >, Pack>,
                                Pair<TSetPos, Tuple<TTupleValue, WEIGHT<TShape>::VALUE, BitPacked<> >, Pack>,
                                TResult>
{
    typedef Tuple<TTupleValue, WEIGHT<TShape>::VALUE, BitPacked<> > TTuple;
    typedef Pair<TSetPos, TTuple, Pack>                             TValue;
    typedef typename Value<TSetPos, 2>::Type                        TSize;

    enum {
        _span = TShape::span,
        _weight = WEIGHT<TShape>::VALUE
    };
    TSize realLengths[_span];
    PositionToLengthTransformMulti_<TLimitString, TSetPos> posToLen;

    DislexTupleCompMulti_(TLimitString const & strSetLimits) : posToLen(strSetLimits)
    {
        #ifdef DISLEX_EXTERNAL_RUNNING_TIMES
        std::cout << "   | String Set, BitPacked version"  << std::endl;
        #endif
        cyclicShapeToSuffixLengths(realLengths, TShape());
    }


    inline TResult operator() (const TValue &a, const TValue &b) const
    {
        // compare Tuples right away (filled with 0s in the rear)
        if (a.i2 < b.i2) return -1;
        if (a.i2 > b.i2) return 1;

        TSetPos la = posToLen(a.i1);
        TSetPos lb = posToLen(b.i1);

        // find out the real lengths of the gapped strings
        TSize rla = (la.i2 < static_cast<TSize>(_span) ? realLengths[la.i2] : static_cast<TSize>(_weight));
        TSize rlb = (lb.i2 < static_cast<TSize>(_span) ? realLengths[lb.i2] : static_cast<TSize>(_weight));
        
        // both cyclic shapes are more than "full"
        if (la.i2 > static_cast<TSize>(_span) && lb.i2 > static_cast<TSize>(_span))
            return 0;

        // if they are NOT equally long
        if (rla != rlb)
            return (rla < rlb ? -1 : 1);

        // if they are equally long
        if (la.i2 == lb.i2)
        {
            if (la.i1 < lb.i1) return 1;
            if (la.i1 > lb.i1) return -1;
        } else {
            if (la.i2 > lb.i2) return 1;
            if (la.i2 < lb.i2) return -1;
        }

        // only occurs when q-grams from the exact same position are passed
        SEQAN_ASSERT_EQ(a.i1, b.i1);
        return 0;
    }
};


// --------------------------------------------------------------------------
// Mapping functor from text to lexText                              [String]
// --------------------------------------------------------------------------

// wrapper for the dislex Pipe
// takes a tuple <l, ACGACA> where p is the suffix position
// and returns L(N-l)
template <typename TValue,
          typename TResult = typename Value<TValue, 1>::Type>
struct DislexMap_ :
    public std::unary_function<TValue, TResult>
{
    DislexTransform_<TResult> formula;

    DislexMap_(TResult S_, TResult N_) : formula(S_, N_)
    {}

    inline TResult operator() (const TValue & x) const
    {
        return formula(x.i1);
    }
};


// --------------------------------------------------------------------------
// Mapping functor from text to lexText                           [StringSet]
// --------------------------------------------------------------------------

// dislex transformation used in the mapper pool
// takes a Pair( Pair(s,p), ACGATCG), where s is the seq id and p the suffix position,
// returns a global position L(s,p)=pos
template <typename TValue,
          typename Tlimits,
          typename TResult = typename Value<typename Value<TValue, 1>::Type, 2>::Type >
struct DislexMapMulti_ :
    public std::unary_function<TValue, TResult>
{
    typedef typename Value<TValue, 1>::Type TPair;
    typedef typename Value<TPair, 2>::Type TSize;

    DislexTransformMulti_<TPair, Tlimits> formula;
    
    DislexMapMulti_(TResult S_, Tlimits const & stringSetLimits) : formula(S_, stringSetLimits)
    {}
    
    inline TResult operator() (const TValue & x) const
    {
        return formula(x.i1);
    }
};

// --------------------------------------------------------------------------
// Pipe Dislex                                                       [String]
// --------------------------------------------------------------------------

    // TODO(meiers): Define metafunctions, e.g. Value, or do I need them?

template <typename TInput, typename TShape, typename TSACA>
struct Pipe<TInput, DislexExternal<TShape, TSACA> >
{
    typedef typename If<typename Eval< BitsPerValue<TypeOf_(TInput)>::VALUE * WEIGHT<TShape>::VALUE <= 64  >::Type,
                        BitPacked<>, Pack >::Type               TPack;

    typedef Pipe<TInput, GappedTupler<TShape, false, TPack> >   TPipeTupler;
    typedef DislexTupleComp_<TypeOf_(TPipeTupler), TShape>      TTupleComparator;
    typedef Pool<TypeOf_(TPipeTupler), SorterSpec<
            SorterConfigSize<TTupleComparator,
            TSizeOf_(TPipeTupler) > > >                         TPoolSorter;

    typedef Pipe< TPoolSorter, Namer<TTupleComparator> >        TPipeNamer;
    typedef DislexMap_<TypeOf_(TPipeNamer) >                    TDislexMapper;
    typedef Pool< TypeOf_(TPipeNamer), MapperSpec<
            MapperConfigSize< TDislexMapper,
            TSizeOf_(TPipeNamer) > > >                          TPoolMapper;

    typedef Pipe< TPoolMapper, Filter<
            filterI2<TypeOf_(TPoolMapper)> > >                  TPipeFilterI2;
    typedef Pipe<TPipeFilterI2, TSACA>                          TPipeSACA;
    typedef _dislexReverseTransform<TypeOf_(TPipeSACA),
            TypeOf_(Pipe)>                                      TDislexReverse;
    typedef Pipe<TPipeSACA, Filter<TDislexReverse> >            TPipeReverseTransform;


    TPipeSACA pool;                 // last pool (skew); will be filled when calling process().
    TPipeReverseTransform in;       // final Pipe

    Pipe()
    {}

    Pipe(TInput & textIn) : in(pool, TDislexReverse(TShape::span, length(textIn)))
    {
        // fill pool right away
        process(textIn);
    }

    inline typename Value<Pipe>::Type const operator*() {
        return *in;
    }

    inline Pipe& operator++() {
        ++in;
        return *this;
    }

    template < typename TInput_ >
    bool process(TInput_ &textIn)
    {
        // 1. Generate Gapped Tuples
        TPipeTupler                                             tupler(textIn);

        // 2. Sort Tuples by the first few characters
        TTupleComparator                                        _comp(length(textIn));
        TPoolSorter                                             sorter(tupler, _comp);
        
        #ifdef DISLEX_EXTERNAL_RUNNING_TIMES
        double teim = sysTime();
        #endif

        sorter << tupler;

        #ifdef DISLEX_EXTERNAL_RUNNING_TIMES
        std::cout << "   | sorter << tupler: " << sysTime() - teim << "s" << std::endl; teim = sysTime();
        #endif

        // 3. Name tuples by their rank
        TPipeNamer                                              namer(sorter, _comp);

        // 4. Map text Positions to lexText positions
        TDislexMapper                                           _map(TShape::span, length(textIn));
        TPoolMapper                                             mapper(namer, _map);
                
        mapper << namer;

        #ifdef DISLEX_EXTERNAL_RUNNING_TIMES
        std::cout << "   | mapper << namer:  " << sysTime() - teim << "s\tsigma = " << (namer.tmp.i2 +1)<< std::endl; teim = sysTime();
        #endif

        // 5. Discard positions, keep rank
        TPipeFilterI2                                           filter(mapper);

        // 6. Run SACA on lex text
        pool << filter;

        #ifdef DISLEX_EXTERNAL_RUNNING_TIMES
        std::cout << "   | pool << filter:   " << sysTime() - teim << "s (len = " << length(textIn) << ")" << std::endl; teim = sysTime();
        #endif

        // 7. Reverse Transform is done during the reading process
        return true;
    }
};

// --------------------------------------------------------------------------
// Pipe Dislex                                                    [StringSet]
// --------------------------------------------------------------------------

template <typename TInput, typename TShape, typename TSACA, typename TPair, typename TLimits>
struct Pipe<TInput, Multi<DislexExternal<TShape, TSACA>, TPair, TLimits> >
{
    typedef typename If<typename Eval< BitsPerValue<TypeOf_(TInput)>::VALUE * WEIGHT<TShape>::VALUE <= 64  >::Type,
                        BitPacked<>, Pack >::Type               TPack;

    typedef Pipe<TInput, Multi<GappedTupler<TShape, false, TPack>,
            TPair, TLimits> >                                   TPipeTupler;
    typedef DislexTupleCompMulti_<TypeOf_(TPipeTupler),
            TShape, TLimits>                                    TTupleComparator;
    typedef Pool<TypeOf_(TPipeTupler), SorterSpec<
            SorterConfigSize<TTupleComparator,
            TSizeOf_(TPipeTupler) > > >                         TPoolSorter;

    typedef Pipe< TPoolSorter, Namer<TTupleComparator> >        TPipeNamer;
    typedef DislexMapMulti_<TypeOf_(TPipeNamer), TLimits>       TDislexMapper;
    typedef Pool< TypeOf_(TPipeNamer), MapperSpec<
            MapperConfigSize< TDislexMapper,
            TSizeOf_(TPipeNamer) > > >                          TPoolMapper;

    typedef Pipe< TPoolMapper, Filter<
            filterI2<TypeOf_(TPoolMapper)> > >                  TPipeFilterI2;
    typedef Pipe<TPipeFilterI2, TSACA>                          TPipeSACA;
    typedef _dislexReverseTransformMulti<TypeOf_(TPipeSACA),
            TLimits, TypeOf_(Pipe)>                             TDislexReverse;
    typedef Pipe<TPipeSACA, Filter<TDislexReverse> >            TPipeReverseTransform;


    TLimits const & _limits;         // StringSetLimits
    TPipeSACA pool;                 // last pool (skew); will be filled when calling process().
    TPipeReverseTransform in;       // final Pipe


    template <typename TLimits_>
    Pipe(TLimits_ const & strSetLimits, SEQAN_CTOR_ENABLE_IF(IsSameType<TLimits, TLimits_>)) :
        _limits(strSetLimits)
    {
        ignoreUnusedVariableWarning(dummy);
    }

    template <typename TLimits_>
    Pipe(TInput& _textIn, TLimits_ const & strSetLimits, SEQAN_CTOR_ENABLE_IF(IsSameType<TLimits, TLimits_>)) :
        _limits(strSetLimits), in(pool, TDislexReverse(TShape::span, strSetLimits))
    {
        // fill pool right away
        process(_textIn);
        ignoreUnusedVariableWarning(dummy);
    }

    inline typename Value<Pipe>::Type const operator*()
    {
        return *in;
    }

    inline Pipe& operator++()
    {
        ++in;
        return *this;
    }

    template < typename TInput_ >
    bool process(TInput_ &textIn)
    {
        // 1. Generate Gapped Tuples
        TPipeTupler                                             tupler(textIn, _limits);

        // 2. Sort Tuples by the first few characters
        TTupleComparator                                        _comp(_limits);
        TPoolSorter                                             sorter(tupler, _comp);

        #ifdef DISLEX_EXTERNAL_RUNNING_TIMES
        double teim = sysTime();
        #endif

        sorter << tupler;

        #ifdef DISLEX_EXTERNAL_RUNNING_TIMES
        std::cout << "   | sorter << tupler: " << sysTime() - teim << "s" << std::endl; teim = sysTime();
        #endif

        // 3. Name tuples by their rank
        TPipeNamer                                              namer(sorter, _comp);

        // 4. Map text Positions to lexText positions
        TDislexMapper                                           _map(TShape::span, _limits);
        TPoolMapper                                             mapper(namer, _map);
        mapper << namer;

        #ifdef DISLEX_EXTERNAL_RUNNING_TIMES
        std::cout << "   | mapper << namer:  " << sysTime() - teim << "s\tsigma = " << (namer.tmp.i2 +1) << std::endl; teim = sysTime();
        #endif

        // 5. Discard positions, keep rank
        TPipeFilterI2                                           filter(mapper);

        // 6. Run SACA on lex text
        pool << filter;

        #ifdef DISLEX_EXTERNAL_RUNNING_TIMES
        std::cout << "   | pool << filter:   " << sysTime() - teim << "s (len = " << length(textIn) << ")" << std::endl; teim = sysTime();
        #endif
        
        // 7. Reverse Transform is done during the reading process
        return true;
    }
};


// ==========================================================================
// Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Value
// --------------------------------------------------------------------------

template < typename TInput, typename TA, typename TB >
struct Value< Pipe< TInput, DislexExternal<TA,TB> > > :
    public Size<TInput>
{};

template < typename TInput, typename TA, typename TB, typename TPair, typename TLimits>
struct Value< Pipe< TInput, Multi<DislexExternal<TA, TB>, TPair, TLimits> > >
{
    typedef TPair Type;
};


// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Operator << for Pipes                                             [String]
// --------------------------------------------------------------------------

template < typename TInput, typename TShape, typename TSACA, typename TObject >
inline bool operator<<(Pipe< TInput, DislexExternal<TShape, TSACA> > &me, TObject &textIn)
{
    typedef Pipe< TInput, DislexExternal<TShape, TSACA> > TPipe;
    me.in = TPipe::TPipeReverseTransform(me.pool, TPipe::TDislexReverse(TShape::span, length(textIn)));
    return me.process(textIn);
}

// --------------------------------------------------------------------------
// Operator << for Pipes                                             [String]
// --------------------------------------------------------------------------

template < typename TInput, typename TShape, typename TSACA, typename TPair, typename TLimits, typename TObject >
inline bool operator<<(
    Pipe< TInput, Multi<DislexExternal<TShape, TSACA>, TPair, TLimits> > &me,
    TObject &textIn)
{
    typedef Pipe< TInput, Multi<DislexExternal<TShape, TSACA>, TPair, TLimits> > TPipe;
    me._limits = stringSetLimits(textIn);
    me.in = TPipe::TPipeReverseTransform(me.pool, TPipe::TDislexReverse(TShape::span, me._limits));
    return me.process(textIn);
}


// --------------------------------------------------------------------------
// function createGappedSuffixArray()                                [String]
// --------------------------------------------------------------------------

template < typename TSA, typename TText, typename TCyclicShape, typename TSACA>
inline void createGappedSuffixArray(
    TSA &SA, // must already be resized already
    TText const &s,
    TCyclicShape const &,
    ModCyclicShape<TCyclicShape> const &,
    DislexExternal<TCyclicShape, TSACA> const &)
{
    _createSuffixArrayPipelining(SA, s, DislexExternal<TCyclicShape, TSACA>());
}

// --------------------------------------------------------------------------
// function createGappedSuffixArray()                             [StringSet]
// --------------------------------------------------------------------------

template < typename TSA, typename TText, typename TSpec, typename TCyclicShape, typename TSACA>
inline void _createGappedSuffixArrayPipelining(
    TSA &SA, // must already be resized already
    StringSet<TText, TSpec> const &s,
    TCyclicShape const &,
    ModCyclicShape<TCyclicShape> const &,
    DislexExternal<TCyclicShape, TSACA> const &)
{
    _createSuffixArrayPipelining(SA, s, DislexExternal<TCyclicShape, TSACA>());
}


}
#endif  // #ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_DISLEX_EXTERNAL_H_
