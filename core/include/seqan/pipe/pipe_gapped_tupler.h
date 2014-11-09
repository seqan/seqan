// ==========================================================================
//                            pipe_gapped_tupler.h
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

#ifndef CORE_INCLUDE_SEQAN_PIPE_PIPE_GAPPED_TUPLER_H_
#define CORE_INCLUDE_SEQAN_PIPE_PIPE_GAPPED_TUPLER_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// GappedTupler
// --------------------------------------------------------------------------

template <typename TShape, bool omitLast = false, typename TPack = void>
struct GappedTupler;

// --------------------------------------------------------------------------
// Helper struct to fill the tuple from the text buffer
// --------------------------------------------------------------------------

template <typename TShape, typename TTuple>
struct _gappedTuplerFillTupleFromBuffer
{
    enum { BufferSize = TShape::span };
    enum { TupleSize = WEIGHT<TShape>::VALUE };

    template <typename TSize, typename TChar>
    inline void operator()(TSize const carePos[TupleSize], TChar const buffer[BufferSize], TTuple & tuple)
    {
        for(TSize i = 0; i < TupleSize; ++i)
            tuple[i] = buffer[carePos[i]];
    }
};

template <typename TShape, typename TValue>
struct _gappedTuplerFillTupleFromBuffer <TShape, Tuple<TValue, WEIGHT<TShape>::VALUE, BitPacked<> > >
{
    typedef Tuple<TValue, WEIGHT<TShape>::VALUE, BitPacked<> > TTuple;
    enum { BufferSize = TShape::span };
    enum { TupleSize = WEIGHT<TShape>::VALUE };

    template <typename TSize, typename TChar>
    inline void operator()(TSize const carePos[TupleSize], TChar const buffer[BufferSize], TTuple & tuple)
    {
        clear(tuple);
        for(TSize i = 0; i < TupleSize; ++i)
        {
            tuple <<= 1;
            tuple |= buffer[ carePos[i] ];
        }
    }
};


/*  Commented out because Ringbuffer version is slower than memmove verrsion :(


template <typename TSize, typename TShape>
void _createGappedTuplerMap(TSize map[TShape::span][WEIGHT<TShape>::VALUE], TShape const &)
{
    TSize weight = static_cast<TSize>(WEIGHT<TShape>::VALUE);
    TSize span   = static_cast<TSize>(TShape::span);

    // for index 0 fill in th usual carePositions
    for (TSize j=0; j < weight; ++j)
        map[0][j] = TShape::carePos[j];

    for (TSize i=1; i < span; ++i)
    {
        for (TSize j=0; j < weight; ++j)
        {
            TSize next = map[i-1][j] +1;
            if (next == span)
                map[i][j] = 0;
            else
                map[i][j] = next;
        }
    }
}

// --------------------------------------------------------------------------
// Pipe < TInput, GappedTupler >                                     [String]
// --------------------------------------------------------------------------

template <typename TInput, typename TShape, bool omitLast, typename TPack>
struct Pipe< TInput, GappedTupler<TShape, omitLast, TPack> >
{
    typedef typename Value<Pipe>::Type          TOutput;
    typedef typename Value<TOutput, 2 >::Type	TTuple;
    typedef typename Value<TTuple>::Type		TValue;
    typedef typename Size<TInput>::Type         TSize;

    enum { BufferSize = TShape::span };
    enum { TupleSize = LENGTH<TTuple>::VALUE };

    TInput                      &in;
    TOutput                     tmp;
    TSize                       lastTuples;

    TValue                      buffer[BufferSize];
    TSize                       buffIndex; // position in buffer to put next *in
    TSize                       map [BufferSize][TupleSize];

    Pipe(TInput& _in): in(_in), buffer()
    {
        _createGappedTuplerMap(map, TShape());
    }

    inline TOutput const & operator*() const {
        return tmp;
    }

    inline Pipe& operator++()
    {
        if (eof(in)) --lastTuples;

        if (lastTuples < TuplerNumberOfLastTuples_<BufferSize, omitLast>::VALUE)
        {
            buffer[buffIndex] = TValue();
            ++buffIndex;
            if (buffIndex >= BufferSize)
                buffIndex = 0;
        }
        else
        {
            buffer[buffIndex] = *in;
            ++in;
            ++buffIndex;
            if (buffIndex >= BufferSize)
                buffIndex = 0;
        }

        ++tmp.i1;
        _fillTmp2();
        return *this;
    }


    inline void fill() {

        for(buffIndex = 0; buffIndex < BufferSize && !eof(in); ++buffIndex, ++in)
            buffer[buffIndex] = *in;

        // set lastTuples depending on the omitLast flag
        if (TuplerNumberOfLastTuples_<BufferSize, omitLast>::VALUE > BufferSize - buffIndex)
            lastTuples = TuplerNumberOfLastTuples_<BufferSize, omitLast>::VALUE - (BufferSize - buffIndex);
        else
            lastTuples = 0; // this will cause eof() of this pipe

        // fill remaining buffer, if it hasn't been filled yet
        for (; buffIndex < BufferSize; ++buffIndex)
            buffer[buffIndex] = TValue();

        // fill tmp
        tmp.i1 = 0;
        buffIndex = 0;
        _fillTmp2();
    }

    inline void _fillTmp2()
    {
        for(TSize i = 0; i < TupleSize; ++i)
            tmp.i2[i] = buffer[ map[buffIndex][i] ];
    }
};


// --------------------------------------------------------------------------
// Pipe < TInput, GappedTupler >                                      [Multi]
// --------------------------------------------------------------------------

template <typename TInput, typename TShape, bool omitLast, typename TPack, typename TPair, typename TLimitsString >
struct Pipe< TInput, Multi<GappedTupler<TShape, omitLast, TPack>, TPair, TLimitsString> >
{
    typedef typename Value<Pipe>::Type                              TOutput;
    typedef typename Value< typename Value<Pipe>::Type, 2 >::Type	TTuple;
    typedef typename Value<TTuple>::Type							TValue;
    typedef typename Size<TInput>::Type                             TSize;
    typedef PairIncrementer_<TPair, TLimitsString>                  Incrementer;

    enum { BufferSize = TShape::span };
    enum { TupleSize = LENGTH<TTuple>::VALUE };

    TInput                     &in;
    Incrementer					localPos;
    TOutput                     tmp;
    TSize                       seqLength, lastTuples;
    TLimitsString const        &limits;

    TValue                      buffer[BufferSize];
    TSize                       buffIndex; // position in buffer to put next *in
    TSize                       map [BufferSize][TupleSize];


    template <typename TLimitsString_>
    // const &_limits is intentionally omitted to suppress implicit casts (if types mismatch) and taking refs of them
    Pipe(TInput& _in, TLimitsString_ &_limits):  in(_in), limits(_limits)
    {
        _createGappedTuplerMap(map, TShape());
    }

    inline TOutput const & operator*() const
    {
        return tmp;
    }

    inline Pipe& operator++()
    {
        // process next sequence
        if (eos())
            if (--lastTuples == 0)
            {
                assignValueI1(tmp.i1, getValueI1(static_cast<TPair>(localPos)));
                fill();
                return *this;
            }

        assignValueI2(tmp.i1, getValueI2(tmp.i1) + 1);

        if (lastTuples < TuplerNumberOfLastTuples_<BufferSize, omitLast>::VALUE)
        {
            buffer[buffIndex] = TValue();
            ++buffIndex;
            if (buffIndex >= BufferSize)
                buffIndex = 0;
        }
        else
        {
            buffer[buffIndex] = *in;
            ++localPos;
            ++in;
            ++buffIndex;
            if (buffIndex >= BufferSize)
                buffIndex = 0;
        }

        _fillTmp2();
        return *this;
    }

    inline void fill()
    {
        do {
            buffIndex = 0;
            if (!eof(in))
                do {
                    buffer[buffIndex] = *in;
                    ++in;
                    ++buffIndex;
                    ++localPos;
                } while ((buffIndex < BufferSize) && !eos());

            // lastTuples = 1 (omitLast true) or span (omitLast false)
            lastTuples = TuplerNumberOfLastTuples_<BufferSize, omitLast>::VALUE;

            // eventually, reduce the number of half-filled tuples
            if (lastTuples <= BufferSize - buffIndex)
                lastTuples = 0;
            else
            {
                lastTuples -= BufferSize - buffIndex;

                // fill up with null chars
                for(; buffIndex < BufferSize; ++buffIndex)
                    buffer[buffIndex] = TValue();
            }

            if (lastTuples == 0)
                assignValueI1(tmp.i1, getValueI1(static_cast<TPair>(localPos)));

        } while ((lastTuples == 0) && !eof(in));

        buffIndex = 0;
        assignValueI2(tmp.i1, 0);
        _fillTmp2();
    }

    inline bool eos() const
    {
        return (getValueI1(static_cast<TPair>(localPos)) > 0) && (getValueI2(static_cast<TPair>(localPos)) == 0);
    }

    inline void _fillTmp2()
    {
        // TODO: Use Loop struct?
        for(TSize i = 0; i < TupleSize; ++i)
            tmp.i2[i] = buffer[ map[buffIndex][i] ];
    }
};
*/

// --------------------------------------------------------------------------
// Pipe < TInput, GappedTupler >                                     [String]
// --------------------------------------------------------------------------

template <typename TInput, typename TShape, bool omitLast, typename TPack>
struct Pipe< TInput, GappedTupler<TShape, omitLast, TPack> >
{
    typedef typename Value<Pipe>::Type          TOutput;
    typedef typename Value<TOutput, 2 >::Type	TTuple;
    typedef typename Value<TTuple>::Type		TValue;
    typedef typename Size<TInput>::Type         TSize;

    enum { BufferSize = TShape::span };
    enum { TupleSize = LENGTH<TTuple>::VALUE };

    _gappedTuplerFillTupleFromBuffer<TShape,TTuple> _fillTmp2;

    TInput      &in;
    TOutput     tmp;
    TSize       lastTuples;

    // all elements will be shifted in ++ (ring buffer is complicated
    // due to many if(p > TShape::span) queries, maybe I will try that later)
    TValue      buffer[BufferSize];
    TSize       carePos[TupleSize];

    Pipe(TInput& _in): in(_in), buffer()
    {
        // TODO(meiers): These care positions of the shape are known at compile time
        //       They should be computed at compile time
        String<TSize> cpos;
        carePositions(cpos, TShape());
        for(TSize i=0; i< TupleSize; ++i)
            carePos[i] = cpos[i];
    }

    inline TOutput const & operator*() const {
        return tmp;
    }

    inline Pipe& operator++() {
        if (eof(in)) --lastTuples;

        // it's just a jump to the left
        // memmove is probably faster than unrolled loop:
        // Loop<ShiftLeftWorker2_, BufferSize - 1>::run(this->buffer);
        memmove(buffer, buffer+1, (BufferSize-1)*sizeof(TValue) );

        if (lastTuples < TuplerNumberOfLastTuples_<BufferSize, omitLast>::VALUE)
            buffer[BufferSize - 1] = TValue();
        else
        {
            buffer[BufferSize - 1] = *in;
            ++in;
        }

        ++tmp.i1;
        _fillTmp2(carePos, buffer, tmp.i2); // a bit expensive, but easier to implement
        return *this;
    }


    inline void fill() {

        unsigned i;
        for(i = 0; i < BufferSize && !eof(in); ++i, ++in)
            buffer[i] = *in;

        // set lastTuples depending on the omitLast flag
        if (TuplerNumberOfLastTuples_<BufferSize, omitLast>::VALUE > BufferSize - i)
            lastTuples = TuplerNumberOfLastTuples_<BufferSize, omitLast>::VALUE - (BufferSize - i);
        else
            lastTuples = 0; // this will cause eof() of this pipe

        // fill remaining buffer, if it hasn't been filled yet
        for (; i < BufferSize; ++i)
            buffer[i] = TValue();

        // fill tmp
        tmp.i1 = 0;
        _fillTmp2(carePos, buffer, tmp.i2);
    }

};


// --------------------------------------------------------------------------
// Pipe < TInput, GappedTupler >                                      [Multi]
// --------------------------------------------------------------------------

template <typename TInput, typename TShape, bool omitLast, typename TPack, typename TPair, typename TLimitsString >
struct Pipe< TInput, Multi<GappedTupler<TShape, omitLast, TPack>, TPair, TLimitsString> >
{
    typedef typename Value<Pipe>::Type                              TOutput;
    typedef typename Value< typename Value<Pipe>::Type, 2 >::Type	TTuple;
    typedef typename Value<TTuple>::Type							TValue;
    typedef typename Size<TInput>::Type                             TSize;
    typedef PairIncrementer_<TPair, TLimitsString>                  Incrementer;

    enum { BufferSize = TShape::span };
    enum { TupleSize = LENGTH<TTuple>::VALUE };

    _gappedTuplerFillTupleFromBuffer<TShape,TTuple> _fillTmp2;

    TInput                     &in;
    Incrementer					localPos;
    TOutput                     tmp;
    TSize                       seqLength, lastTuples;
    TLimitsString const        &limits;

    TValue                      buffer[BufferSize];
    TSize                       carePos[TupleSize];

    template <typename TLimitsString_>
    // const &_limits is intentionally omitted to suppress implicit casts (if types mismatch) and taking refs of them
    Pipe(TInput& _in, TLimitsString_ &_limits):  in(_in), limits(_limits)
    {
        /// TODO(meiers): These care positions of the shape are known at compile time
        //       They should be computed at compile time
        String<TSize> cpos;
        carePositions(cpos, TShape());
        for(TSize i=0; i< TupleSize; ++i)
            carePos[i] = cpos[i];
    }

    inline TOutput const & operator*() const
    {
        return tmp;
    }

    inline Pipe& operator++()
    {
        // process next sequence
        if (eos())
            if (--lastTuples == 0)
            {
                assignValueI1(tmp.i1, getValueI1(static_cast<TPair>(localPos)));
                fill();
                return *this;
            }

        // shift left 1 character
        memmove(buffer, buffer+1, (BufferSize-1)*sizeof(TValue) );
        assignValueI2(tmp.i1, getValueI2(tmp.i1) + 1);

        if (lastTuples < TuplerNumberOfLastTuples_<BufferSize, omitLast>::VALUE)
        {
            buffer[BufferSize - 1] = TValue();
        } else
        {
            buffer[BufferSize - 1] = *in;
            ++localPos;
            ++in;
        }

        _fillTmp2(carePos, buffer, tmp.i2);
        return *this;
    }

    inline void fill()
    {
        do {
            unsigned i = 0;
            if (!eof(in))
                do {
                    buffer[i] = *in;
                    ++in;
                    ++i;
                    ++localPos;
                } while ((i < BufferSize) && !eos());
            lastTuples = TuplerNumberOfLastTuples_<BufferSize, omitLast>::VALUE;

            // eventually, reduce the number of half-filled tuples
            if (lastTuples <= BufferSize - i)
                lastTuples = 0;
            else
            {
                lastTuples -= BufferSize - i;
                
                // fill up with null chars
                for(; i < BufferSize; ++i)
                    buffer[i] = TValue();
            }
            
            if (lastTuples == 0)
                assignValueI1(tmp.i1, getValueI1(tmp.i1) + 1);
            
        } while ((lastTuples == 0) && !eof(in));
        
        assignValueI2(tmp.i1, 0);
        _fillTmp2(carePos, buffer, tmp.i2);
    }
    
    inline bool eos() const
    {
        return (getValueI1(static_cast<TPair>(localPos)) > 0) && (getValueI2(static_cast<TPair>(localPos)) == 0);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction Value< Pipe < TInput, GappedTupler > >
// --------------------------------------------------------------------------

template <typename TInput, typename TShape, bool omitLast, typename TPack >
struct Value< Pipe< TInput, GappedTupler< TShape, omitLast, TPack > > >
{
    typedef Tuple<typename Value<TInput>::Type, WEIGHT<TShape>::VALUE, TPack>	TTuple;
    typedef Pair<typename Size<TInput>::Type, TTuple, Pack>         Type;
};

template <typename TInput, typename TShape, bool omitLast, typename TPack, typename TPair, typename TLimitsString >
struct Value< Pipe< TInput, Multi< GappedTupler< TShape, omitLast, TPack >, TPair, TLimitsString > > >
{
    typedef Tuple<typename Value<TInput>::Type, WEIGHT<TShape>::VALUE, TPack>	TTuple;
    typedef Pair<TPair, TTuple, Pack>                                           Type;
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Pipe control functions                                          [Sequence]
// --------------------------------------------------------------------------

template <typename TInput, typename TShape, bool omitLast, typename TPack >
inline bool
control(
        Pipe< TInput, GappedTupler< TShape, omitLast, TPack > > &me,
        ControlBeginRead const &command)
{
    if (!control(me.in, command))
        return false;
    me.fill();
    return true;
}

template <typename TInput, typename TShape, bool omitLast, typename TPack >
inline bool
control(
        Pipe< TInput, GappedTupler< TShape, omitLast, TPack > > &me,
        ControlEof const &)
{
    return me.lastTuples == 0;
}

template <typename TInput, typename TShape, bool omitLast, typename TPack >
inline bool
control(Pipe< TInput, GappedTupler< TShape, omitLast, TPack > > &me,
        ControlEos const &)
{
    return control(me, ControlEof());
}

// Note(meiers): work-around "unsigned expression >= 0 is always true"
template <typename TInput, typename TShape, bool omitLast, typename TPack >
inline typename Size< Pipe< TInput, GappedTupler<TShape, omitLast, TPack> > >::Type
_length(Pipe<TInput, GappedTupler< TShape, omitLast, TPack > > const &me, True const &)
{
    return length(me.in);
}
template <typename TInput, typename TShape, bool omitLast, typename TPack >
inline typename Size< Pipe< TInput, GappedTupler<TShape, omitLast, TPack> > >::Type
_length(Pipe<TInput, GappedTupler< TShape, omitLast, TPack > > const &me, False const &)
{
    typedef Pipe< TInput, GappedTupler< TShape, omitLast, TPack > >	TPipe;
    typedef TuplerNumberOfLastTuples_<TPipe::BufferSize, omitLast>  TLast;
    if (length(me.in) >= (TPipe::BufferSize - TLast::VALUE))
        return length(me.in) - (TPipe::BufferSize - TLast::VALUE);
    else
        return 0;
}
template <typename TInput, typename TShape, bool omitLast, typename TPack >
inline typename Size< Pipe< TInput, GappedTupler<TShape, omitLast, TPack> > >::Type
length(Pipe<TInput, GappedTupler< TShape, omitLast, TPack > > const &me)
{
    typedef Pipe< TInput, GappedTupler< TShape, omitLast, TPack > >	TPipe;
    typedef TuplerNumberOfLastTuples_<TPipe::BufferSize, omitLast>  TLast;
    typedef typename If<
        typename Eval<static_cast<unsigned>(TPipe::BufferSize) == static_cast<unsigned>(TLast::VALUE)>::Type, 
        True, 
        False>::Type TSwitch;
    return _length(me, TSwitch());
}

template <typename TInput, typename TShape, bool omitLast, typename TPack >
inline unsigned
countSequences(Pipe< TInput, GappedTupler< TShape, omitLast, TPack > > const &)
{
    return 1;
}

// --------------------------------------------------------------------------
// Pipe control functions                                             [Multi]
// --------------------------------------------------------------------------


template <typename TInput, typename TShape, bool omitLast, typename TPack, typename TPair, typename TLimitsString >
inline bool
control(Pipe< TInput, Multi<GappedTupler< TShape, omitLast, TPack >, TPair, TLimitsString> > &me,
        ControlBeginRead const &command)
{
    if (!control(me.in, command))
        return false;
    setHost(me.localPos, me.limits);
    assignValueI1(me.tmp.i1, 0);
    me.fill();
    return true;
}

template <typename TInput, typename TShape, bool omitLast, typename TPack, typename TPair, typename TLimitsString >
inline bool
control(Pipe< TInput, Multi<GappedTupler< TShape, omitLast, TPack >, TPair, TLimitsString> > &me,
        ControlEof const &)
{
    return (me.lastTuples == 0 && getValueI1(static_cast<TPair>(me.localPos)) >= length(me.limits) -1);
}

template <typename TInput, typename TShape, bool omitLast, typename TPack, typename TPair, typename TLimitsString >
inline bool
control(Pipe< TInput, Multi<GappedTupler< TShape, omitLast, TPack >, TPair, TLimitsString> > &me,
        ControlEos const &)
{
    return me.eos();
}

// Note(meiers): work-around "unsigned expression >= 0 is always true"
template <typename TInput, typename TShape, bool omitLast, typename TPack, typename TPair, typename TLimitsString >
inline typename Size< Pipe< TInput, Multi<GappedTupler< TShape, omitLast, TPack >, TPair, TLimitsString> > >::Type
_length(Pipe< TInput, Multi<GappedTupler< TShape, omitLast, TPack >, TPair, TLimitsString> > const &me, True const &)
{
    return length(me.in);
}
template <typename TInput, typename TShape, bool omitLast, typename TPack, typename TPair, typename TLimitsString >
inline typename Size< Pipe< TInput, Multi<GappedTupler< TShape, omitLast, TPack >, TPair, TLimitsString> > >::Type
_length(Pipe< TInput, Multi<GappedTupler< TShape, omitLast, TPack >, TPair, TLimitsString> > const &me, False const &)
{
    typedef Pipe< TInput, GappedTupler< TShape, omitLast, TPack > >	TPipe;
    typedef TuplerNumberOfLastTuples_<TPipe::BufferSize, omitLast>  TLast;
    unsigned seqs = countSequences(me);

    if (length(me.in) >= seqs * (TPipe::BufferSize - TLast::VALUE))
        return length(me.in) - seqs * (TPipe::BufferSize - TLast::VALUE);
    else
        return 0;
}
template <typename TInput, typename TShape, bool omitLast, typename TPack, typename TPair, typename TLimitsString >
inline typename Size< Pipe< TInput, Multi<GappedTupler< TShape, omitLast, TPack >, TPair, TLimitsString> > >::Type
length(Pipe< TInput, Multi<GappedTupler< TShape, omitLast, TPack >, TPair, TLimitsString> > const &me)
{
    typedef Pipe< TInput, GappedTupler< TShape, omitLast, TPack > >	TPipe;
    typedef TuplerNumberOfLastTuples_<TPipe::BufferSize, omitLast>  TLast;
    typedef typename If<
        typename Eval<static_cast<unsigned>(TPipe::BufferSize) == static_cast<unsigned>(TLast::VALUE)>::Type, 
        True, 
        False>::Type TSwitch;
    return _length(me, TSwitch());
}

template <typename TInput, typename TShape, bool omitLast, typename TPack, typename TPair, typename TLimitsString >
inline typename Size< Pipe< TInput, Multi<GappedTupler< TShape, omitLast, TPack >, TPair, TLimitsString> > >::Type
countSequences(Pipe< TInput, Multi<GappedTupler< TShape, omitLast, TPack >, TPair, TLimitsString> > const &me)
{
    return length(me.limits) - 1;
}

}

#endif  // #ifndef CORE_INCLUDE_SEQAN_PIPE_PIPE_GAPPED_TUPLER_H_
