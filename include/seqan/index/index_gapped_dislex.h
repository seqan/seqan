// ==========================================================================
//                          index_gapped_sa_dislex.h
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

#ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_SA_DISLEX_H_
#define CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_SA_DISLEX_H_

//#define DISLEX_INTERNAL_RUNNING_TIMES

namespace SEQAN_NAMESPACE_MAIN
{

// ==========================================================================
// DisLex mapping functors
// ==========================================================================

// these functors implement the logic of mapping text positions to lexText
// positions and vice versa. They are used by the internal as well as the
// external Dislex algorithm.

// TODO(meiers): Write some documentation(?)

    
// --------------------------------------------------------------------------
// Map: Text -> LexText                                              [String]
// --------------------------------------------------------------------------

// real dislex transformation:
// takes a position x and returns L(x)
template <typename TInput, typename TResult = TInput>
struct DislexTransform_ :
    public std::unary_function<TInput, TResult>
{
    TInput const _s, _n, _b, _c;

    DislexTransform_(TInput s, TInput n) : _s(s), _n(n), _b(n / s), _c(n % s)
    {}

    inline TResult operator() (TInput const & x) const
    {
        TInput r = x / _s;
        TInput bb = (_n - x) % _s;
        TInput ret = (_s - 1 - bb) * _b + r;
        if (bb <= _c)
            ret += _c - bb;
        return static_cast<TResult>(ret);
    }
};


// --------------------------------------------------------------------------
// Map: Text -> LexText                                           [StringSet]
// --------------------------------------------------------------------------

// dislex transformation for multiple strings
// takes a Pair (seq,pos) and returns L(seq,pos) = pos
template <typename TValue,
          typename TString, typename TResult = typename Value<TValue, 2>::Type>
struct DislexTransformMulti_ :
    public std::unary_function<TValue, TResult>
{
    TString const & _limits;
    TResult const _s;

    DislexTransformMulti_(TResult s, TString const & strSetLimits) : _limits(strSetLimits), _s(s)
    {}

    inline TResult operator() (const TValue & x) const
    {
        TResult seq = x.i1;
        TResult pos = x.i2;
        TResult n = _limits[seq+1] - _limits[seq];

        // TODO: Store the following values for each String?
        TResult b = n / _s;
        TResult c = n % _s;

        TResult r = pos / _s;
        TResult bb = (n - pos) % _s;
        TResult ret = _limits[seq] + (_s - 1 - bb) * b + r;
        if (bb > c)
            return ret;
        else
            return ret + c - bb;
    }
};


// --------------------------------------------------------------------------
// Map: LexText -> Text                                              [String]
// --------------------------------------------------------------------------

template <typename TInput, typename TResult = TInput>
struct _dislexReverseTransform :
    public std::unary_function<TInput, TResult>
{
    const TInput _s, _n, _b, _c;

    _dislexReverseTransform(TInput s, TInput n) : _s(s), _n(n), _b(n / s), _c(n % s)
    {}

    inline TResult operator() (const TInput & x) const
    {
        TInput bb, r, ret;

        if (x < (_s - _c - 1) * _b) {
            bb = _s - 1 - (x / _b);
            r = x - (_s - bb - 1) * _b;
            ret = (r * _s) + _c - bb + _s ;
        } else {
            bb = _c - (x - (_s - _c - 1) * _b) / (_b + 1);
            r = x - (_s - bb - 1) * _b - _c + bb;
            ret = (r * _s) + _c - bb;
        }
        return static_cast<TResult>(ret);
    }
};

// --------------------------------------------------------------------------
// Map: LexText -> Text                                           [StringSet]
// --------------------------------------------------------------------------

template <typename TInput,                    // global pos
          typename TLimits,                   // limits
          typename TResult = Pair<TInput> >   // local pos
struct _dislexReverseTransformMulti :
    public std::unary_function<TInput, TResult>
{
    TLimits const & _limits;
    TInput const _s;

    _dislexReverseTransformMulti(TInput s, TLimits const & strSetLimits) : _limits(strSetLimits), _s(s)
    {}

    inline TResult operator() (const TInput & x) const
    {
        TResult ret;
        // binary search to find the corresponding sequence
        posLocalize(ret, x, _limits);
        return local2local(ret);
    }
    
    inline TResult local2local(TResult ret) const
    {
        TInput n = _limits[ret.i1 + 1] - _limits[ret.i1];
        TInput b = n / _s;
        TInput c = n % _s;
        TInput bb, r, i = ret.i2;
        if (i < (_s - c - 1) * b) {
            bb = _s - 1 - (i / b);
            r = i - (_s - bb - 1) * b;
            ret.i2 = (r * _s) + c - bb + _s ;
        } else {
            bb = c - (i - (_s - c - 1) * b) / (b + 1);
            r = i - (_s - bb - 1) * b - c + bb;
            ret.i2 = (r * _s) + c - bb;
        }
        return ret;
    }
}; 
} //namespace

#endif  // #ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_SA_DISLEX_H_
