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
    TInput const S,N,B,C;

    DislexTransform_(TInput S_, TInput N_) : S(S_), N(N_), B(N_/S_), C(N_%S_)
    {}

    inline TResult operator() (TInput const & x) const
    {
        TInput r = x/S;
        TInput b = (N-x)%S;
        TInput ret = (S-1-b)*B + r;
        if (b <= C)
            ret += C-b;
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
    TString const & limits;
    TResult const S;

    DislexTransformMulti_(TResult S_, TString const & stringSetLimits) : limits(stringSetLimits), S(S_)
    {}

    inline TResult operator() (const TValue & x) const
    {
        TResult seq = x.i1;
        TResult pos = x.i2;
        TResult N = limits[seq+1] - limits[seq];

        // TODO: Store the following values for each String?
        TResult B = N/S;
        TResult C = N%S;

        TResult r = (pos)/S;
        TResult b = (N-pos)%S;
        TResult ret = limits[seq] + (S-1-b)*B + r;
        if (b > C)  return ret;
        else        return ret + C-b;
    }
};


// --------------------------------------------------------------------------
// Map: LexText -> Text                                              [String]
// --------------------------------------------------------------------------

template <typename TInput, typename TResult = TInput>
struct _dislexReverseTransform :
    public std::unary_function<TInput, TResult>
{
    const TInput  S,N,B,C;

    _dislexReverseTransform(TInput S_, TInput N_) : S(S_), N(N_), B(N_/S_), C(N_%S_)
    {}

    inline TResult operator() (const TInput & x) const
    {
        TInput b,r,ret;

        if (x < (S-C-1)*B) {
            b = S-1- x/B;
            r = x -(S-b-1)*B;
            ret = r*S + C - b + S ;
        } else {
            b = C -(x-(S-C-1)*B)/(B+1);
            r = x-(S-b-1)*B -C + b;
            ret = r*S + C - b;
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
    TLimits const & limits;
    TInput const S;

    _dislexReverseTransformMulti(TInput S_, TLimits const & stringSetLimits) : limits(stringSetLimits), S(S_)
    {}

    inline TResult operator() (const TInput & x) const
    {
        TResult ret;
        // binary search to find the corresponding sequence
        posLocalize(ret, x, limits);
        return local2local(ret);
    }
    
    inline TResult local2local(TResult ret) const
    {
        TInput N = limits[ret.i1+1] - limits[ret.i1];
        TInput B = N/S;
        TInput C = N%S;
        TInput b,r, i = ret.i2;
        if (i < (S-C-1)*B) {
            b = S-1- i/B;
            r = i -(S-b-1)*B;
            ret.i2 = r*S + C - b + S ;
        } else {
            b = C -(i-(S-C-1)*B)/(B+1);
            r = i-(S-b-1)*B -C + b;
            ret.i2 = r*S + C - b;
        }
        return ret;
    }
};

 
} //namespace

#endif  // #ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_SA_DISLEX_H_
