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



template <typename TString>
void printStr(TString const & str) { for(unsigned i=0; i< length(str); ++i)
    std::cout << str[i] << ","; std::cout << std::endl; }

template <typename T, typename A>
void printStrSet(StringSet<T,A> const & strSet) { for(unsigned i=0; i< length(strSet); ++i) {
    std::cout << "[" << i << "]"; printStr(strSet[i]); std::cout << std::endl; } }

template <typename TSACA = Skew7>
struct Dislex {};


// ==========================================================================
// DisLex Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction LexText
// --------------------------------------------------------------------------

template <typename TText, typename TMod>
struct LexText {
    typedef typename SAValue<TText>::Type TSize;
    typedef String<TSize> Type;
};
template <typename TText, typename TSpec, typename TMod>
struct LexText<StringSet<TText,TSpec>, TMod> : LexText<TText, TMod> {};


// ==========================================================================
// DisLex mapping functors
// ==========================================================================

// --------------------------------------------------------------------------
// Map: Text -> LexText                                              [String]
// --------------------------------------------------------------------------

// real dislex transformation:
// takes a position x and returns L(x)
template <
typename TInput,
typename TResult = TInput>
struct _dislexTransform :
public std::unary_function<TInput, TResult>
{
    TInput const S,N,B,C;

    _dislexTransform(TInput S_, TInput N_) : S(S_), N(N_), B(N_/S_), C(N_%S_)
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
template <
typename TValue,
typename TString,
typename TResult = typename Value<TValue, 2>::Type>
struct _dislexTransformMulti :
public std::unary_function<TValue, TResult>
{
    TString const & limits;
    TResult const S;

    _dislexTransformMulti(TResult S_, TString const & stringSetLimits) : limits(stringSetLimits), S(S_)
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

template <
    typename TInput,
    typename TResult = TInput>
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

template <
    typename TInput,                    // global pos
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




// ==========================================================================
// DisLex comparsion functors
// ==========================================================================

// NOTE(meiers):
//  - For a full suffix comparison use _SuffixLess (index_sa_qsort.h)
//  - For the comparison of suffixes of which we know are equal up to their
//    last character, use _ZeroBucketComparator (radix_inplace.h)
//  - Inside _dislex() we need to compare gapped k-mers. We use
//    GappedSuffixQgramLess_ functors for that.
//  - The external algorithm saves a part of the sequence into a tuple and
//    then uses _dislexTupleComp/Multi to compare them. These have overloads
//    for bitpacked tuples (index_gapped_dislex_external.h)


// --------------------------------------------------------------------------
// struct GappedSuffixQgramLess_                                     [String]
// --------------------------------------------------------------------------

template <typename TSAValue, typename TShape, typename TText, typename TResult=int>
struct GappedSuffixQgramLess_;

template <typename TSAValue, typename TShape, typename TText, typename TResult>
struct GappedSuffixQgramLess_ : public std::binary_function<TSAValue, TSAValue, TResult>
{
    typedef typename Size<TText>::Type                                          TSize;
    typedef ModifiedString<typename Suffix<TText const>::Type, ModCyclicShape<TShape> >      TSuffix;
    typedef typename Iterator<TSuffix,Standard>::Type                           TSuffIter;

    TText const &   _text;
    TShape const &  _shape;
    TSize const     _weight, _span, N;

    GappedSuffixQgramLess_(TText const &text, TShape const & shape, TSize weight):
    _text(text), _shape(shape), _weight(weight), _span(shape.span), N(length(text))
    {}

    inline int operator() (TSAValue a, TSAValue b) const
    {
        if (a == b) return 0;

        TSuffix sa(suffix(_text, a), _shape);
        TSuffix sb(suffix(_text, b), _shape);

        TSuffIter saIt = begin(sa, Standard());
        TSuffIter sbIt = begin(sb, Standard());

        TSuffIter saEnd = end(sa, Standard());
        TSuffIter sbEnd = end(sb, Standard());

        TSize p = 0;

        // lexicographical comparison
        for (; saIt < saEnd && sbIt < sbEnd && p < _weight; ++saIt, ++sbIt, ++p)
        {
            if (ordValue(*saIt) == ordValue(*sbIt)) continue;
            return (ordValue(*saIt) < ordValue(*sbIt) ? -1 : 1);
        }

        typename Size<TText>::Type lena = N - a;
        typename Size<TText>::Type lenb = N - b;

        // both cyclic shapes are "full"
        if (lena >= _shape.span && lenb >= _shape.span)
            return 0;

        // if they are equally long, start position decides
        if (sbEnd - sbIt == saEnd - saIt)
        {
            if (lena > lenb) return 1;
            if (lena < lenb) return -1;
            // should never occur:
            SEQAN_ASSERT_EQ(true,false);
            return 0;
        }

        // they are not equally long
        if (!(sbIt < sbEnd)) return 1;
        if (!(saIt < saEnd)) return -1;

        // should never occur:
        SEQAN_ASSERT_EQ(true, false);
        return 0;
    }
};

// --------------------------------------------------------------------------
// struct GappedSuffixQgramLess_                                  [StringSet]
// --------------------------------------------------------------------------
// Note(meiers):
// For string sets the comparsion has to consider more than weight characters!
// See my Master thesis for an explanation

template <typename TSAValue, typename TShape, typename TText, typename TSpec, typename TResult>
struct GappedSuffixQgramLess_ <TSAValue, TShape, StringSet<TText, TSpec>, TResult> :
    public std::binary_function<TSAValue, TSAValue, TResult>
{
    typedef StringSet<TText, TSpec>                         TSet;
    typedef typename Size<TText>::Type                      TSize;
    typedef ModifiedString<typename Suffix<TSet const>::Type, ModCyclicShape<TShape> >      TSuffix;
    typedef typename Iterator<TSuffix,Standard>::Type       TSuffIter;

    TSet const &   _text;
    TShape const &  _shape;
    TSize const     _weight;

    GappedSuffixQgramLess_(TSet const &text, TShape const & shape, TSize weight):
    _text(text), _shape(shape), _weight(weight)
    {}

    inline int operator() (TSAValue a, TSAValue b) const
    {
        if (a == b) return 0;

        TSuffix sa(suffix(_text, a), _shape);
        TSuffix sb(suffix(_text, b), _shape);

        TSuffIter saIt = begin(sa, Standard());
        TSuffIter sbIt = begin(sb, Standard());

        TSuffIter saEnd = end(sa, Standard());
        TSuffIter sbEnd = end(sb, Standard());

        TSize p = 0;

        // lexicographical comparison
        for (; saIt < saEnd && sbIt < sbEnd && p < _weight; ++saIt, ++sbIt, ++p)
        {
            if (ordValue(*saIt) == ordValue(*sbIt)) continue;
            return (ordValue(*saIt) < ordValue(*sbIt) ? -1 : 1);
        }

        typename Size<TSet>::Type lena = length(_text[getSeqNo(a)]) - getSeqOffset(a);
        typename Size<TSet>::Type lenb = length(_text[getSeqNo(b)]) - getSeqOffset(b);

        // both cyclic shapes are more than "full"
        if (lena > _shape.span && lenb > _shape.span)
            return 0;

        // if they are equally long, start position and seqId decide
        if (sbEnd - sbIt == saEnd - saIt)
        {
            if (lena == lenb)
            {
                if (getSeqNo(a) < getSeqNo(b)) return 1;
                if (getSeqNo(a) > getSeqNo(b)) return -1;
            } else {
                if (lena > lenb) return 1;
                if (lena < lenb) return -1;
            }
            // should never occur:
            SEQAN_ASSERT_EQ(true,false);
            return 0;
        }

        // they are not equally long
        if (!(sbIt < sbEnd)) return 1;
        if (!(saIt < saEnd)) return -1;

        // should never occur:
        SEQAN_ASSERT_EQ(true, false);
        return 0;
    }
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// helper function _maxSigma: determine an upper bound for the alphabet size
// --------------------------------------------------------------------------

template <typename TText>
typename Size<TText>::Type _maxSigma(TText const & text, typename Size<TText>::Type weight, typename Size<TText>::Type span)
{
    typedef typename Size<TText>::Type TSize;
    typedef typename Value<TText const>::Type TAlph;
    typedef typename Iterator<TText const>::Type TIter;
    
    // short texts
    if (length(text) == 0)
        return 1;
    if (length(text) < span)
        return length(text);
    
    String<TSize> alphabet;
    resize(alphabet, valueSize<TAlph>());
    std::fill(begin(alphabet), end(alphabet), 0u);    
    
    for(TIter it=begin(text); it!=end(text); ++it)
        ++alphabet[ordValue(*it)];
        
    TSize sigma = 0;
    for(unsigned i=0; i<length(alphabet); ++i)
        if(alphabet[i]) sigma++;
    
    TSize fullShapes = std::min( static_cast<TSize>(std::pow(double(sigma), double(weight))), length(text)-span+1);
    TSize incompleteShapes = span-1;
    
    return fullShapes + incompleteShapes;
}


// NOTE(meiers): The function MaxValue or  valueSize will fail on Finite alphabets
template <typename TText, typename TSpec>
typename Size<TText>::Type _maxSigma(StringSet<TText, TSpec> const & text, typename Size<TText>::Type weight, typename Size<TText>::Type span)
{
    typedef typename Size<TText>::Type TSize;
    typedef typename Value<typename Value<TText const>::Type>::Type TAlph;
    
    String<TSize> alphabet;
    resize(alphabet, valueSize<TAlph>());
    std::fill(begin(alphabet), end(alphabet), static_cast<TSize>(0));
    
    for(TSize i=0; i<length(text); ++i)
        for (TSize j=0; j<length(text[i]); ++j)
            ++alphabet[ordValue(text[i][j])];

    TSize sigma = 0;
    for(unsigned i=0; i<length(alphabet); ++i)
        if(alphabet[i]) sigma++;
     
    TSize incompleteShapes = 0;
    for (TSize i=0; i< length(text); ++i)
        incompleteShapes += std::min(length(text[i]), static_cast<TSize>(span));
        
    TSize fullShapes = std::min( static_cast<TSize>(std::pow(double(sigma), double(weight))), lengthSum(text)-incompleteShapes);
    
    return fullShapes + incompleteShapes;
}


// --------------------------------------------------------------------------
// function _dislex()                                                [String]
// --------------------------------------------------------------------------

template <typename TLexText, typename TSA, typename TText, typename TCyclicShape>
inline typename Value<TLexText>::Type _dislex(
    TLexText & lexText,                         // random access (write)
    TSA const & partiallyOrderedSA,             // sequential scan
    TText const & origText,                     // random access
    TCyclicShape const & cyclic)
{
    typedef typename Size<TSA>::Type                        TSize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIter;
    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Value<TLexText>::Type                  TRank;
    typedef ModifiedString<typename Suffix<TText const
    >::Type, ModCyclicShape<TCyclicShape> >      TModText;

    if (empty(partiallyOrderedSA)) return 0; // otherwise *sa will fail

    // dislex position calculator
    _dislexTransform<TSize> dislex(cyclic.span, length(origText));

    // q-gram comparator to determine rank
    GappedSuffixQgramLess_<TSAValue, TCyclicShape,TText> comp(origText, cyclic, weight(cyclic));

    resize(lexText, length(origText), Exact());

    TSAIter sa    = begin(partiallyOrderedSA, Standard());
    TSAIter saEnd = end(partiallyOrderedSA, Standard());
    TRank       rank    = 0;
    TSAValue    txtPos  = *sa;

    // scan along the SA
    for(++sa; sa != saEnd; ++sa)
    {
        // write rank to position in lexText
        lexText[dislex(txtPos)] = rank;

        // compare two consecutive values: this is probably slow
        if(comp(txtPos, *sa))
            ++rank;

        //std::cout << txtPos << "  ...   " << *sa << " -> " << comp(txtPos, *sa) << "\t\"" << TModText(suffix(origText, txtPos)) << "\", \"" << TModText(suffix(origText, *sa)) << "\"" << std::endl;
        SEQAN_ASSERT_GEQ(0, comp(txtPos,*sa));
        txtPos = *sa;
    }
    lexText[dislex(txtPos)] = rank;
    return rank;
}


// --------------------------------------------------------------------------
// function _dislex()                                             [StringSet]
// --------------------------------------------------------------------------

template < typename TLexText, typename TSA, typename TText, typename TTextSpec, typename TCyclicShape>
inline typename Value<typename Concatenator<TLexText>::Type>::Type _dislex(
    TLexText & lexText,                             // random access
    TSA const & partiallyOrderedSA,                 // sequential scan
    StringSet<TText, TTextSpec> const & origText,   // random access
    TCyclicShape const & cyclic)
{
    typedef typename Size<TSA>::Type                                    TSize;
    typedef typename Value<TSA>::Type                                   TSAValue;
    typedef typename Value<typename Concatenator<TLexText>::Type>::Type TRank;
    typedef typename Suffix<StringSet<TText, TTextSpec> const>::Type    TSuffix; // NOTE: This is != Suffix<TText>
    typedef ModifiedString<TSuffix, ModCyclicShape<TCyclicShape> >      TModText;
    typedef typename Iterator<TSA const, Standard>::Type                TSAIter;

    // position calculator
    typedef StringSet<TText, TTextSpec> const               TStringSet;
    typedef typename StringSetLimits<TStringSet>::Type      TStringSetLimits;   // expected: String<unsigned>
    typedef Pair<typename Size<TText>::Type>                TSetPosition;       // expected: Pair<unsigned>

    if (empty(partiallyOrderedSA)) return 0; // otherwise *sa will fail
    
    TStringSetLimits xxx = stringSetLimits(origText);
    _dislexTransformMulti<TSetPosition, TStringSetLimits>
    dislex(cyclic.span, xxx);


    // q-gram comparator to determine rank
    GappedSuffixQgramLess_<TSAValue, TCyclicShape, StringSet<TText, TTextSpec> > comp(origText, cyclic, static_cast<TSize>(weight(cyclic)));

    resize(lexText, lengthSum(origText));

    // scan along the SA
    TSAIter sa    = begin(partiallyOrderedSA, Standard());
    TSAIter saEnd = end(partiallyOrderedSA, Standard());

    TRank       rank = 0;
    TSAValue    txtPos = *sa;

    for(++sa; sa != saEnd; ++sa)
    {
        // write rank to position in lexText
        lexText[dislex(txtPos)] = rank;

        // compare two consecutive values, this is probably slow
        if(comp(txtPos, *sa))
            ++rank;

        //std::cout << txtPos << "  ...   " << *sa << " -> " << comp(txtPos, *sa) << "\t\"" << TModText(suffix(origText, txtPos)) << "\", \"" << TModText(suffix(origText, *sa)) << "\"" << std::endl;
        SEQAN_ASSERT_GEQ(0, comp(txtPos,*sa));
        txtPos = *sa;
    }
    lexText[dislex(txtPos)] = rank;
    return rank;
}


// --------------------------------------------------------------------------
// function _dislexReverse()                                         [String]
// --------------------------------------------------------------------------

template <typename TSA, typename TCyclicShape, typename TText, typename TLexSA>
void _dislexReverse(
    TSA & finalSA,                                  // random access
    TLexSA &,                                       // not needed
    TLexSA const & lexSA,                           // sequential scan
    TText const &,                                  // not needed
    TCyclicShape const & cyclic)
{
    typedef typename Iterator<TSA const, Standard>::Type    TLexSAIter;
    typedef typename Iterator<TSA, Standard>::Type          TSAIter;
    typedef typename Size<TSA>::Type                        TSize;

    _dislexReverseTransform<TSize> dislexRev(cyclic.span, length(lexSA));

    TLexSAIter sa       = begin(lexSA, Standard());
    TLexSAIter saEnd    = end(lexSA, Standard());
    TSAIter insert      = begin(finalSA, Standard());

    for(; sa < saEnd; ++sa, ++insert)
        *insert = dislexRev (*sa);
}

// --------------------------------------------------------------------------
// function _dislexReverse()                                      [StringSet]
// --------------------------------------------------------------------------

// non-linear-time version
template <typename TSA, typename TLexSA, typename TCyclicShape, typename TText, typename TTextSpec>
inline void _dislexReverse(
    TSA & finalSA,                                  // random access
    TLexSA &,                                       // only needed in linear time variant
    TLexSA const & lexSA,                           // sequential scan
    StringSet<TText, TTextSpec> const & origText,
    TCyclicShape const & cyclic,
    True const &)
{
    typedef typename Iterator<TLexSA const, Standard>::Type TLexSAIter;
    typedef typename Iterator<TSA, Standard>::Type          TSAIter;
    typedef typename Size<TSA>::Type                        TSize;

    typedef StringSet<TText, TTextSpec> const               TStringSet;
    typedef typename StringSetLimits<TStringSet>::Type      TStringSetLimits;   // expected: String<unsigned>
 
    _dislexReverseTransformMulti<TSize, TStringSetLimits>
    dislexRev(cyclic.span, stringSetLimits(origText));

    TLexSAIter sa       = begin(lexSA, Standard());
    TLexSAIter saEnd    = end(lexSA, Standard());
    TSAIter insert      = begin(finalSA, Standard());

    for(; sa < saEnd; ++sa, ++insert)
        *insert = dislexRev (*sa);
}

// linear-time version
template <typename TSA, typename TLexSA, typename TCyclicShape, typename TText, typename TTextSpec>
inline void _dislexReverse(
    TSA & finalSA,                                  // random access
    TLexSA & lexText,                               // only needed in linear time variant
    TLexSA const & lexSA,                           // sequential scan
    StringSet<TText, TTextSpec> const & origText,
    TCyclicShape const & cyclic,
    False const &)
{
    typedef typename Iterator<TLexSA const, Standard>::Type TLexSAIter;
    typedef typename Iterator<TSA, Standard>::Type          TSAIter;
    typedef typename Size<TSA>::Type                        TSize;

    typedef StringSet<TText, TTextSpec> const               TStringSet;
    typedef typename StringSetLimits<TStringSet>::Type      TStringSetLimits;

    // Prepare reverse mapping functor
    _dislexReverseTransformMulti<TSize, TStringSetLimits>
    dislexRev(cyclic.span, stringSetLimits(origText));
    
    // inverse lexSA into lexText
    for (TSize i=0; i< length(lexSA); ++i)
        lexText[lexSA[i]] = i;
    
    // Iterate along inverse suffix array using the Pairincrementer!
    PairIncrementer_<typename Value<TSA>::Type, TStringSetLimits> localPos;
    setHost(localPos, stringSetLimits(origText));
    TLexSAIter isaIt       = begin(lexText, Standard());
    TLexSAIter isaEnd      = end(lexText, Standard());
    
    for (; isaIt != isaEnd; ++isaIt, ++localPos) 
        finalSA[*isaIt] = dislexRev.local2local(value(localPos));
}

template <typename TSA, typename TLexSA, typename TCyclicShape, typename TText, typename TTextSpec>
inline void _dislexReverse(
    TSA & finalSA,                                  // random access
    TLexSA & lexText,                               // only needed in linear time variant
    TLexSA const & lexSA,                           // sequential scan
    StringSet<TText, TTextSpec> const & origText,
    TCyclicShape const & cyclic)
{
    // if the set consists only out of a single string, use pos localize
    if (length(origText) == 1u)
        _dislexReverse(finalSA, lexText, lexSA, origText, cyclic, True());

    // otherwise build the inverse suffix array
    else
        _dislexReverse(finalSA, lexText, lexSA, origText, cyclic, False());
}



// --------------------------------------------------------------------------
// function createGappedSuffixArray()                                [String]
// --------------------------------------------------------------------------
// Only defined for CyclicShapes

template < typename TSA, typename TText, typename TCyclicShape, typename TSACA>
inline void createGappedSuffixArray(
    TSA &SA, // must already be resized already
    TText const &s,
    TCyclicShape const & shape,
    ModCyclicShape<TCyclicShape> const &,
    Dislex<TSACA> const &)
{
    typedef typename LexText<TText, TCyclicShape>::Type         TLexText;

    #ifdef DISLEX_INTERNAL_RUNNING_TIMES
    double teim = sysTime();
    #endif

    // insert positions into SA
    _initializeSA(SA, s);

    #ifdef DISLEX_INTERNAL_RUNNING_TIMES
    std::cout << "   |     init: " << sysTime() - teim << "s" << std::endl; teim = sysTime();
    #endif

    // sort newSA according to the Shape
    inplaceRadixSort(SA, s, weight(shape)+1, shape, ModCyclicShape<TCyclicShape>());

    #ifdef DISLEX_INTERNAL_RUNNING_TIMES
    std::cout << "   | radix[" << (int)weight(shape) << "]: " << sysTime() - teim << "s" << std::endl; teim = sysTime();
    #endif

    // disLexTransformation
    TLexText lexText;
    typename Value<TLexText>::Type sigma = _dislex(lexText, SA, s, shape)+1u;

    #ifdef DISLEX_INTERNAL_RUNNING_TIMES
    typename Value<TLexText>::Type maxSigma = _maxSigma(s,  static_cast<typename Value<TLexText>::Type>(weight(shape)),
                                                            static_cast<typename Value<TLexText>::Type>(shape.span));
    std::cout << "   |   dislex: " << sysTime() - teim << "s\t\tsigma = " << sigma << " (" << maxSigma << ")" << std::endl; teim = sysTime();
    SEQAN_ASSERT_GEQ_MSG(maxSigma, sigma, "Alphabet size of lecical names exceeds the theoretical upper bound");
    #endif

    // Build Index using TSACA into the memory of SA.
    // TODO(meiers): Some createSuffixArray() functions do not accept 5 arguments
    createSuffixArray(SA, lexText, TSACA(), sigma, 0);

    #ifdef DISLEX_INTERNAL_RUNNING_TIMES
    std::cout << "   |     saca: " << sysTime() - teim << "s (len = " << length(concat(lexText)) << ")" << std::endl; teim = sysTime();
    #endif

    // reverse Transform of Index:
    // Note: 2nd argument not used in the specialization for String
    _dislexReverse(SA, SA, SA, s, shape) ;

    #ifdef DISLEX_INTERNAL_RUNNING_TIMES
    std::cout << "   |  reverse: " << sysTime() - teim << "s (len = " << length(SA) << ")" << std::endl; teim = sysTime();
    #endif
}

// --------------------------------------------------------------------------
// function createGappedSuffixArray()                             [StringSet]
// --------------------------------------------------------------------------

template < typename TSA, typename TText, typename TTextSpec, typename TCyclicShape, typename TSACA>
inline void createGappedSuffixArray(
    TSA &SA,                             	
    StringSet<TText, TTextSpec> const & s,
    TCyclicShape const & shape,
    ModCyclicShape<TCyclicShape> const &,
    Dislex<TSACA> const &)
{
    typedef typename LexText<TText, TCyclicShape>::Type         TLexText;

    #ifdef DISLEX_INTERNAL_RUNNING_TIMES
    double teim = sysTime();
    #endif

    // insert positions into SA
    _initializeSA(SA, s);

    #ifdef DISLEX_INTERNAL_RUNNING_TIMES
    std::cout << "   |     init: " << sysTime() - teim << "s" << std::endl; teim = sysTime();
    #endif

    // sort newSA according to the Shape
    inplaceRadixSort(SA, s, weight(shape)+1, shape, ModCyclicShape<TCyclicShape>());

    #ifdef DISLEX_INTERNAL_RUNNING_TIMES
    std::cout << "   | radix[" << (int)weight(shape) << "]: " << sysTime() - teim << "s" << std::endl; teim = sysTime();
    #endif

    // disLexTransformation
    TLexText lexText;
    typename Value<TLexText>::Type sigma = _dislex(lexText, SA, s, shape)+1u;

    #ifdef DISLEX_INTERNAL_RUNNING_TIMES
    typename Value<TLexText>::Type maxSigma = _maxSigma(s, weight(shape), shape.span);
    std::cout << "   |   dislex: " << sysTime() - teim << "s\t\tsigma = " << sigma << " (" << maxSigma << ")" << std::endl; teim = sysTime();
    SEQAN_ASSERT_GEQ_MSG(maxSigma, sigma, "Alphabet size of lecical names exceeds the theoretical upper bound");
    #endif

    // Build Index using Skew7
    TLexText innerSa;
    resize(innerSa, length(SA), Exact());
    createSuffixArray(innerSa, lexText, TSACA(), sigma, 0);

    #ifdef DISLEX_INTERNAL_RUNNING_TIMES
    std::cout << "   |     saca: " << sysTime() - teim << "s (len = " << length(concat(lexText)) << ")" << std::endl; teim = sysTime();
    #endif

    // reverse Transform of Index:
    _dislexReverse(SA, lexText, innerSa, s, shape) ;

    #ifdef DISLEX_INTERNAL_RUNNING_TIMES
    std::cout << "   |  reverse: " << sysTime() - teim << "s (len = " << length(innerSa) << ")" << std::endl; teim = sysTime();
    #endif
}

} //namespace

#endif  // #ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_GAPPED_SA_DISLEX_H_
