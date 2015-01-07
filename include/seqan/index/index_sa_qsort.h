/// ==========================================================================
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
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_SA_QSORT_H
#define SEQAN_HEADER_INDEX_SA_QSORT_H

namespace SEQAN_NAMESPACE_MAIN
{

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

struct SAQSort{};
    
template < typename TSAValue, typename TText, typename TMod = void >
struct SuffixLess_;

// --------------------------------------------------------------------------
// Class SuffixLess_                                     [String], [Standard]
// --------------------------------------------------------------------------

template < typename TSAValue, typename TText >
struct SuffixLess_<TSAValue, TText, void> :
    public ::std::binary_function < TSAValue, TSAValue, bool >
{
    typedef typename Iterator<TText const, Standard>::Type TIter;
    TIter _begin, _end;

    SuffixLess_(TText const &text):
        _begin(begin(text, Standard())),
        _end(end(text, Standard()))
    {}

    // skip the first <offset> characters
    template <typename TSize>
    SuffixLess_(TText const &text, TSize offset):
        _begin(begin(text, Standard()) + offset),
        _end(end(text, Standard()))
    {}

    inline bool operator() (TSAValue const a, TSAValue const b) const 
    {
        if (a == b) return false;
        TIter itA = _begin + a;
        TIter itB = _begin + b;
        if (a <= b)
        {
            for(; itB != _end; ++itB, ++itA) {
                if (ordLess(*itA, *itB)) return true;
                if (ordLess(*itB, *itA)) return false;
            }
            return false;
        } else {
            for(; itA != _end; ++itA, ++itB) {
                if (ordLess(*itA, *itB)) return true;
                if (ordLess(*itB, *itA)) return false;
            }
            return true;
        }
    }
};

// --------------------------------------------------------------------------
// Class SuffixLess_                                  [StringSet], [Standard]
// --------------------------------------------------------------------------

template < typename TSAValue, typename TString, typename TSetSpec >
struct SuffixLess_<TSAValue, StringSet<TString, TSetSpec> const, void > :
    public ::std::binary_function < TSAValue, TSAValue, bool >
{
    typedef StringSet<TString, TSetSpec> const TText;
    
    TText &_text;
    typename Size<TString>::Type _offset;

    SuffixLess_(TText &text):
        _text(text), _offset(0)
    {}

    // skip the first <offset> characters
    template <typename TSize>
    SuffixLess_(TText &text, TSize offset):
        _text(text),
        _offset(offset)
    {}
    
    inline bool operator() (TSAValue const a, TSAValue const b) const 
    {
        typedef typename Iterator<TString const, Standard>::Type TIter;
        if (a == b) return false;
        TIter itA = begin(getValue(_text, getSeqNo(a)), Standard()) + getSeqOffset(a) + _offset;
        TIter itB = begin(getValue(_text, getSeqNo(b)), Standard()) + getSeqOffset(b) + _offset;
        TIter itAEnd = end(getValue(_text, getSeqNo(a)), Standard());
        TIter itBEnd = end(getValue(_text, getSeqNo(b)), Standard());
        if (itAEnd - itA <= itBEnd - itB)
        {
            // a is shorter or equal to b
            for(; itA != itAEnd; ++itA, ++itB) {
                if (ordLess(*itA, *itB)) return true;
                if (ordLess(*itB, *itA)) return false;
            }
            // if a is really shorter than b, return true
            if (itB != itBEnd) return true;
            // else the higher sequence ID makes the smaller suffix:
            else               return getSeqNo(a) > getSeqNo(b);
        }
        else
        {
            // b is shorter than a
            for(; itB < itBEnd; ++itB, ++itA) {
                if (ordLess(*itA, *itB)) return true;
                if (ordLess(*itB, *itA)) return false;
            }
            return false;
        }
    }
};
    
/* old version. If mine breaks, I might need that code
 
		typename Size<TString>::Type _offset;
		TText &_text;

		SuffixLess_(TText &text):
			_text(text) {}
			
		// skip the first <offset> characters
		template <typename TSize>
		SuffixLess_(TText &text, TSize offset):
			_text(text),
			_offset(offset) {}
		
		inline bool operator() (TSAValue const a, TSAValue const b) const 
		{
			typedef typename Iterator<TString const, Standard>::Type TIter;
			if (a == b) return false;
			TIter itA = begin(getValue(_text, getSeqNo(a)), Standard()) + getSeqOffset(a);
			TIter itB = begin(getValue(_text, getSeqNo(b)), Standard()) + getSeqOffset(b);
			TIter itAEnd = end(getValue(_text, getSeqNo(a)), Standard());
			TIter itBEnd = end(getValue(_text, getSeqNo(b)), Standard());
			if (itAEnd - itA < itBEnd - itB) 
            {
				for(; itA != itAEnd; ++itA, ++itB) {
					if (ordLess(*itA, *itB)) return true;
					if (ordLess(*itB, *itA)) return false;
				}
				return true;
			} else {
				for(; itB != itBEnd; ++itB, ++itA) {
					if (ordLess(*itA, *itB)) return true;
					if (ordLess(*itB, *itA)) return false;
				}
                if (itA != itAEnd)
                    return false;
				return getSeqNo(a) > getSeqNo(b);
			}
		}	
	};
*/


// --------------------------------------------------------------------------
// SuffixLess_                                           [String], [Modified]
// --------------------------------------------------------------------------

template < typename TSAValue, typename TText, typename TSuffixMod>
struct SuffixLess_ :
    public ::std::binary_function < TSAValue, TSAValue, bool >
{
    typedef ModifiedString<typename Suffix<TText const>::Type, TSuffixMod>      TSuffix;
    typedef typename Cargo<TSuffix>::Type                                       TModCargo;
    typedef typename Iterator<TSuffix,Standard>::Type                           TSuffIter;

    TText const &               _text;
    TModCargo                   _modifier;
    typename Size<TText>::Type  _offset;

    SuffixLess_(TText const &text, TModCargo const & modifier):
        _text(text), _modifier(modifier), _offset(0)
    {}

    // skip the first <offset> characters
    template <typename TSize>
    SuffixLess_(TText const &text, TModCargo const & modifier, TSize offset):
        _text(text), _modifier(modifier), _offset(offset)
    {}


    inline bool operator() (TSAValue a, TSAValue b) const
    {
        if (a == b) return false;

        TSuffix sa(suffix(_text, a), _modifier);
        TSuffix sb(suffix(_text, b), _modifier);

        TSuffIter saIt = begin(sa, Standard()) + _offset;
        TSuffIter sbIt = begin(sb, Standard()) + _offset;

        TSuffIter saEnd = end(sa, Standard());
        TSuffIter sbEnd = end(sb, Standard());

        for (; saIt < saEnd && sbIt < sbEnd; ++saIt, ++sbIt)
        {
            if (ordLess(*saIt, *sbIt)) return true;
            if (ordLess(*sbIt, *saIt)) return false;
        }

        // if both suffixes are empty, the suff pos decides
        if (!(saIt < saEnd) && !(sbIt < sbEnd)) // NOTE: cannot write == here because of _offset
            return a > b;

        return sbIt < sbEnd; // a < b only if sb is not empty

    }
};

// --------------------------------------------------------------------------
// SuffixLess_                                        [StringSet], [Modified]
// --------------------------------------------------------------------------

template <typename TSAValue, typename TString, typename TSetSpec, typename TSuffixMod>
struct SuffixLess_<TSAValue, StringSet<TString, TSetSpec> const, TSuffixMod> :
    public ::std::binary_function<TSAValue, TSAValue, bool>
{
    typedef StringSet<TString, TSetSpec>                                    TText;
    typedef ModifiedString<typename Suffix<TText const>::Type, TSuffixMod>  TSuffix;
    typedef typename Cargo<TSuffix>::Type                                   TModCargo;
    typedef typename Iterator<TSuffix,Standard>::Type                       TSuffIter;

    TText const &               _text;
    TModCargo                   _modifier;
    typename Size<TText>::Type  _offset;

    SuffixLess_(TText const &text, TModCargo const & modifier):
        _text(text), _modifier(modifier), _offset(0)
    {}

    // skip the first <offset> characters
    template <typename TSize>
    SuffixLess_(TText const & text, TModCargo const & modifier, TSize offset):
        _text(text), _modifier(modifier), _offset(offset)
    {}

    inline bool operator() (TSAValue a, TSAValue b) const
    {
        if (a == b) return false;

        TSuffix sa(suffix(_text, a), _modifier);
        TSuffix sb(suffix(_text, b), _modifier);

        TSuffIter saIt = begin(sa, Standard()) + _offset;
        TSuffIter sbIt = begin(sb, Standard()) + _offset;

        TSuffIter saEnd = end(sa, Standard());
        TSuffIter sbEnd = end(sb, Standard());

        for (; saIt < saEnd && sbIt < sbEnd; ++saIt, ++sbIt)
        {
            if (ordLess(*saIt, *sbIt)) return true;
            if (ordLess(*sbIt, *saIt)) return false;
        }

        // if both suffixes are empty, at first the underlying length and then the seqId decide
        if (!(saIt < saEnd) && !(sbIt < sbEnd)) // NOTE: cannot write == here because of _offset
        {
            typename Size<TText>::Type lena = suffixLength(a, _text);
            typename Size<TText>::Type lenb = suffixLength(b, _text);

            if (lena == lenb)
                return getSeqNo(a) > getSeqNo(b);
            else
                return lena < lenb;
        }

        // a < b  if sb is not yet empty
        return sbIt < sbEnd;
    }
};


// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function _sortBucketQuickSort
// --------------------------------------------------------------------------

// TODO(meiers): Old function, does somebody use it anywhere?

template < typename TSA, 
        typename TText,
        typename TSize>
void _sortBucketQuickSort(
    TSA &sa,
    TText &text,
    TSize lcp)
{
SEQAN_CHECKPOINT
    // sort bucket with quicksort
    ::std::sort(
        begin(sa, Standard()), 
        end(sa, Standard()), 
        SuffixLess_<typename Value<TSA>::Type, TText, void>(text, lcp));
}

// --------------------------------------------------------------------------
// function _initializeSA()
// --------------------------------------------------------------------------

template <typename TSA, typename TString>
inline void _initializeSA(TSA & sa, TString const & /**/)
{
    typedef typename Iterator<TSA, Standard>::Type TIter;
    typedef typename Value<TSA>::Type TSAVal;
    TIter it = begin(sa, Standard());
    TIter itEnd = end(sa, Standard());
    for(TSAVal i = 0; it != itEnd; ++it, ++i)
        *it = i;
}

template <typename TSA, typename TString, typename TSpec>
inline void _initializeSA(TSA & sa, StringSet<TString, TSpec> const & strSet)
{
    typedef typename Iterator<TSA, Standard>::Type TIter;
    typedef typename Value<TSA>::Type TSAVal;
    typedef typename Size<TSA>::Type TSize; // TODO: derive Size Type from TSAVal, but how?

    TIter it = begin(sa, Standard());
    TIter itEnd = end(sa, Standard());
    TSize const setLen = length(strSet);
    for(TSize j = 0; j < setLen && it != itEnd; ++j)
    {
        TSize const len = length(strSet[j]);
        for(TSize i = 0; i < len && it != itEnd; ++i, ++it)
            *it = TSAVal(j, i);
    }
}

// --------------------------------------------------------------------------
// Function createSuffixArray
// --------------------------------------------------------------------------

template < typename TSA, typename TText>
inline void createSuffixArray(
    TSA &SA,
    TText const &s,
    SAQSort const &)
{
    typedef typename Size<TSA>::Type TSize;

    // 1. Fill suffix array with a permutation (the identity)
    _initializeSA(SA, s);

    // 2. Sort suffix array with quicksort
    ::std::sort(
        begin(SA, Standard()), end(SA, Standard()),
        SuffixLess_<typename Value<TSA>::Type, TText const, void>(s));
}

// --------------------------------------------------------------------------
// function createGappedSuffixArray
// --------------------------------------------------------------------------

/*!
 * @fn createGappedSuffixArray
 * @headerfile <seqan/index.h>
 * @brief Creates a gapped suffix array from a given text and a suffix modifier (e.g. a cyclic shape).
 *
 * @signature void createGappedSuffixArray(suffixArray, text, modCargo, modifierTag[, algoTag]);
 *
 * @param[out] suffixArray The resulting gapped suffix array.
 * @param[in]  text    A given text. Types: @link ContainerConcept @endlink
 * @param[in]  modCargo The Cargo of the Suffix modifier, i.e. a CyclicShape.
 * @param[in]  modfierTag The Suffix modifier class, e.g. ModCyclicShape<TShape>.
 * @param[in]  algoTag A tag that identifies the algorithm which is used for creation, e.g. SAQSort.
 *
 * This function should not be called directly and the documentation will be changed accordingly later on.
 * The suffix array must have the correct size. Here is the typical usage:
 *
 * @code{.cpp}
 * typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,2> >,0> > TShape; // 1101
 * typedef ModifiedString<typename Suffix<CharString>::Type, ModCyclicShape<TShape> > TSuffix;
 *
 * CharString string = ".a..b..c....a.";
 * String<int> sa;
 * resize(sa, length(string));
 *
 * createGappedSuffixArray(sa,
 *                         string,
 *                         TShape(),
 *                         ModCyclicShape<TShape>(),
 *                         SAQSort());
 * for (unsigned i=0; i< length(sa); ++i)
 * {
 *     std::cout << i << "\t" << sa[i] << "\t" <<
 *     TSuffix(suffix(string, sa[i]), TShape()) <<
 *     std::endl;
 * }
 * @endcode
 * @code{.console}
 * 0	13	.
 *  1	10	...
 *  2	5	.....a.
 *  3	2	....c....
 *  4	8	...a.
 *  5	9	..a.
 *  6	11	.a
 *  7	0	.a.b.c...a.
 *  8	3	.b.c...a
 *  9	6	.c....
 *  10	12	a.
 *  11	1	a.b.....a.
 *  12	4	b.c...a.
 *  13	7	c...a
 * @endcode
 */

// todo(meiers): Once the GappedIndex class is there, write a documentation for requireIndex etc.
// Moreover, write a tutorial for gapped suffix arrays.

template <typename TSA, typename TText, typename TCargo, typename TMod>
inline void createGappedSuffixArray(
    TSA &SA,
    TText const &s,
    TCargo const & modifierCargo,
    TMod const &,
    SAQSort const &)
{
    typedef typename Size<TSA>::Type TSize;
    typedef typename Iterator<TSA, Standard>::Type TIter;

    // 1. Fill suffix array with a permutation (the identity)
    _initializeSA(SA, s);

    // 2. Sort suffix array with quicksort
    ::std::sort(
        begin(SA, Standard()), end(SA, Standard()),
        SuffixLess_<typename Value<TSA>::Type, TText const, TMod>(s, modifierCargo));
}


// Old stuff:

    //////////////////////////////////////////////////////////////////////////////
    // suffix quicksort pipe
    template < typename TInput >
    struct Pipe< TInput, SAQSort >
    {
		typedef typename Value<TInput>::Type	TValue;
		typedef typename SAValue<TInput>::Type	TSAValue;

		typedef String<TValue, Alloc<> >		TText;
		typedef String<TSAValue, Alloc<> >		TSA;
		typedef Pipe<TSA, Source<> >			TSource;

		TSA		sa;
		TSource	in;

		Pipe(TInput &_textIn):
			in(sa)
		{
			TText text;
			text << _textIn;

			resize(sa, length(_textIn), Exact());
			createSuffixArray(sa, text, SAQSort());
		}

		inline typename Value<TSource>::Type const & operator*()
        {
            return *in;
        }
        
        inline Pipe& operator++()
        {
            ++in;
            return *this;
        }        
	};

}

#endif
