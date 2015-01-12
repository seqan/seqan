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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// Adaptions for STL arrays to SeqAn strings.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS: No forwards are generated for this file.

#ifndef SEQAN_SEQUENCE_ADAPT_STD_ARRAY_H_
#define SEQAN_SEQUENCE_ADAPT_STD_ARRAY_H_

#include <array>

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ===========================================================================
// Concepts
// ===========================================================================

// ----------------------------------------------------------------------------
// Concept StringConcept
// ----------------------------------------------------------------------------

template <typename TChar, std::size_t N>
SEQAN_CONCEPT_IMPL((std::array<TChar, N>), (StringConcept));          // resizable container

template <typename TChar, std::size_t N>
SEQAN_CONCEPT_IMPL((std::array<TChar, N> const), (ContainerConcept)); // read-only container

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename TChar, std::size_t N>
struct IsContiguous< std::array<TChar, N> >
{
    enum { VALUE = true };
};

template <typename  TChar, std::size_t N>
struct IsContiguous< std::array<TChar, N> const>
        : IsContiguous< std::array<TChar, N> > {};

template <typename TChar, std::size_t N>
struct Value< std::array<TChar, N> >
{
    typedef typename std::array<TChar, N>::reference Type;
};

template <typename TChar, std::size_t N>
struct Value< std::array<TChar, N> const>
{
    typedef typename std::array<TChar, N>::const_reference Type;
};

template <typename TChar, std::size_t N>
struct GetValue< std::array<TChar, N> >
{
    typedef typename std::array<TChar, N>::value_type Type;
};

template <typename TChar, std::size_t N>
struct GetValue< std::array<TChar, N> const>
{
    typedef typename std::array<TChar, N>::value_type const Type ;
};

template <typename TChar, std::size_t N>
struct Reference< std::array<TChar, N> >
{
    typedef typename std::array<TChar, N>::reference Type;
};

template <typename TChar, std::size_t N>
struct Reference< std::array<TChar, N> const>
{
    typedef typename std::array<TChar, N>::const_reference Type;
};

template <typename TChar, std::size_t N>
struct Iterator< std::array<TChar, N>, Rooted>
{
    typedef std::array<TChar, N> TArray_;
    typedef Iter<TArray_, StdIteratorAdaptor> TIterator_;
    typedef Iter<TArray_, AdaptorIterator<TIterator_> > Type;
};

template <typename TChar, std::size_t N>
struct Iterator< std::array<TChar, N> const, Rooted>
{
    typedef std::array<TChar, N> const TArray_;
    typedef Iter<TArray_, StdIteratorAdaptor> TIterator_;
    typedef Iter<TArray_, AdaptorIterator<TIterator_> > Type;
};

template <typename TChar, std::size_t N>
struct Iterator< std::array<TChar, N>, Standard >
{
    typedef Iter< std::array<TChar, N>, StdIteratorAdaptor > Type;
};

template <typename TChar, std::size_t N>
struct Iterator< std::array<TChar, N> const, Standard>
{
    typedef Iter< std::array<TChar, N> const, StdIteratorAdaptor > Type;
};

template <typename TChar, std::size_t N>
struct Position< std::array<TChar, N> >
{
    typedef typename std::array<TChar, N>::size_type Type;
};

template <typename TChar, std::size_t N>
struct Position< std::array<TChar, N> const>
        : Position< std::array<TChar, N> > {};

template <typename TChar, std::size_t N>
struct Size< std::array<TChar, N> >
{
    typedef typename std::array<TChar, N>::size_type Type;
};

template <typename TChar, std::size_t N>
struct Size< std::array<TChar, N> const>
        : Size< std::array<TChar, N> > {};

template <typename TChar, std::size_t N>
struct DefaultOverflowImplicit< std::array<TChar, N> >
{
    typedef Generous Type;
};

template <typename TChar, std::size_t N>
struct StdContainerIterator< std::array<TChar, N> >
{
    typedef std::array<TChar, N> TContainer_;
    typedef typename TContainer_::iterator Type;
};

template <typename TChar, std::size_t N>
struct StdContainerIterator< std::array<TChar, N> const>
{
    typedef std::array<TChar, N> TContainer_;
    typedef typename TContainer_::const_iterator Type;
};

template <typename TChar, std::size_t N>
struct IsSequence<std::array<TChar, N> > : True {};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TChar, std::size_t N>
inline void const *
getObjectId(std::array<TChar, N> const & me)
{
    SEQAN_CHECKPOINT;
    if (me.empty())
        return NULL;
    else
        return (& *(me.end() - 1)) + 1;
}

template <typename TChar, std::size_t N>
inline typename Iterator< std::array<TChar, N>, Standard>::Type
begin(std::array<TChar, N> & me,
      Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< std::array<TChar, N>, Standard>::Type(me.begin());
}
template <typename TChar, std::size_t N>
inline typename Iterator< std::array<TChar, N> const, Standard>::Type
begin(std::array<TChar, N> const & me,
      Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< std::array<TChar, N> const, Standard>::Type(me.begin());
}

template <typename TChar, std::size_t N>
inline typename Iterator< std::array<TChar, N>, Standard>::Type
end(std::array<TChar, N> & me,
    Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< std::array<TChar, N>, Standard>::Type(me.end());
}
template <typename TChar, std::size_t N>
inline typename Iterator< std::array<TChar, N> const, Standard>::Type
end(std::array<TChar, N> const & me,
    Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< std::array<TChar, N> const, Standard>::Type(me.end());
}

template <typename TChar, std::size_t N, typename TPos>
inline typename GetValue< std::array<TChar, N> >::Type
value(std::array<TChar, N> & me,
      TPos pos)
{
    SEQAN_CHECKPOINT;
    return me[pos];
}
template <typename TChar, std::size_t N, typename TPos>
inline typename GetValue< std::array<TChar, N> const>::Type
value(std::array<TChar, N> const & me,
      TPos pos)
{
    SEQAN_CHECKPOINT;
    return me[pos];
}

template <typename TChar, std::size_t N>
inline typename Size< std::array<TChar, N> >::Type
length(std::array<TChar, N> const & me)
{
    SEQAN_CHECKPOINT;
    return me.size();
}

template <typename TChar, std::size_t N>
inline typename Size< std::array<TChar, N> >::Type
capacity(std::array<TChar, N> const & me)
{
    SEQAN_CHECKPOINT;
    return me.max_size();
}

template <typename TChar, std::size_t N>
inline bool
empty(std::array<TChar, N> const & me)
{
    SEQAN_CHECKPOINT;
    return me.empty();
}

// template <typename TChar, std::size_t N>
// inline void
// clear(std::array<TChar, N> & me)
// {
//     SEQAN_CHECKPOINT;
//     me.clear();
// }

template <typename TChar, std::size_t N>
inline typename Reference<std::array<TChar, N> >::Type
front(std::array<TChar, N> & list)
{
    SEQAN_CHECKPOINT;
    return list.front();
}

template <typename TChar, std::size_t N>
inline typename Reference<std::array<TChar, N> const>::Type
front(std::array<TChar, N> const & list)
{
    SEQAN_CHECKPOINT;
    return list.front();
}

template <typename TChar, std::size_t N>
inline typename Reference<std::array<TChar, N> >::Type
back(std::array<TChar, N> & list)
{
    SEQAN_CHECKPOINT;
    return list.back();
}

template <typename TChar, std::size_t N>
inline typename Reference<std::array<TChar, N> const>::Type
back(std::array<TChar, N> const & list)
{
    SEQAN_CHECKPOINT;
    return list.back();
}

//////////////////////////////////////////////////////////////////////////////
//assign to std::array

template <typename TChar, std::size_t N, typename TSource>
inline void
assignImpl(std::array<TChar, N> & target,
           TSource const & source)
{
    std::size_t range = std::min(length(target), length(source));
    copy(begin(target), begin(target) + range,
         begin(source), begin(source) + range);
}

template <typename TChar, std::size_t N, typename TSize, typename TSource>
inline void
assignImplLimit(std::array<TChar, N> & target,
                TSource const & source,
                TSize limit)
{
    std::size_t range = std::min({static_cast<TSize>(length(target)),
                                  static_cast<TSize>(length(source)),
                                  limit});
    copy(begin(target), begin(target) + range,
         begin(source), begin(source) + range);
}

template <typename TChar, std::size_t N, typename TSource>
inline void
assign(std::array<TChar, N> & target,
       TSource & source)
{
    assignImpl(target, source);
}

template <typename TChar, std::size_t N, typename TSource>
inline void
assign(std::array<TChar, N> & target,
       TSource const & source)
{
    assignImpl(target, source);
}

template <typename TChar, std::size_t N, typename TSource, typename TSize>
inline void
assign(std::array<TChar, N> & target,
       TSource & source,
       TSize limit)
{
    assignImplLimit(target, source, limit);
}

template <typename TChar, std::size_t N, typename TSource, typename TSize>
inline void
assign(std::array<TChar, N> & target,
       TSource const & source,
       TSize limit)
{
    assignImplLimit(target, source, limit);
}

//____________________________________________________________________________

template <typename TChar, std::size_t N, typename TSource>
inline void
assign(std::array<TChar, N> & target,
       TSource & source,
       Generous)
{
    assignImpl(target, source);
}

template <typename TChar, std::size_t N, typename TSource>
inline void
assign(std::array<TChar, N> & target,
       TSource const & source,
       Generous)
{
    assignImpl(target, source);
}

template <typename TChar, std::size_t N, typename TSource>
inline void
assign(std::array<TChar, N> & target,
       TSource & source,
       typename Size< std::array<TChar, N> >::Type limit,
       Generous)
{
    assignImplLimit(target, source, limit);
}
template <typename TChar, std::size_t N, typename TSource>
inline void
assign(std::array<TChar, N> & target,
       TSource const & source,
       typename Size< std::array<TChar, N> >::Type limit,
       Generous)
{
    assignImplLimit(target, source, limit);
}

//____________________________________________________________________________

template <typename TChar, std::size_t N, typename TSource>
inline void
assign(std::array<TChar, N> & target,
       TSource & source,
       Limit)
{
    assignImpl(target, source);
}

template <typename TChar, std::size_t N, typename TSource>
inline void
assign(std::array<TChar, N> & target,
       TSource const & source,
       Limit)
{
    assignImpl(target, source);
}

template <typename TChar, std::size_t N, typename TSource>
inline void
assign(std::array<TChar, N> & target,
       TSource & source,
       typename Size< std::array<TChar, N> >::Type limit,
       Limit)
{
    assignImplLimit(target, source, limit);
}

template <typename TChar, std::size_t N, typename TSource>
inline void
assign(std::array<TChar, N> & target,
       TSource const & source,
       typename Size< std::array<TChar, N> >::Type limit,
       Limit)
{
    assignImplLimit(target, source, limit);
}


//////////////////////////////////////////////////////////////////////////////
//append to std::array
/*
template <typename TChar, std::size_t N, typename TSource>
inline void
append(std::array<TChar, N> & target,
       TSource const & source,
       Generous)
{
    SEQAN_CHECKPOINT;
    target.insert(target.end(), begin(source, Standard()), end(source, Standard()));
}

template <typename TChar, std::size_t N, typename TSource>
inline void
append(std::array<TChar, N> & target,
       TSource const & source,
       typename Size< std::array<TChar, N> >::Type limit,
       Generous)
{
    SEQAN_CHECKPOINT;
    typename Size< std::array<TChar, N> >::Type target_length = target.length();
    if (target_length > limit)
    {
        target.resize(limit);
    }
    else
    {
        limit -= target_length;
        typename Iterator<TSource const, Standard>::Type source_begin = begin(source, Standard());
        typename Size<TSource const>::Type source_length = length(source);
        if (source_length > limit)
        {
            source_length = limit;
        }

        target.insert(target.end(), source_begin, source_begin + source_length);
    }
}

//____________________________________________________________________________

template <typename TChar, std::size_t N, typename TSource>
inline void
append(std::array<TChar, N> & target,
       TSource const & source,
       Limit)
{
    SEQAN_CHECKPOINT;
    append(target, source, target.capacity(), Generous());
}

template <typename TChar, std::size_t N, typename TSource>
inline void
append(std::array<TChar, N> & target,
       TSource const & source,
       typename Size< std::array<TChar, N> >::Type limit,
       Limit)
{
    SEQAN_CHECKPOINT;
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }

    append(target, source, limit, Generous());
}

//////////////////////////////////////////////////////////////////////////////
template <typename TChar, std::size_t N, typename TValue, typename TTag>
inline void
appendValue(std::array<TChar, N> & me,
            TValue const & _value,
            TTag)
{
    SEQAN_CHECKPOINT;
    me.push_back(_value);
}

template <typename TChar, std::size_t N, typename TValue>
inline void
appendValue(std::array<TChar, N> & me,
            TValue const & _value,
            Limit)
{
    SEQAN_CHECKPOINT;
    if (capacity(me) > length(me)) me.push_back(_value);
}
*/
//////////////////////////////////////////////////////////////////////////////
//replace to std::array

template <typename TChar, std::size_t N, typename TSource>
inline void
replace(std::array<TChar, N> & target,
        typename Position< std::array<TChar, N> >::Type pos_begin,
        typename Position< std::array<TChar, N> >::Type pos_end,
        TSource const & source,
        Generous)
{
    std::size_t range = std::min(static_cast<std::size_t>(pos_end - pos_begin),
                                 static_cast<std::size_t>(length(source)));
    copy(begin(target) + pos_begin, begin(target) + pos_begin + range,
         begin(source), begin(source) + range);
}

template <typename TChar, std::size_t N, typename TSource>
inline void
replace(std::array<TChar, N> & target,
        typename Position< std::array<TChar, N> >::Type pos_begin,
        typename Position< std::array<TChar, N> >::Type pos_end,
        TSource const & source,
        typename Size< std::array<TChar, N> >::Type limit,
        Generous)
{
    std::size_t range = std::min({static_cast<std::size_t>(pos_end - pos_begin),
                                  static_cast<std::size_t>(length(source)),
                                  limit});
    copy(begin(target) + pos_begin, begin(target) + pos_begin + range,
         begin(source), begin(source) + range);
}

template <typename TChar, std::size_t N, typename TSource>
inline void
replace(std::array<TChar, N> & target,
        typename Position< std::array<TChar, N> >::Type pos_begin,
        typename Position< std::array<TChar, N> >::Type pos_end,
        TSource const & source,
        Limit)
{
    replace(target, pos_begin, pos_end, source, target.capacity(), Generous());
}

template <typename TChar, std::size_t N, typename TSource>
inline void
replace(std::array<TChar, N> & target,
        typename Position< std::array<TChar, N> >::Type pos_begin,
        typename Position< std::array<TChar, N> >::Type pos_end,
        TSource const & source,
        typename Size< std::array<TChar, N> >::Type limit,
        Limit)
{
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }
    replace(target, pos_begin, pos_end, source, limit, Generous());
}


//////////////////////////////////////////////////////////////////////////////
// handling of iterators as begin and end

template<typename TChar, typename TCharTraits, std::size_t N, typename TSource, typename TExpand>
inline void
replace(std::array<TChar, N> & target,
        typename Iterator< std::array<TChar, N>, Rooted>::Type pos_begin,
        typename Iterator< std::array<TChar, N>, Rooted>::Type pos_end,
        TSource & source,
        Tag<TExpand> tag)
{
    replace(target, position(pos_begin), position(pos_end), source, tag);
}

/*
template<typename TChar, std::size_t N, typename TSource, typename TExpand>
inline void
replace(std::array<TChar, N> & target,
        typename Iterator< std::array<TChar, N>, Rooted>::Type pos_begin,
        typename Iterator< std::array<TChar, N>, Rooted>::Type pos_end,
        TSource & source,
        typename Size< std::array<TChar, N> >::Type limit,
        Tag<TExpand> tag)
{
    replace(target,  position(pos_begin),  position(pos_end), source, tag);
}
*/


template <typename TChar, std::size_t N, typename TSize, typename TExpand>
inline typename Size< std::array<TChar, N> >::Type
reserve(std::array<TChar, N> & seq,
        TSize,
        Tag<TExpand>)
{
    // do nothing
    return capacity(seq);
}

template <typename TChar, std::size_t N, typename TSize>
inline typename Size< std::array<TChar, N> >::Type
reserve(std::array<TChar, N> & seq,
        TSize,
        Insist const &)
{
    // do nothing
    return capacity(seq);
}

template <typename TChar, std::size_t N, typename TSize>
inline typename Size< std::array<TChar, N> >::Type
reserve(std::array<TChar, N> & seq,
        TSize new_capacity,
        Limit const &)
{
    // do nothing
    return capacity(seq);
}

template <typename TChar, std::size_t N, typename TSize, typename TExpand>
inline typename Size< std::array<TChar, N> >::Type
resize(std::array<TChar, N> & me,
       TSize new_length,
       Tag<TExpand>)
{
    // do nothing
    return length(me);
}

template <typename TChar, std::size_t N, typename TSize, typename TExpand>
inline typename Size< std::array<TChar, N> >::Type
fill(std::array<TChar, N> & me,
    TSize new_length,
    TChar const & val,
    Tag<TExpand>)
{
    me.fill(val);
    return length(me);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_ADAPT_STD_ARRAY_H_
