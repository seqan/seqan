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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_SEQUENCE_CONTAINER_VIEW_ZIPPED_H_
#define INCLUDE_SEQAN_SEQUENCE_CONTAINER_VIEW_ZIPPED_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Tag to subclass basic container view.
template <typename TSpec = void>
struct ZipContainer;

// Use tuple to store variadic container types.
template <typename... TContPack, typename TSpec>
class ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> >
{
public:
    typedef typename Iterator<ContainerView, Standard>::Type TIterator;

    TIterator   _begin;
    TIterator   _end;

    ContainerView() : _begin(), _end()
    {
        static_assert(sizeof...(TContPack) > 0, "Requries at least one argument!");
    }

    template <typename... TOtherContPack>
    ContainerView(TOtherContPack && ...contArgs)
    {
        static_assert(sizeof...(TOtherContPack) > 0, "Requires at least one argument!");
        static_assert(sizeof...(TOtherContPack) == sizeof...(TContPack), "Number of parameters don't match!");

        // TODO: We need to assert if lenght of containers differ!

        _begin = makeZippedIterator(begin(contArgs, Standard())...);
        _end   = makeZippedIterator(  end(contArgs, Standard())...);
    }

    // Custom c'tor setting the iterators.
    ContainerView(TIterator const & begin, TIterator const & end) :
        _begin(begin),
        _end(end)
    {}

    // ------------------------------------------------------------------------
    // Operator =
    // ------------------------------------------------------------------------

    template <typename TOtherContainer>
    SEQAN_HOST_DEVICE inline
    ContainerView &
    operator= (TOtherContainer & other)
    {
        assign(*this, other);
        return *this;
    }

    // ------------------------------------------------------------------------
    // Operator []
    // ------------------------------------------------------------------------

    template <typename TPos>
    SEQAN_HOST_DEVICE inline
    typename Reference<ContainerView>::Type
    operator[] (TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    SEQAN_HOST_DEVICE inline
    typename GetValue<ContainerView>::Type
    operator[] (TPos pos) const
    {
        return getValue(*this, pos);
    }
};

// TODO(rrahn): Implement resizable ContainerViewZipper!

// ============================================================================
// Metafunctions
// ============================================================================

// TODO(rrahn): View
// TODO(rrahn): RemoveView
// TODO(rrahn): IsView

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

// Delegate to the Value MF of the ZipIterator.
template <typename... TContPack, typename TSpec>
struct Value<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> > >
{
    typedef typename Iterator<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> >, Standard>::Type TIter_;
    typedef typename Value<TIter_>::Type                                                                     Type;
};

template <typename... TContPack, typename TSpec>
struct Value<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> > const>
{
    typedef typename Iterator<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> > const, Standard>::Type TIter_;
    typedef typename Value<TIter_>::Type                                                                           Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename... TContPack, typename TSpec>
struct GetValue<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> > >
{
    typedef typename Iterator<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> >, Standard>::Type TIter_;
    typedef typename GetValue<TIter_>::Type                                                                  Type;
};

template <typename... TContPack, typename TSpec>
struct GetValue<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> > const>
{
    typedef typename Iterator<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> > const, Standard>::Type TIter_;
    typedef typename GetValue<TIter_>::Type                                                                        Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename... TContPack, typename TSpec>
struct Reference<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> > >
{
    typedef typename Iterator<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> >, Standard>::Type TIter_;
    typedef typename Reference<TIter_>::Type                                                                 Type;
};

template <typename... TContPack, typename TSpec>
struct Reference<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> > const>
{
    typedef typename Iterator<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> > const, Standard>::Type TIter_;
    typedef typename Reference<TIter_>::Type                                                                       Type;
};


// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <typename... TContPack, typename TSpec>
struct Difference<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> > >
{
    typedef typename Iterator<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> >, Standard>::Type TIter_;
    typedef typename Difference<TIter_>::Type                                                                Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename... TContPack, typename TSpec>
struct Size<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> > >
{
    typedef typename Iterator<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> >, Standard>::Type TIter_;
    typedef typename Size<TIter_>::Type                                                                      Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename... TContPack, typename TSpec>
struct Position<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> > >
{
    typedef typename Iterator<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> >, Standard>::Type TIter_;
    typedef typename Position<TIter_>::Type                                                                  Type;
};

// TODO(rrahn): IsSequence
// TODO(rrahn): Infix
// TODO(rrahn): Prefix
// TODO(rrahn): Suffix
// TODO(rrahn): View (Segment)

// Create Zipped Iterator.
template <typename... TContPack, typename TSpec>
struct Iterator<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> >, Standard>
{
    typedef Iter<std::tuple<typename Iterator<TContPack, Standard>::Type...>, ZipIterator> Type;
};

template <typename... TContPack, typename TSpec>
struct Iterator<ContainerView<std::tuple<TContPack...>, ZipContainer<TSpec> > const, Standard>
{
    typedef Iter<std::tuple<typename Iterator<TContPack const, Standard>::Type...>, ZipIterator> Type;
};

// ============================================================================
// Functions
// ============================================================================

// Uses default functions from base ContainerView.

// ----------------------------------------------------------------------------
// Function makeZippedView()
// ----------------------------------------------------------------------------

// Helper function to create zipped view from container pack.
template <typename... TContPack, typename TSpec>
inline ContainerView<std::tuple<typename std::remove_reference<TContPack>::type...>, ZipContainer<TSpec> >
makeZippedView(TContPack && ...contArgs,
               TSpec const)
{
    return ContainerView<std::tuple<typename std::remove_reference<TContPack>::type...>, ZipContainer<TSpec> >(std::forward<TContPack>(contArgs)...);
}

template <typename... TContPack>
inline ContainerView<std::tuple<typename std::remove_reference<TContPack>::type...>, ZipContainer<> >
makeZippedView(TContPack && ...contArgs)
{
    return ContainerView<std::tuple<typename std::remove_reference<TContPack>::type...>, ZipContainer<> >(std::forward<TContPack>(contArgs)...);
}

}

#endif  // #ifndef INCLUDE_SEQAN_SEQUENCE_CONTAINER_VIEW_ZIPPED_H_
