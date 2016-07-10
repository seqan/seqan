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

#ifndef INCLUDE_SEQAN_SEQUENCE_CONTAINER_VIEW_ZIP_H_
#define INCLUDE_SEQAN_SEQUENCE_CONTAINER_VIEW_ZIP_H_

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

/*!
 * @class ZipContainerView
 * @extends ContainerView
 * @headerfile <seqan/sequence.h>
 *
 * @brief Ties one or more containers into a non-resizable container view together.
 *
 * @signature template <typename... TContTypes, typename TSpec>
              class ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> >
 *
 * @tparam TContTypes One or more types of all containers tied into the zip container view.
 * @tparam TSpec The specialization type.
 *
 * @section Example
 * 
 * The following example demonstrates the functionality of the @link ZipContainerView @endlink:
 *
 * @include demos/dox/sequence/container_view_zip.cpp
 * 
 * This outputs the following to the console:
 * @include demos/dox/sequence/container_view_zip.cpp.stdout
 */

// Use tuple to store variadic container types.
template <typename... TContTypes, typename TSpec>
class ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> >
{
public:
    typedef typename Iterator<ContainerView, Standard>::Type TIterator;

    TIterator   _begin;
    TIterator   _end;

    ContainerView() : _begin(), _end()
    {
        static_assert(sizeof...(TContTypes) > 0, "Requries at least one argument!");
    }

    template <typename... TOtherContPack>
    ContainerView(TOtherContPack && ...contArgs)
    {
        static_assert(sizeof...(TOtherContPack) > 0, "Requires at least one argument!");
        static_assert(sizeof...(TOtherContPack) == sizeof...(TContTypes), "Number of parameters don't match!");

        // TODO(rrahn): Check if all containers contain same number of elements.

        _begin = makeZipIterator(begin(contArgs, Standard())...);
        _end   = makeZipIterator(  end(contArgs, Standard())...);
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
    inline ContainerView &
    operator= (TOtherContainer & other)
    {
        assign(*this, other);
        return *this;
    }

    // ------------------------------------------------------------------------
    // Operator []
    // ------------------------------------------------------------------------

    template <typename TPos>
    inline typename Reference<ContainerView>::Type
    operator[] (TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename GetValue<ContainerView>::Type
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
template <typename... TContTypes, typename TSpec>
struct Value<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> > >
{
    typedef typename Iterator<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> >, Standard>::Type TIter_;
    typedef typename Value<TIter_>::Type                                                                     Type;
};

template <typename... TContTypes, typename TSpec>
struct Value<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> > const>
{
    typedef typename Iterator<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> > const, Standard>::Type TIter_;
    typedef typename Value<TIter_>::Type                                                                           Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename... TContTypes, typename TSpec>
struct GetValue<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> > >
{
    typedef typename Iterator<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> >, Standard>::Type TIter_;
    typedef typename GetValue<TIter_>::Type                                                                  Type;
};

template <typename... TContTypes, typename TSpec>
struct GetValue<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> > const>
{
    typedef typename Iterator<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> > const, Standard>::Type TIter_;
    typedef typename GetValue<TIter_>::Type                                                                        Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename... TContTypes, typename TSpec>
struct Reference<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> > >
{
    typedef typename Iterator<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> >, Standard>::Type TIter_;
    typedef typename Reference<TIter_>::Type                                                                 Type;
};

template <typename... TContTypes, typename TSpec>
struct Reference<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> > const>
{
    typedef typename Iterator<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> > const, Standard>::Type TIter_;
    typedef typename Reference<TIter_>::Type                                                                       Type;
};

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <typename... TContTypes, typename TSpec>
struct Difference<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> > >
{
    typedef typename Iterator<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> >, Standard>::Type TIter_;
    typedef typename Difference<TIter_>::Type                                                                Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename... TContTypes, typename TSpec>
struct Size<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> > >
{
    typedef typename Iterator<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> >, Standard>::Type TIter_;
    typedef typename Size<TIter_>::Type                                                                      Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename... TContTypes, typename TSpec>
struct Position<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> > >
{
    typedef typename Iterator<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> >, Standard>::Type TIter_;
    typedef typename Position<TIter_>::Type                                                                  Type;
};

// TODO(rrahn): IsSequence
// TODO(rrahn): Infix
// TODO(rrahn): Prefix
// TODO(rrahn): Suffix
// TODO(rrahn): View (Segment)

// Create Zipped Iterator.
template <typename... TContTypes, typename TSpec>
struct Iterator<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> >, Standard>
{
    typedef Iter<std::tuple<typename Iterator<TContTypes, Standard>::Type...>, ZipIterator> Type;
};

template <typename... TContTypes, typename TSpec>
struct Iterator<ContainerView<std::tuple<TContTypes...>, ZipContainer<TSpec> > const, Standard>
{
    typedef Iter<std::tuple<typename Iterator<TContTypes const, Standard>::Type...>, ZipIterator> Type;
};

// ============================================================================
// Functions
// ============================================================================

// Uses default functions from base ContainerView.

// ----------------------------------------------------------------------------
// Function makeZipView()
// ----------------------------------------------------------------------------

/*!
 * @fn makeZipView
 * @headerfile <seqan/sequence.h>
 * @brief Creates a @link ZipContainerView @endlink, deducing the container types from the arguments.
 *
 * @signature cont makeZipView(TContainerTypes... args)
 *
 * @param [in] args One or more container instance to construct the @link ZipContainerView @endlink from.
 *
 * @return cont A @link ZipContainerView @endlink containing the views over the given containers..
 */

// Helper function to create zipped view from container pack.

template <typename... TContTypes>
inline ContainerView<std::tuple<typename std::remove_reference<TContTypes>::type...>, ZipContainer<> >
makeZipView(TContTypes && ...contArgs)
{
#ifdef SEQAN_CLANG35_FREEBSD_BUG
    // the condition always evaluates to false, but ensures that the assertion
    // only fires if the function is actually instantiated
    static_assert(sizeof...(contArgs) == 0, "The Zip Container triggers a bug on FreeBSD+clang-3.5, please upgrade you compiler!");
#endif
    return ContainerView<std::tuple<typename std::remove_reference<TContTypes>::type...>, ZipContainer<> >(std::forward<TContTypes>(contArgs)...);
}

}

#endif  // #ifndef INCLUDE_SEQAN_SEQUENCE_CONTAINER_VIEW_ZIP_H_
