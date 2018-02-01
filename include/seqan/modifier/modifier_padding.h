// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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

#ifndef INCLUDE_SEQAN_MODIFIER_MODIFIER_PADDING_H_
#define INCLUDE_SEQAN_MODIFIER_MODIFIER_PADDING_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================
// --------------------------------------------------------------------------
// Class ModPad Iterator
// --------------------------------------------------------------------------

struct ModPadding_;
using ModPadding = Tag<ModPadding_>;

/*!
 * @class ModPaddingIterator
 * @extends ModifiedIterator
 * @headerfile <seqan/modifier.h>
 * @brief Adds padding characters beyond the end of a string, without modifing the original string.
 *
 * @signature template <typename THost>
 *            class ModifiedIterator<THost, ModPadding>;
 *
 * @tparam THost original iterator.
 */

/*!
 * @class ModPaddingString
 * @extends ModifiedString
 * @headerfile <seqan/modifier.h>
 * @brief Pad characters beyond the end of a string with default value.
 *
 * @signature template <typename THost>
 *            class ModifiedString<THost, ModPadding>;
 *
 * @tparam THost original string.
 */

template <typename THost>
struct ModPaddingCargo
{
    using TSize  = typename Size<THost>::Type;
    using TValue = typename std::remove_const<typename Value<THost>::Type>::type;

    TSize   _numPaddedChar  = 0;
    TSize   _remainingSteps = 0;
    TValue  _paddedValue    = TValue();
};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction Reference
// --------------------------------------------------------------------------

template <typename THost>
struct Reference<ModifiedString<THost, ModPadding> > :
    Reference<THost>
{};

template < typename THost>
struct Reference<ModifiedString<THost, ModPadding> const> :
    Reference<THost const>
{};

// --------------------------------------------------------------------------
// Metafunction Cargo                           [ModPadding ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
struct Cargo<ModifiedString<THost, ModPadding> >
{
    using Type = ModPaddingCargo<THost>;
};

template <typename THost>
struct Cargo<ModifiedIterator<THost, ModPadding> >
{
    using Type = ModPaddingCargo<typename Container<THost>::Type>;
};

// --------------------------------------------------------------------------
// Metafunction Iterator                          [ModPadding ModifiedString]
// --------------------------------------------------------------------------

template <typename THost>
struct Iterator<ModifiedString<THost, ModPadding>, Standard>
{
    typedef ModifiedIterator<typename Iterator<THost, Rooted>::Type, ModPadding> Type;
};

template <typename THost>
struct Iterator<ModifiedString<THost, ModPadding> const, Standard>
{
    typedef ModifiedIterator<typename Iterator<THost, Rooted>::Type, ModPadding> Type;
};

// --------------------------------------------------------------------------
// Metafunction DefaultIteratorSpec               [ModReverse ModifiedString]
// --------------------------------------------------------------------------

template <typename THost>
struct DefaultIteratorSpec< ModifiedString<THost, ModPadding> >
{
    typedef Rooted Type;
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function expand()
// --------------------------------------------------------------------------

/*!
 * @fn ModPaddingString#expand
 * @headerfile <seqan/modifier.h>
 * @brief Expands the original string by the given size.
 *
 * @signature void expand(str, size[, pad])
 *
 * @param [in,out] str  The modified string to be padded.
 * @param [in]     size The number of padded characters.
 * @param [in]     pad  The character to pad the sequence with.
 *
 * @datarace Not thread-safe.
 */

template <typename T>
bool _isValid(T* value)
{
    return value != nullptr;
}

template <typename T>
bool _isValid(T const &)
{
    return true;
}

template <typename THost, typename TSize, typename TPadding>
inline void expand(ModifiedString<THost, ModPadding> & me,
                   TSize const newSize,
                   TPadding const & _padding)
{
    SEQAN_ASSERT(_isValid(me._host));

    cargo(me)._numPaddedChar = newSize;
    cargo(me)._paddedValue = _padding;
}

template <typename THost, typename TSize>
inline void expand(ModifiedString<THost, ModPadding> & me,
                   TSize const newSize)
{
    expand(me, newSize, typename Value<THost>::Type());
}

// --------------------------------------------------------------------------
// Function length()
// --------------------------------------------------------------------------

template <typename THost>
inline auto
length(ModifiedString<THost, ModPadding> const & me)
{
    return length(host(me)) + cargo(me)._numPaddedChar;
}

// ----------------------------------------------------------------------------
// Function cargoValue()
// ----------------------------------------------------------------------------

template <typename THost>
inline typename Reference<ModifiedString<THost, ModPadding> >::Type
cargoValue(ModifiedString<THost, ModPadding> & me)
{
    return cargo(me)._paddedValue;
}

// NOTE(rrahn): The problem with the padding symbol is, that it is always stored as a member
// of the modifier class. Hence, if the modifier is const all it's members are const.
// Now, the cargo could be either defined mutable or, and this what we did right now, the
// the const is cast-away. However, we use SFINAE to only apply this hack to the Host types,
// for which this becomes relevant. That are Host types like the Segment class who copy pointer semantics, i.e.
// the constness of the object is not propagated to the underlying source.

// The default version, where Reference<THost const>::Type gives back a const reference.
template <typename THost,
          std::enable_if_t<std::is_same<std::remove_reference_t<
                                            typename Reference<ModifiedString<THost, ModPadding>>::Type>,
                                        std::add_const_t<std::remove_reference_t<
                                            typename Reference<ModifiedString<THost, ModPadding>>::Type>>>::value,
                           int> = 0>
inline typename Reference<ModifiedString<THost, ModPadding> >::Type
cargoValue(ModifiedString<THost, ModPadding> const & me)
{
    return cargo(me)._paddedValue;
}

// The version, where Reference<THost const>::Type gives back a non-const reference.
template <typename THost,
          std::enable_if_t<!std::is_same<std::remove_reference_t<
                                            typename Reference<ModifiedString<THost, ModPadding>>::Type>,
                                         std::add_const_t<std::remove_reference_t<
                                            typename Reference<ModifiedString<THost, ModPadding>>::Type>>>::value,
                           int> = 0>
inline typename Reference<ModifiedString<THost, ModPadding> >::Type
cargoValue(ModifiedString<THost, ModPadding> const & me)
{
    using TTargetType = typename Reference<ModifiedString<THost, ModPadding> >::Type;
    return const_cast<TTargetType>(cargo(me)._paddedValue);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename THost, typename TPosition>
inline typename Reference<ModifiedString<THost, ModPadding> >::Type
value(ModifiedString<THost, ModPadding> & me, TPosition const pos)
{
    SEQAN_ASSERT_LT(pos, static_cast<TPosition>(length(me)));
    return (SEQAN_LIKELY(pos < static_cast<TPosition>(length(host(me))))) ? host(me)[pos] : cargoValue(me);
}

template <typename THost, typename TPosition>
inline typename Reference<ModifiedString<THost, ModPadding> const>::Type
value(ModifiedString<THost, ModPadding> const & me, TPosition const pos)
{
    SEQAN_ASSERT_LT(pos, static_cast<TPosition>(length(me)));
        return (SEQAN_LIKELY(pos < static_cast<TPosition>(length(host(me))))) ? value(host(me), pos) : cargoValue(me);
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

template <typename THost, typename TPosition>
inline typename GetValue<THost>::Type
getValue(ModifiedString<THost, ModPadding> const & me, TPosition const pos)
{
    SEQAN_ASSERT_LT(pos, static_cast<TPosition>(length(me)));
    return (SEQAN_LIKELY(pos < static_cast<TPosition>(length(host(me))))) ? host(me)[pos] : cargo(me)._paddedValue;
}

// --------------------------------------------------------------------------
// Function begin()                               [ModReverse ModifiedString]
// --------------------------------------------------------------------------

template < typename THost, typename TTagSpec>
inline typename Iterator< ModifiedString<THost, ModPadding> const>::Type
begin(ModifiedString<THost, ModPadding> const & me, Tag<TTagSpec> const & /*tag*/)
{
    typename Iterator<ModifiedString<THost, ModPadding> const, Standard>::Type temp(begin(host(me), Rooted()));

    _copyCargo(temp, me);
    cargo(temp)._remainingSteps = cargo(me)._numPaddedChar;
    return temp;
}

template < typename THost, typename TTagSpec>
inline typename Iterator< ModifiedString<THost, ModPadding> >::Type
begin(ModifiedString<THost, ModPadding> & me, Tag<TTagSpec> const & /*tag*/)
{
    typename Iterator<ModifiedString<THost, ModPadding>, Standard>::Type temp(begin(host(me), Rooted()));

    _copyCargo(temp, me);
    cargo(temp)._remainingSteps = cargo(me)._numPaddedChar;
    return temp;
}

// --------------------------------------------------------------------------
// Function end()                                 [ModReverse ModifiedString]
// --------------------------------------------------------------------------

template <typename THost, typename TTagSpec>
inline auto
end(ModifiedString<THost, ModPadding> const & me, Tag<TTagSpec> const)
{
    typename Iterator<ModifiedString<THost, ModPadding> const, Standard>::Type temp(end(host(me), Rooted()));

    _copyCargo(temp, me);
    cargo(temp)._remainingSteps = 0;
    return temp;
}

template <typename THost, typename TTagSpec>
inline auto
end(ModifiedString<THost, ModPadding> & me, Tag<TTagSpec> const)
{
    typename Iterator<ModifiedString<THost, ModPadding>, Standard>::Type temp(end(host(me), Rooted()));

    _copyCargo(temp, me);
    cargo(temp)._remainingSteps = 0;
    return temp;
}

// ----------------------------------------------------------------------------
// Function operator*()
// ----------------------------------------------------------------------------

template <typename THost>
inline auto
operator*(ModifiedIterator<THost, ModPadding> & me)
{
    return (SEQAN_UNLIKELY(atEnd(host(me)))) ? cargo(me)._paddedValue : *host(me);
}

template <typename THost>
inline auto
operator*(ModifiedIterator<THost, ModPadding> const & me)
{
    return (SEQAN_UNLIKELY(atEnd(host(me)))) ? (cargo(me)._paddedValue) : (*host(me));
}

// ----------------------------------------------------------------------------
// Function operator++()
// ----------------------------------------------------------------------------

template <typename THost>
inline ModifiedIterator<THost, ModPadding> &
operator++(ModifiedIterator<THost, ModPadding> & me)
{
    if (SEQAN_UNLIKELY(atEnd(host(me))))
        --cargo(me)._remainingSteps;
    else
        ++host(me);
    return me;
}

// ----------------------------------------------------------------------------
// Function operator+=()
// ----------------------------------------------------------------------------

template <typename THost, typename TSize>
inline ModifiedIterator<THost, ModPadding> &
operator+=(ModifiedIterator<THost, ModPadding> & me, TSize const steps)
{
    if (SEQAN_UNLIKELY(atEnd(host(me))))
    {
        cargo(me)._remainingSteps -= steps;   // Remove 'steps' from remaining steps.
    }
    else
    {
        auto rem = (end(container(host(me)), Rooted()) - host(me));
        if (SEQAN_LIKELY(static_cast<decltype(rem)>(steps) <= rem))  // Move host by 'steps' forward.
        {
            host(me) += steps;
        }
        else  // Move host by 'rem' forward and remove diff from cargo.
        {
            host(me) += rem;
            cargo(me)._remainingSteps -= (steps - rem);
        }
    }
    return me;
}

// ----------------------------------------------------------------------------
// Function operator--()
// ----------------------------------------------------------------------------

template <typename THost>
inline ModifiedIterator<THost, ModPadding> &
operator--(ModifiedIterator<THost, ModPadding> & me)
{
    if (SEQAN_LIKELY(cargo(me)._remainingSteps == cargo(me)._numPaddedChar))
        --host(me);
    else
        ++cargo(me)._remainingSteps;
    return me;
}

// ----------------------------------------------------------------------------
// Function operator-=()
// ----------------------------------------------------------------------------

template <typename THost, typename TSize>
inline ModifiedIterator<THost, ModPadding> &
operator-=(ModifiedIterator<THost, ModPadding> & me, TSize const steps)
{
    if (SEQAN_UNLIKELY(atEnd(host(me))))
    {
        auto rem = cargo(me)._numPaddedChar - cargo(me)._remainingSteps;
        if (static_cast<decltype(rem)>(steps) <= rem)
        {
            cargo(me)._remainingSteps += steps;
        }
        else
        {
            cargo(me)._remainingSteps = cargo(me)._numPaddedChar;
            host(me) -= steps - rem;
        }
    }
    else
    {
        host(me) -= steps;
    }
    return me;
}

// ----------------------------------------------------------------------------
// Function operator-()
// ----------------------------------------------------------------------------

template <typename THost>
inline typename Difference<ModifiedIterator<THost, ModPadding> >::Type
operator-(ModifiedIterator<THost, ModPadding> const & a,
          ModifiedIterator<THost, ModPadding> const & b)
{
    return host(a) - host(b) + cargo(b)._remainingSteps - cargo(a)._remainingSteps;
}

// --------------------------------------------------------------------------
// Function operator==()
// --------------------------------------------------------------------------

template <typename THost>
inline bool
operator == (ModifiedIterator<THost, ModPadding> const & a, ModifiedIterator<THost, ModPadding> const & b)
{
    return host(a) == host(b) && cargo(a)._remainingSteps == cargo(b)._remainingSteps;
}

// --------------------------------------------------------------------------
// Function operator!=()
// --------------------------------------------------------------------------

template <typename THost>
inline bool
operator != (ModifiedIterator<THost, ModPadding> const & a, ModifiedIterator<THost, ModPadding> const & b)
{
    return !(a == b);
}

}

#endif  // #ifndef INCLUDE_SEQAN_MODIFIER_MODIFIER_PADDING_H_
