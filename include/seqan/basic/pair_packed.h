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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Packed pair specialization.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_PAIR_PACKED_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_PAIR_PACKED_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Specialization Packed Pair
// ----------------------------------------------------------------------------

/*!
 * @class PackedPair
 * @extends Pair
 * @headerfile <seqan/basic.h>
 * @brief Stores two arbitrary objects. Saves memory by disabling memory alignment.
 *
 * @signature template <typename T1, typename T2>
 *            class Pair<T1, T2, Pack>;
 *
 * @tparam T1 The type of the first object.
 * @tparam T2 The type of the second object.
 *
 * Useful for external storage.
 *
 * Memory access could be slower.  Direct access to members by pointers is not allowed on all platforms.
 *
 * Functions <tt>value()</tt> is not implemented yet since there it would require using a proxy.  Use
 * <tt>getValue()</tt>, <tt>assignValue()</tt>, <tt>moveValue()</tt>, <tt>setValue()</tt> instead.
 */

#pragma pack(push,1)
template <typename T1, typename T2>
struct Pair<T1, T2, Pack>
{
    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------
    T1 i1{};
    T2 i2{};

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    // Pair() = default; does not work on gcc4.9, it issues warnings if T1/T2
    // have no proper default constructor. >=gcc5.0 reports no warnings.
    // Caused by yara_indexer build, demo_tutorial_indices_base and
    // demo_tutorial_index_iterators_index_bidirectional_search.
#if defined(COMPILER_GCC) && (__GNUC__ <= 4)
    Pair() : i1(T1()), i2(T2()) {};
#else
    Pair() = default;
#endif

    // NOTE(marehr) intel compiler bug in 17.x and 18.x: defaulted copy-constructor
    // in classes with `#pragma pack(push, 1)` seg-faults. This leads to a
    // seg-fault in yara-mapper (app test case yara).
#if defined(COMPILER_LINTEL) || defined(COMPILER_WINTEL)
    Pair(Pair const & p) : i1(p.i1), i2(p.i2) {};
#else
    Pair(Pair const &) = default;
#endif
    Pair(Pair &&) = default;
    ~Pair() = default;
    Pair & operator=(Pair const &) = default;
    Pair & operator=(Pair &&) = default;

    Pair(T1 const & _i1, T2 const & _i2) : i1(_i1), i2(_i2) {}

    template <typename T1_, typename T2_, typename TSpec__>
    // TODO(holtgrew): explicit?
    Pair(Pair<T1_, T2_, TSpec__> const &_p) :
        i1(getValueI1(_p)), i2(getValueI2(_p))
    {}
};
#pragma pack(pop)

// ============================================================================
// Metafunctions
// ============================================================================

template <typename T1, typename T2>
struct MakePacked< Pair<T1, T2> >
{
    typedef Pair<T1, T2, Pack> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function set().
// ----------------------------------------------------------------------------

// References to members to packed structs do not work.  Always copy.

template <typename T1, typename T2>
inline void
set(Pair<T1, T2, Pack> & p1, Pair<T1, T2, Pack> & p2)
{
    p1 = p2;
}

// ----------------------------------------------------------------------------
// Function move().
// ----------------------------------------------------------------------------

// References to members to packed structs do not work.  Always copy.

template <typename T1, typename T2>
inline void
move(Pair<T1, T2, Pack> & p1, Pair<T1, T2, Pack> & p2)
{
    p1 = p2;
}

// -----------------------------------------------------------------------
// Function setValueIX()
// -----------------------------------------------------------------------

// References to members to packed structs do not work.  Always copy.

template <typename T1, typename T2, typename T>
inline void setValueI1(Pair<T1, T2, Pack> & pair, T const & _i)
{
    pair.i1 = _i;
}

template <typename T1, typename T2, typename T>
inline void setValueI2(Pair<T1, T2, Pack> & pair, T const & _i)
{
    pair.i2 = _i;
}

// -----------------------------------------------------------------------
// Function moveValueIX()
// -----------------------------------------------------------------------

// References to members to packed structs do not work.  Always copy.

template <typename T1, typename T2, typename T>
inline void moveValueI1(Pair<T1, T2, Pack> & pair, T & _i)
{
    pair.i1 = _i;
}

template <typename T1, typename T2, typename T>
inline void moveValueI2(Pair<T1, T2, Pack> & pair, T & _i)
{
    pair.i2 = _i;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_PAIR_PACKED_H_
