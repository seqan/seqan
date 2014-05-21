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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements the coverage store to hold the data of the coverage.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_COVERAGE_STORE_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_COVERAGE_STORE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DeltaCoverageStore
// ----------------------------------------------------------------------------

template <typename TSpec = void>
class DeltaCoverageStore
{
public:
    typedef typename Value<DeltaCoverageStore>::Type TBitVector;
    typedef typename Size<TBitVector>::Type TSize;

    TSize              _coverageSize;
    String<TBitVector> _coverageData;

    DeltaCoverageStore() : _coverageSize(0)
    {}

    template <typename TSize>
    DeltaCoverageStore(TSize const & newSize) : _coverageSize(newSize)
    {}

    template <typename TPosition>
    inline TBitVector &
    operator[](TPosition const & pos)
    {
        return value(*this, pos);
    }

    template <typename TPosition>
    inline TBitVector const &
    operator[](TPosition const & pos) const
    {
        return value(*this, pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Value<DeltaCoverageStore<TSpec> >
{
    typedef String<bool, Packed<> > Type;
};

template <typename TSpec>
struct Value<DeltaCoverageStore<TSpec> const>
{
    typedef String<bool, Packed<> > const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Reference<DeltaCoverageStore<TSpec> >
{
    typedef String<bool, Packed<> > & Type;
};

template <typename TSpec>
struct Reference<DeltaCoverageStore<TSpec> const>
{
    typedef String<bool, Packed<> > const & Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Size<DeltaCoverageStore<TSpec> >
{
    typedef typename Value<DeltaCoverageStore<TSpec> >::Type TValue_;
    typedef typename Size<TValue_>::Type Type;
};

template <typename TSpec>
struct Size<DeltaCoverageStore<TSpec> const> :
    Size<DeltaCoverageStore<TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Position<DeltaCoverageStore<TSpec> > :
    Size<DeltaCoverageStore<TSpec> >{};

template <typename TSpec>
struct Position<DeltaCoverageStore<TSpec> const> :
    Size<DeltaCoverageStore<TSpec> const>{};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Iterator<DeltaCoverageStore<TSpec>, Standard>
{
    typedef typename Value<DeltaCoverageStore<TSpec> >::Type TBitVector;
    typedef String<TBitVector> TCoverage;
    typedef typename Iterator<TCoverage, Standard>::Type Type;
};

template <typename TSpec>
struct Iterator<DeltaCoverageStore<TSpec> const, Standard>
{
    typedef typename Value<DeltaCoverageStore<TSpec> >::Type TBitVector;
    typedef String<TBitVector> TCoverage;
    typedef typename Iterator<TCoverage const, Standard>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline typename Iterator<DeltaCoverageStore<TSpec>, Standard>::Type
begin(DeltaCoverageStore<TSpec> & store, Standard const & /*tag*/)
{
    return begin(store._coverageData, Standard());
}
template <typename TSpec>
inline typename Iterator<DeltaCoverageStore<TSpec> const, Standard>::Type
begin(DeltaCoverageStore<TSpec> const & store, Standard const & /*tag*/)
{
    return begin(store._coverageData, Standard());
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline typename Iterator<DeltaCoverageStore<TSpec>, Standard>::Type
end(DeltaCoverageStore<TSpec> & store, Standard const & /*tag*/)
{
    return end(store._coverageData, Standard());
}

template <typename TSpec>
inline typename Iterator<DeltaCoverageStore<TSpec> const, Standard>::Type
end(DeltaCoverageStore<TSpec> const & store, Standard const & /*tag*/)
{
    return end(store._coverageData, Standard());
}

// ----------------------------------------------------------------------------
// Function setCoverageSize()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TSize>
inline void
setCoverageSize(DeltaCoverageStore<TSpec> & store, TSize newSize)
{
    typedef typename Iterator<DeltaCoverageStore<TSpec>, Standard>::Type TIter;

    store._coverageSize = newSize;

    // Update all coverages to the new size.
    if (empty(store._coverageData))
        return;

    for (TIter it = begin(store, Standard()); it != end(store, Standard()); ++it)
        if (length(*it) != store._coverageSize)
            resize(*it, store._coverageSize, false, Exact());
}

// ----------------------------------------------------------------------------
// Function coverageSize()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline typename Size<DeltaCoverageStore<TSpec> >::Type
coverageSize(DeltaCoverageStore<TSpec> const & store)
{
    return store._coverageSize;
}

// ----------------------------------------------------------------------------
// Function addCoverage()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPosition>
inline void
addCoverage(DeltaCoverageStore<TSpec> & store,
            typename Value<DeltaCoverageStore<TSpec> >::Type const & coverage,
            TPosition const & pos)
{
    if (length(coverage) != coverageSize(store))
        setCoverageSize(store, _max(length(coverage), coverageSize(store)));

    insertValue(store._coverageData, pos, coverage);
}

template <typename TSpec>
inline void
addCoverage(DeltaCoverageStore<TSpec> & store,
            typename Value<DeltaCoverageStore<TSpec> >::Type const & coverage)
{
    addCoverage(store, coverage, length(store._coverageData));
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPosition>
inline typename Reference<DeltaCoverageStore<TSpec> >::Type
value(DeltaCoverageStore<TSpec> & store, TPosition const & pos)
{
    return value(store._coverageData, pos);
}

template <typename TSpec, typename TPosition>
inline typename Reference<DeltaCoverageStore<TSpec> const>::Type
value(DeltaCoverageStore<TSpec> const & store, TPosition const & pos)
{
    return value(store._coverageData, pos);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
clear(DeltaCoverageStore<TSpec> & store)
{
    return clear(store._coverageData);
    setCoverageSize(store, 0);
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TSize, typename TExpand>
inline void
resize(DeltaCoverageStore<TSpec> & store, TSize const & newSize, Tag<TExpand> const & tag)
{
    resize(store._coverageData, newSize, tag);
}

template <typename TSpec, typename TSize>
inline void
resize(DeltaCoverageStore<TSpec> const & store, TSize const & newSize)
{
    typedef typename Value<DeltaCoverageStore<TSpec> >::Type TValue;
    typedef String<TValue> TData;

    resize(store, newSize, typename DefaultOverflowExplicit<TData>::Type());
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline typename Size<DeltaCoverageStore<TSpec> >::Type
length(DeltaCoverageStore<TSpec> const & store)
{
    return length(store._coverageData);
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_COVERAGE_STORE_H_
