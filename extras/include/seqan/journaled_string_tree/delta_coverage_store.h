// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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

struct DeltaCoverageString_;
typedef Tag<DeltaCoverageString_> DeltaCoverageString;

// ----------------------------------------------------------------------------
// Class DeltaCoverageStore
// ----------------------------------------------------------------------------

template <typename TSpec = void>
class DeltaCoverageStore
{
public:
    typedef typename Member<DeltaCoverageStore, DeltaCoverageString>::Type TCoverageString;
    typedef typename Value<DeltaCoverageStore>::Type TBitVector;
    typedef typename Size<DeltaCoverageStore>::Type TSize;

    TSize           _coverageSize;
    TCoverageString _coverageData;

    DeltaCoverageStore() : _coverageSize(0)
    {}

    template <typename TSize>
    DeltaCoverageStore(TSize newSize) : _coverageSize(newSize)
    {}

    template <typename TPosition>
    inline TBitVector &
    operator[](TPosition pos)
    {
        return value(*this, pos);
    }

    template <typename TPosition>
    inline TBitVector const &
    operator[](TPosition pos) const
    {
        return value(*this, pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Member
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Member<DeltaCoverageStore<TSpec>, DeltaCoverageString>
{
    typedef typename Value<DeltaCoverageStore<TSpec> >::Type TValue_;
    typedef String<TValue_> Type;
};

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
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Size<DeltaCoverageStore<TSpec> >
{
    typedef typename Member<DeltaCoverageStore<TSpec>, DeltaCoverageString>::Type TCoverageString_;
    typedef typename Size<TCoverageString_>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Iterator<DeltaCoverageStore<TSpec>, Standard>
{
    typedef typename Member<DeltaCoverageStore<TSpec>, DeltaCoverageString>::Type TCoverageString_;
    typedef typename Iterator<TCoverageString_, Standard>::Type Type;
};

template <typename TSpec>
struct Iterator<DeltaCoverageStore<TSpec> const, Standard>
{
    typedef typename Member<DeltaCoverageStore<TSpec>, DeltaCoverageString>::Type const TCoverageString_;
    typedef typename Iterator<TCoverageString_, Standard>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _syncCoverageSizes()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void _syncCoverageSize(DeltaCoverageStore<TSpec> & store)
{
    typedef typename Iterator<DeltaCoverageStore<TSpec>, Standard>::Type TIter;

    for (TIter it = begin(store, Standard()); it != end(store, Standard()); ++it)
        resize(*it, store._coverageSize, false, Exact());
}

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
    _syncCoverageSize(store);
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
// Function insert()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPosition>
inline void
insert(DeltaCoverageStore<TSpec> & store,
       TPosition pos,
       typename Value<DeltaCoverageStore<TSpec> >::Type const & coverage)
{
    // Adapt to coverage size of store.
    if (length(coverage) < coverageSize(store))
    {
        typename Value<DeltaCoverageStore<TSpec> >::Type tmp(coverage);
        resize(tmp, coverageSize(store), false, Exact());
        insertValue(store._coverageData, pos, tmp);
    }
    // Set new coverage size if necessary.
    if (length(coverage) > coverageSize(store))
        setCoverageSize(store, length(coverage));
    insertValue(store._coverageData, pos, coverage);
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
    typedef typename Value<DeltaCoverageStore<TSpec> >::Type TCoverage;
    TCoverage tmp;
    resize(tmp, store._coverageSize, false, Exact());
    resize(store._coverageData, newSize, tmp, tag);
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
