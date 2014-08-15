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
// Implements a data structure to store deltas efficiently.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_STORE_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_STORE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

template <typename TDeltaStore, typename TDeltaType>
struct DeltaValue;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Struct DeltaType
// ----------------------------------------------------------------------------

/*!
 * @enum DeltaType
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Keys for specifying the delta type to be accessed.
 *
 * @val DeltaType DELTA_TYPE_SNP
 * @brief Id to denote SNP events.
 *
 * @val DeltaType DELTA_TYPE_DEL
 * @brief Id to denote deletion events.
 *
 * @val DeltaType DELTA_TYPE_INS
 * @brief Id to denote insertion events.
 *
 * @val DeltaType DElTA_TYPE_INDEL
 * @brief Id to denote replacement events.
 */

enum DeltaType
{
    DELTA_TYPE_SNP = 0,
    DELTA_TYPE_DEL = 1,
    DELTA_TYPE_INS = 2,
    DELTA_TYPE_SV = 3
};

struct DeltaTypeSnp_;
typedef Tag<DeltaTypeSnp_> DeltaTypeSnp;

struct DeltaTypeDel_;
typedef Tag<DeltaTypeDel_> DeltaTypeDel;

struct DeltaTypeIns_;
typedef Tag<DeltaTypeIns_> DeltaTypeIns;

struct DeltaTypeSV_;
typedef Tag<DeltaTypeSV_> DeltaTypeSV;

typedef TagList<DeltaTypeSnp,
        TagList<DeltaTypeDel,
        TagList<DeltaTypeIns,
        TagList<DeltaTypeSV> > > > DeltaTypes;

typedef TagSelector<DeltaTypes> DeltaTypeSelector;

// ----------------------------------------------------------------------------
// Class DeltaStore
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
class DeltaStore
{
public:
    typedef typename Member<DeltaStore, DeltaTypeSnp>::Type TSnpData;
    typedef typename Member<DeltaStore, DeltaTypeDel>::Type TDelData;
    typedef typename Member<DeltaStore, DeltaTypeIns>::Type TInsData;
    typedef typename Member<DeltaStore, DeltaTypeSV>::Type  TSVData;

    // TODO(rmaerker): Elaborate on these ideas!
    // Idea a) Use ConcatStringSet for insertions. Use as global insertion buffer for all journal sequences.
    // Idea b) Use bit encoding for DNA alphabet.
    // Idea c) Instead of insertion buffer, we append the inserted strings to the reference and only use original nodes.

    TSnpData    _snpData;
    TDelData    _delData;
    TInsData    _insData;
    TSVData     _svData;

    DeltaStore()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction BitsPerValue
// ----------------------------------------------------------------------------

template <>
struct BitsPerValue<DeltaType>
{
    static const unsigned VALUE = 2;
};

// ----------------------------------------------------------------------------
// Metafunction Member
// ----------------------------------------------------------------------------

// DeltaTypeSnp.
template <typename TSize, typename TAlphabet>
struct Member<DeltaStore<TSize, TAlphabet>, DeltaTypeSnp>
{
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaTypeSnp>::Type TSnpValue;
    typedef String<TSnpValue> Type;
};

// DeltaTypeDel.
template <typename TSize, typename TAlphabet>
struct Member<DeltaStore<TSize, TAlphabet>, DeltaTypeDel>
{
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaTypeDel>::Type TDelValue;
    typedef String<TDelValue> Type;
};

// DeltaTypeIns.
template <typename TSize, typename TAlphabet>
struct Member<DeltaStore<TSize, TAlphabet>, DeltaTypeIns>
{
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaTypeIns>::Type TInsValue;
    typedef String<TInsValue> Type;
};

// DeltaTypeSV.
template <typename TSize, typename TAlphabet>
struct Member<DeltaStore<TSize, TAlphabet>, DeltaTypeSV>
{
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaTypeSV>::Type TSVValue;
    typedef String<TSVValue> Type;
};

// Generic const version.
template <typename TSize, typename TAlphabet, typename TDeltaType>
struct Member<DeltaStore<TSize, TAlphabet> const, TDeltaType>
{
    typedef DeltaStore<TSize, TAlphabet> TStore_;
    typedef typename Member<TStore_, TDeltaType>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                     [DELTA_TYPE_SNP]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaTypeSnp>
{
    typedef TAlphabet Type;
};

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaTypeSnp>
{
    typedef TAlphabet const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                     [DELTA_TYPE_DEL]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaTypeDel>
{
    typedef TSize Type;
};

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaTypeDel>
{
    typedef TSize const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                     [DELTA_TYPE_INS]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaTypeIns>
{
    typedef String<TAlphabet> Type;
};

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaTypeIns>
{
    typedef String<TAlphabet> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                   [DELTA_TYPE_SV]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaTypeSV>
{
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaTypeIns>::Type TIns_;
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaTypeDel>::Type TDel_;
    typedef Pair<TDel_, TIns_> Type;
};

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaTypeSV>
{
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaTypeIns>::Type TIns_;
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaTypeDel>::Type TDel_;
    typedef Pair<TDel_, TIns_> const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function isDeltaType()                                        [DeltaTypeSnp]
// ----------------------------------------------------------------------------

inline bool
isDeltaType(DeltaType dType, DeltaTypeSnp const & /*tag*/)
{
    return dType == DELTA_TYPE_SNP;
}

// ----------------------------------------------------------------------------
// Function isDeltaType()                                        [DeltaTypeDel]
// ----------------------------------------------------------------------------

inline bool
isDeltaType(DeltaType dType, DeltaTypeIns const & /*tag*/)
{
    return dType == DELTA_TYPE_INS;
}

// ----------------------------------------------------------------------------
// Function isDeltaType()                                        [DeltaTypeIns]
// ----------------------------------------------------------------------------

inline bool
isDeltaType(DeltaType dType, DeltaTypeDel const & /*tag*/)
{
    return dType == DELTA_TYPE_DEL;
}

// ----------------------------------------------------------------------------
// Function isDeltaType()                                         [DeltaTypeSV]
// ----------------------------------------------------------------------------

inline bool
isDeltaType(DeltaType dType, DeltaTypeSV const & /*tag*/)
{
    return dType == DELTA_TYPE_SV;
}

// ----------------------------------------------------------------------------
// Function isDeltaType()
// ----------------------------------------------------------------------------

template <typename TTag>
inline bool
isDeltaType(DeltaType /*type*/, TTag const & /*tag*/)
{
    return false; // Unknown tag.
}

// ----------------------------------------------------------------------------
// Function selectDeltaType()                                    [DeltaTypeSnp]
// ----------------------------------------------------------------------------

inline DeltaType
selectDeltaType(DeltaTypeSnp /*tag*/)
{
    return DELTA_TYPE_SNP;
}

// ----------------------------------------------------------------------------
// Function selectDeltaType()                                    [DeltaTypeDel]
// ----------------------------------------------------------------------------

inline DeltaType
selectDeltaType(DeltaTypeDel /*tag*/)
{
    return DELTA_TYPE_DEL;
}

// ----------------------------------------------------------------------------
// Function selectDeltaType()                                    [DeltaTypeIns]
// ----------------------------------------------------------------------------

inline DeltaType
selectDeltaType(DeltaTypeIns /*tag*/)
{
    return DELTA_TYPE_INS;
}

// ----------------------------------------------------------------------------
// Function selectDeltaType()                                     [DeltaTypeSV]
// ----------------------------------------------------------------------------

inline DeltaType
selectDeltaType(DeltaTypeSV /*tag*/)
{
    return DELTA_TYPE_SV;
}

// ----------------------------------------------------------------------------
// Function getDeltaStore()                                     [DeltaTypeSnp]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
inline typename Member<DeltaStore<TSize, TAlphabet>, DeltaTypeSnp>::Type &
getDeltaStore(DeltaStore<TSize, TAlphabet> & store, DeltaTypeSnp /*tag*/)
{
    return store._snpData;
}

template <typename TSize, typename TAlphabet>
inline typename Member<DeltaStore<TSize, TAlphabet> const, DeltaTypeSnp>::Type &
getDeltaStore(DeltaStore<TSize, TAlphabet> const & store, DeltaTypeSnp /*tag*/)
{
    return store._snpData;
}

// ----------------------------------------------------------------------------
// Function getDeltaStore()                                     [DeltaTypeDel]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
inline typename Member<DeltaStore<TSize, TAlphabet>, DeltaTypeDel>::Type &
getDeltaStore(DeltaStore<TSize, TAlphabet> & store, DeltaTypeDel /*tag*/)
{
    return store._delData;
}

template <typename TSize, typename TAlphabet>
inline typename Member<DeltaStore<TSize, TAlphabet> const, DeltaTypeDel>::Type &
getDeltaStore(DeltaStore<TSize, TAlphabet> const & store, DeltaTypeDel /*tag*/)
{
    return store._delData;
}

// ----------------------------------------------------------------------------
// Function getDeltaStore()                                     [DeltaTypeIns]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
inline typename Member<DeltaStore<TSize, TAlphabet>, DeltaTypeIns>::Type &
getDeltaStore(DeltaStore<TSize, TAlphabet> & store, DeltaTypeIns /*tag*/)
{
    return store._insData;
}

template <typename TSize, typename TAlphabet>
inline typename Member<DeltaStore<TSize, TAlphabet> const, DeltaTypeIns>::Type &
getDeltaStore(DeltaStore<TSize, TAlphabet> const & store, DeltaTypeIns /*tag*/)
{
    return store._insData;
}

// ----------------------------------------------------------------------------
// Function getDeltaStore()                                      [DeltaTypeSV]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
inline typename Member<DeltaStore<TSize, TAlphabet>, DeltaTypeSV>::Type &
getDeltaStore(DeltaStore<TSize, TAlphabet> & store, DeltaTypeSV /*tag*/)
{
    return store._svData;
}

template <typename TSize, typename TAlphabet>
inline typename Member<DeltaStore<TSize, TAlphabet> const, DeltaTypeSV>::Type &
getDeltaStore(DeltaStore<TSize, TAlphabet> const & store, DeltaTypeSV /*tag*/)
{
    return store._svData;
}

// ----------------------------------------------------------------------------
// Function addDeltaValue()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TTag>
inline typename Size<DeltaStore<TSize, TAlphabet> >::Type
addDeltaValue(DeltaStore<TSize, TAlphabet> & store,
              typename DeltaValue<DeltaStore<TSize, TAlphabet>, TTag>::Type const & value,
              TTag /*deltaType*/)
{
    appendValue(getDeltaStore(store, TTag()), value);
    return length(getDeltaStore(store, TTag())) - 1;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename Reference<DeltaStore<TSize, TAlphabet> >::Type
value(DeltaStore<TSize, TAlphabet> & deltaStore,
      TPosition const & pos)
{
    return value(deltaStore._varDataMap, pos);
}

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename Reference<DeltaStore<TSize, TAlphabet> const>::Type
value(DeltaStore<TSize, TAlphabet> const & deltaStore,
      TPosition const & pos)
{
    return value(deltaStore._varDataMap, pos);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline void
clear(DeltaStore<TSize, TAlphabet> & deltaStore)
{
    clear(deltaStore._delData);
    clear(deltaStore._svData);
    clear(deltaStore._insData);
    clear(deltaStore._snpData);
}

// ----------------------------------------------------------------------------
// Function deltaValue()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPos, typename TTag>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet>, TTag>::Type &
deltaValue(DeltaStore<TSize, TAlphabet> & store, TPos pos, TTag const & tag)
{
    return value(getDeltaStore(store, tag), pos);
}

template <typename TSize, typename TAlphabet, typename TPos, typename TTag>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet> const, TTag>::Type &
deltaValue(DeltaStore<TSize, TAlphabet> const & store, TPos pos, TTag const & tag)
{
    return value(getDeltaStore(store, tag), pos);
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_STORE_H_
