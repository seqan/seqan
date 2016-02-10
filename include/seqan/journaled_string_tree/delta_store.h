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

/*!
 * @defgroup DeltaTypeTags Delta Type Tags
 * @brief Tags used for the different delta types.
 */

/*!
 * @tag DeltaTypeTags#DeltaTypeSnp
 * @brief Tag used to select SNPs.
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @signature struct DeltaTypeSnp_;
 *            typedef Tag<DeltaTypeSnp_> DeltaTypeSnp;
 */
struct DeltaTypeSnp_;
typedef Tag<DeltaTypeSnp_> DeltaTypeSnp;

/*!
 * @tag DeltaTypeTags#DeltaTypeDel
 * @brief Tag used to select deletions.
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @signature struct DeltaTypeDel_;
 *            typedef Tag<DeltaTypeDel_> DeltaTypeDel;
 */
struct DeltaTypeDel_;
typedef Tag<DeltaTypeDel_> DeltaTypeDel;

/*!
 * @tag DeltaTypeTags#DeltaTypeIns
 * @brief Tag used to select insertions.
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @signature struct DeltaTypeIns_;
 *            typedef Tag<DeltaTypeIns_> DeltaTypeIns;
 */
struct DeltaTypeIns_;
typedef Tag<DeltaTypeIns_> DeltaTypeIns;

/*!
 * @tag DeltaTypeTags#DeltaTypeSV
 * @brief Tag used to select SVs.
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @signature struct DeltaTypeSV_;
 *            typedef Tag<DeltaTypeSV_> DeltaTypeSV;
 */
struct DeltaTypeSV_;
typedef Tag<DeltaTypeSV_> DeltaTypeSV;

typedef TagList<DeltaTypeSnp,
        TagList<DeltaTypeDel,
        TagList<DeltaTypeIns,
        TagList<DeltaTypeSV> > > > DeltaTypeTagList;

// ----------------------------------------------------------------------------
// Class DeltaTypeSelector
// ----------------------------------------------------------------------------

typedef TagSelector<DeltaTypeTagList> DeltaTypeSelector;

// ----------------------------------------------------------------------------
// Enum DeltaType
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
    DELTA_TYPE_SNP = Find<DeltaTypeSelector, DeltaTypeSnp>::VALUE,
    DELTA_TYPE_DEL = Find<DeltaTypeSelector, DeltaTypeDel>::VALUE,
    DELTA_TYPE_INS = Find<DeltaTypeSelector, DeltaTypeIns>::VALUE,
    DELTA_TYPE_SV = Find<DeltaTypeSelector, DeltaTypeSV>::VALUE
};

namespace impl
{

// ----------------------------------------------------------------------------
// Class DeltaStore
// ----------------------------------------------------------------------------

template <typename TSnp, typename TDel = uint32_t, typename TIns = String<TSnp>, typename TSV = Pair<TDel, TIns> >
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

}  // namespace impl

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
template <typename TSnp, typename TDel, typename TIns, typename TSV>
struct Member<impl::DeltaStore<TSnp, TDel, TIns, TSV>, DeltaTypeSnp>
{
    typedef String<TSnp> Type;
};

// DeltaTypeDel.
template <typename TSnp, typename TDel, typename TIns, typename TSV>
struct Member<impl::DeltaStore<TSnp, TDel, TIns, TSV>, DeltaTypeDel>
{
    typedef String<TDel> Type;
};

// DeltaTypeIns.
template <typename TSnp, typename TDel, typename TIns, typename TSV>
struct Member<impl::DeltaStore<TSnp, TDel, TIns, TSV>, DeltaTypeIns>
{
    typedef String<TIns> Type;
};

// DeltaTypeSV.
template <typename TSnp, typename TDel, typename TIns, typename TSV>
struct Member<impl::DeltaStore<TSnp, TDel, TIns, TSV>, DeltaTypeSV>
{
    typedef String<TSV> Type;
};

// Generic const version.
template <typename TSnp, typename TDel, typename TIns, typename TSV, typename TDeltaType>
struct Member<impl::DeltaStore<TSnp, TDel, TIns, TSV> const, TDeltaType>
{
    typedef impl::DeltaStore<TSnp, TDel, TIns, TSV> TStore_;
    typedef typename Member<TStore_, TDeltaType>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                     [DELTA_TYPE_SNP]
// ----------------------------------------------------------------------------

template <typename TSnp, typename TDel, typename TIns, typename TSV>
struct DeltaValue<impl::DeltaStore<TSnp, TDel, TIns, TSV>, DeltaTypeSnp>
{
    typedef TSnp Type;
};

template <typename TSnp, typename TDel, typename TIns, typename TSV>
struct DeltaValue<impl::DeltaStore<TSnp, TDel, TIns, TSV> const, DeltaTypeSnp>
{
    typedef TSnp const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                     [DELTA_TYPE_DEL]
// ----------------------------------------------------------------------------

template <typename TSnp, typename TDel, typename TIns, typename TSV>
struct DeltaValue<impl::DeltaStore<TSnp, TDel, TIns, TSV>, DeltaTypeDel>
{
    typedef TDel Type;
};

template <typename TSnp, typename TDel, typename TIns, typename TSV>
struct DeltaValue<impl::DeltaStore<TSnp, TDel, TIns, TSV> const, DeltaTypeDel>
{
    typedef TDel const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                     [DELTA_TYPE_INS]
// ----------------------------------------------------------------------------

template <typename TSnp, typename TDel, typename TIns, typename TSV>
struct DeltaValue<impl::DeltaStore<TSnp, TDel, TIns, TSV>, DeltaTypeIns>
{
    typedef TIns Type;
};

template <typename TSnp, typename TDel, typename TIns, typename TSV>
struct DeltaValue<impl::DeltaStore<TSnp, TDel, TIns, TSV> const, DeltaTypeIns>
{
    typedef TIns const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                   [DELTA_TYPE_SV]
// ----------------------------------------------------------------------------

template <typename TSnp, typename TDel, typename TIns, typename TSV>
struct DeltaValue<impl::DeltaStore<TSnp, TDel, TIns, TSV>, DeltaTypeSV>
{
    typedef TSV Type;
};

template <typename TSnp, typename TDel, typename TIns, typename TSV>
struct DeltaValue<impl::DeltaStore<TSnp, TDel, TIns, TSV> const, DeltaTypeSV>
{
    typedef TSV const Type;
};

// ============================================================================
// Punlic Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function isDeltaType()
// ----------------------------------------------------------------------------

template <typename TTag>
inline bool
isDeltaType(DeltaType type, TTag const & /*tag*/)
{
    return type == static_cast<DeltaType>(Find<DeltaTypeSelector, TTag>::VALUE);
}

// ----------------------------------------------------------------------------
// Function selectDeltaType()                                    [DeltaTypeSnp]
// ----------------------------------------------------------------------------

template <typename TTag>
constexpr inline DeltaType
selectDeltaType(TTag const & /*tag*/)
{
    return static_cast<DeltaType>(Find<DeltaTypeSelector, TTag>::VALUE);
}

// ----------------------------------------------------------------------------
// Function setDeltaType()
// ----------------------------------------------------------------------------

inline void
setDeltaType(DeltaTypeSelector & selector, DeltaType const deltaType)
{
    value(selector) = deltaType;
}

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function getDeltaStore()                                     [DeltaTypeSnp]
// ----------------------------------------------------------------------------

template <typename TSnp, typename TDel, typename TIns, typename TSV>
inline typename Member<DeltaStore<TSnp, TDel, TIns, TSV>, DeltaTypeSnp>::Type &
getDeltaStore(DeltaStore<TSnp, TDel, TIns, TSV> & store, DeltaTypeSnp /*tag*/)
{
    return store._snpData;
}

template <typename TSnp, typename TDel, typename TIns, typename TSV>
inline typename Member<DeltaStore<TSnp, TDel, TIns, TSV> const, DeltaTypeSnp>::Type &
getDeltaStore(DeltaStore<TSnp, TDel, TIns, TSV> const & store, DeltaTypeSnp /*tag*/)
{
    return store._snpData;
}

// ----------------------------------------------------------------------------
// Function getDeltaStore()                                     [DeltaTypeDel]
// ----------------------------------------------------------------------------

template <typename TSnp, typename TDel, typename TIns, typename TSV>
inline typename Member<DeltaStore<TSnp, TDel, TIns, TSV>, DeltaTypeDel>::Type &
getDeltaStore(DeltaStore<TSnp, TDel, TIns, TSV> & store, DeltaTypeDel /*tag*/)
{
    return store._delData;
}

template <typename TSnp, typename TDel, typename TIns, typename TSV>
inline typename Member<DeltaStore<TSnp, TDel, TIns, TSV> const, DeltaTypeDel>::Type &
getDeltaStore(DeltaStore<TSnp, TDel, TIns, TSV> const & store, DeltaTypeDel /*tag*/)
{
    return store._delData;
}

// ----------------------------------------------------------------------------
// Function getDeltaStore()                                     [DeltaTypeIns]
// ----------------------------------------------------------------------------

template <typename TSnp, typename TDel, typename TIns, typename TSV>
inline typename Member<DeltaStore<TSnp, TDel, TIns, TSV>, DeltaTypeIns>::Type &
getDeltaStore(DeltaStore<TSnp, TDel, TIns, TSV> & store, DeltaTypeIns /*tag*/)
{
    return store._insData;
}

template <typename TSnp, typename TDel, typename TIns, typename TSV>
inline typename Member<DeltaStore<TSnp, TDel, TIns, TSV> const, DeltaTypeIns>::Type &
getDeltaStore(DeltaStore<TSnp, TDel, TIns, TSV> const & store, DeltaTypeIns /*tag*/)
{
    return store._insData;
}

// ----------------------------------------------------------------------------
// Function getDeltaStore()                                      [DeltaTypeSV]
// ----------------------------------------------------------------------------

template <typename TSnp, typename TDel, typename TIns, typename TSV>
inline typename Member<DeltaStore<TSnp, TDel, TIns, TSV>, DeltaTypeSV>::Type &
getDeltaStore(DeltaStore<TSnp, TDel, TIns, TSV> & store, DeltaTypeSV /*tag*/)
{
    return store._svData;
}

template <typename TSnp, typename TDel, typename TIns, typename TSV>
inline typename Member<DeltaStore<TSnp, TDel, TIns, TSV> const, DeltaTypeSV>::Type &
getDeltaStore(DeltaStore<TSnp, TDel, TIns, TSV> const & store, DeltaTypeSV /*tag*/)
{
    return store._svData;
}

// ----------------------------------------------------------------------------
// Function addDeltaValue()
// ----------------------------------------------------------------------------

template <typename TSnp, typename TDel, typename TIns, typename TSV, typename TTag>
inline typename Size<DeltaStore<TSnp, TDel, TIns, TSV> >::Type
addDeltaValue(DeltaStore<TSnp, TDel, TIns, TSV> & store,
              typename DeltaValue<DeltaStore<TSnp, TDel, TIns, TSV>, TTag>::Type const & value,
              TTag /*deltaType*/)
{
    appendValue(getDeltaStore(store, TTag()), value);
    return length(getDeltaStore(store, TTag())) - 1;
}

// ----------------------------------------------------------------------------
// Function eraseDeltaValue()
// ----------------------------------------------------------------------------

template <typename TSnp, typename TDel, typename TIns, typename TSV, typename TPos, typename TTag>
inline typename Size<DeltaStore<TSnp, TDel, TIns, TSV> >::Type
eraseDeltaValue(DeltaStore<TSnp, TDel, TIns, TSV> & store,
                TPos recordPos,
                TTag /*deltaType*/)
{
    typedef typename Size<DeltaStore<TSnp, TDel, TIns, TSV> >::Type TSize;
    if (SEQAN_LIKELY(static_cast<TSize>(recordPos) < length(getDeltaStore(store, TTag()))))
        erase(getDeltaStore(store, TTag()), recordPos);
    return length(getDeltaStore(store, TTag()));
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSnp, typename TDel, typename TIns, typename TSV>
inline void
clear(DeltaStore<TSnp, TDel, TIns, TSV> & deltaStore)
{
    clear(deltaStore._delData);
    clear(deltaStore._svData);
    clear(deltaStore._insData);
    clear(deltaStore._snpData);
}

// ----------------------------------------------------------------------------
// Function deltaValue()
// ----------------------------------------------------------------------------

template <typename TSnp, typename TDel, typename TIns, typename TSV, typename TPos, typename TTag>
inline typename DeltaValue<DeltaStore<TSnp, TDel, TIns, TSV>, TTag>::Type &
deltaValue(DeltaStore<TSnp, TDel, TIns, TSV> & store, TPos pos, TTag const & tag)
{
    typedef typename Size<DeltaStore<TSnp, TDel, TIns, TSV> >::Type TSize SEQAN_TYPEDEF_FOR_DEBUG;
    SEQAN_ASSERT_LT(static_cast<TSize>(pos), length(getDeltaStore(store, tag)));

    return value(getDeltaStore(store, tag), pos);
}

template <typename TSnp, typename TDel, typename TIns, typename TSV, typename TPos, typename TTag>
inline typename DeltaValue<DeltaStore<TSnp, TDel, TIns, TSV> const, TTag>::Type &
deltaValue(DeltaStore<TSnp, TDel, TIns, TSV> const & store, TPos pos, TTag const & tag)
{
    typedef typename Size<DeltaStore<TSnp, TDel, TIns, TSV> const>::Type TSize SEQAN_TYPEDEF_FOR_DEBUG;
    SEQAN_ASSERT_LT(static_cast<TSize>(pos), length(getDeltaStore(store, tag)));

    return value(getDeltaStore(store, tag), pos);
}

// ----------------------------------------------------------------------------
// Function deletionSize()
// ----------------------------------------------------------------------------

template <typename TDeltaStore>
inline typename Size<TDeltaStore>::Type
deletionSize(typename DeltaValue<TDeltaStore, DeltaTypeSnp>::Type const & /*snp*/)
{
    return 1;
}

template <typename TDeltaStore>
inline typename Size<TDeltaStore>::Type
deletionSize(typename DeltaValue<TDeltaStore, DeltaTypeIns>::Type const & /*ins*/)
{
    return 0;
}

template <typename TDeltaStore>
inline typename Size<TDeltaStore>::Type
deletionSize(typename DeltaValue<TDeltaStore, DeltaTypeDel>::Type const & del)
{
    return del;
}

template <typename TDeltaStore>
inline typename Size<TDeltaStore>::Type
deletionSize(typename DeltaValue<TDeltaStore, DeltaTypeSV>::Type const & sv)
{
    return sv.i1;
}

template <typename TStore, typename TPos, typename TTag>
inline typename Size<TStore>::Type
deletionSize(TStore const & store, TPos const pos, TTag const & /*tag*/)
{
    return deletionSize<TStore>(getDeltaStore(store, TTag())[pos]);
}

// ----------------------------------------------------------------------------
// Function insertionSize()
// ----------------------------------------------------------------------------

template <typename TDeltaStore>
inline typename Size<TDeltaStore>::Type
insertionSize(typename DeltaValue<TDeltaStore, DeltaTypeSnp>::Type const & /*snp*/)
{
    return 1;
}

template <typename TDeltaStore>
inline typename Size<TDeltaStore>::Type
insertionSize(typename DeltaValue<TDeltaStore, DeltaTypeIns>::Type const & ins)
{
    return length(ins);
}

template <typename TDeltaStore>
inline typename Size<TDeltaStore>::Type
insertionSize(typename DeltaValue<TDeltaStore, DeltaTypeDel>::Type const & /*del*/)
{
    return 0;
}

template <typename TDeltaStore>
inline typename Size<TDeltaStore>::Type
insertionSize(typename DeltaValue<TDeltaStore, DeltaTypeSV>::Type const & sv)
{
    return length(sv.i2);
}

template <typename TStore, typename TPos, typename TTag>
inline typename Size<TStore>::Type
insertionSize(TStore const & store, TPos const pos, TTag const & /*tag*/)
{
    return insertionSize<TStore>(getDeltaStore(store, TTag())[pos]);
}

// ----------------------------------------------------------------------------
// Function netSize()
// ----------------------------------------------------------------------------

template <typename TStore, typename TPos, typename TTag>
inline typename MakeSigned<typename Size<TStore>::Type>::Type
netSize(TStore const & store, TPos const pos, TTag const & /*tag*/)
{
    typedef typename MakeSigned<typename Size<TStore>::Type>::Type TSignedSize;
    return static_cast<TSignedSize>(insertionSize(store, pos, TTag())) -
           static_cast<TSignedSize>(deletionSize(store, pos, TTag()));
}

}  // namespace impl

}  // namespace seqan

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_STORE_H_
