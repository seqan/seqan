// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2024, Knut Reinert, FU Berlin
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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Implementation of the StringSet specialization Dependent<Tight>, the
// the default specialization of Dependent<>.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_STRING_SET_DEPENDENT_TIGHT_H_
#define SEQAN_SEQUENCE_STRING_SET_DEPENDENT_TIGHT_H_

namespace seqan2 {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(holtgrew): Change name of specialization to Dependent StringSet.

/*!
 * @class DependentStringSet Dependent StringSet
 * @extends StringSet
 * @headerfile <seqan/sequence.h>
 * @brief StringSet implementation that only stores pointers to strings in other string sets.
 *
 * @signature template <typename TString, typename TSpec>
 *            class StringSet<TString, Depedent<TSpec> >;
 *
 * @tparam TString The type of the string to store in the string set.
 * @tparam TSpec   Tag for further specializing the string set.
 *
 * Important: This is an abstract class. Use one of the following specializations: @link TightDependentStringSet @endlink and
 * @link GenerousDependentStringSet @endlink.
 *
 * A Dependent StringSet (DSS) can be used like a normal StringSet while internally storing only pointers to a source set.
 *
 * @section Modifying a Dependent StringSet
 *
 * (1) Removing a sequence from a DSS only removes the pointer and does not change the source set.
 *
 * (2) Appending a sequence to a DSS appends a pointer to that sequence and does not change the source set.
 *
 * (3) Assigning a sequence to a position (or id) of a DSS, dereferences the pointer first and does thus additionally change the source set.
 *
 * (4) Accessing the DSS at a position (or id) dereferences the pointer, and when stored as a reference, modifications also lead to a change in the source set.
 *
 * @section Position vs. Id
 *
 * When a sequence is removed in a DSS, the positions of pointers shift and do not represent the exact position in the original source set anymore.
 * To distinguish between new and original positions, we introduce the term <tt>id</tt> which refers to the original positions in the source set.
 * Every modification of a DSS can be either based on the id (function ending on <tt>ById()</tt>) or the position depending on the behaviour you want to realize.
 *
 * The following figure illustrates the behaviour when removing a sequence:
 *
 * <img src="position_vs_id.png" title="Impact to positions and ids on removal of an entry in an DSS" width="900">
 *
 * @section Tight vs. Generous
 *
 * The two different specializations <tt>Tight</tt> and <tt>Generous</tt> provide the same functionality but behave slightly different
 * concerning run time and certain errors (e.g. index out of range). See the correspoinding documentation pages @link TightDependentStringSet @endlink
 * and @link GenerousDependentStringSet @endlink for further details.
 */

/*!
 * @class TightDependentStringSet Tight Dependent StringSet
 * @extends DependentStringSet
 * @headerfile <seqan/sequence.h>
 * @brief Very space efficient Dependent StringSet implementation.
 *
 * @signature template <typename TString>
 *            class StringSet<TString, Depedent<Tight> >;
 *
 * @tparam TString The type of the string to store in the string set.
 *
 * The Tight Dependent StringSet stores pointers to a source set, enabling the user to perform deletions and additions to the set without
 * changing the original source set (See @link DependentStringSet @endlink for further details).
 *
 * @section Run time and Memory
 *
 * When a value is removed from the Tight Dependent StringSet, the array of pointers is resized accordingly.
 * Therefore, in order to call sequences by id or position, the stringset keeps a id-to-position map, which affects the run time complexity of the following functions:
 *
 * - value() or operator []: O(1)
 *
 * - getValueById(): O(log(n))
 *
 * - removeValueById(): O(log(n))
 *
 * The memory consumption is linear to the number of pointers.
 *
 * See @link GenerousDependentStringSet @endlink for a Dependent StringSet implementation that allows for more
 * efficient access to strings in the container via ids at the cost of higher memory usage.
 *
 * @section Accessing non-existing entries results in undefined behaviour
 *
 * Because the Tight Dependent StringSet keeps every array "tight", every entry that is being removed, is actually deleted (in contrast to <tt>Generous</tt>)
 * and accessing it will result in undefined behaviour.
 *
 */

// Default id holder string set
template <typename TSpec = Tight>
struct Dependent;

// StringSet with individual sequences in a tight string of string pointers and corr. IDs
template <typename TString>
class StringSet<TString, Dependent<Tight> >
{
public:
    typedef String<TString *>                           TStrings;
    typedef typename Id<StringSet>::Type                TIdType;
    typedef typename Position<StringSet>::Type          TPosition;
    typedef String<TIdType>                             TIds;
    typedef std::map<TIdType, TPosition>                TIdPosMap;
    typedef typename StringSetLimits<StringSet>::Type   TLimits;
    typedef typename Concatenator<StringSet>::Type      TConcatenator;

    TIdType         lastId;
    TStrings        strings;
    TIds            ids;
    TIdPosMap       id_pos_map;
    TLimits         limits;
    bool            limitsValid;        // is true if limits contains the cumulative sum of the sequence lengths
    TConcatenator   concat;

    StringSet() :
        lastId(0),
        limitsValid(true)
    {
        _initStringSetLimits(*this);
    }

    template <typename TDefault>
    StringSet(StringSet<TString, Owner<TDefault> > const & _other) :
        lastId(0),
        limitsValid(true)
    {
        _initStringSetLimits(*this);
        for (unsigned int i = 0; i < length(_other); ++i)
            appendValue(*this, _other[i]);
    }

    // ----------------------------------------------------------------------
    // Subscription operators; have to be defined in class def.
    // ----------------------------------------------------------------------

    template <typename TPos>
    inline typename Reference<StringSet>::Type
    operator[] (TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<StringSet const>::Type
    operator[] (TPos pos) const
    {
        return value(*this, pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function appendValue()
// --------------------------------------------------------------------------

template <typename TString, typename TExpand >
inline void appendValue(
    StringSet<TString, Dependent<Tight> > & me,
    TString const & obj,
    Tag<TExpand> tag)
{
    typedef typename Position<StringSet<TString, Dependent<Tight> > >::Type TPos;
    appendValue(me.limits, lengthSum(me) + length(obj), tag);
    typedef typename StringSet<TString, Dependent<Tight> >::TIdType TIdType;
    appendValue(me.strings, const_cast<TString*>(&obj));
    TIdType last = me.lastId++;
    appendValue(me.ids, last, tag);
    me.id_pos_map.insert(std::make_pair(last, (TPos)(length(me.strings) - 1)));
}

// --------------------------------------------------------------------------
// Function clear()
// --------------------------------------------------------------------------

template <typename TString >
inline void clear(StringSet<TString, Dependent<Tight> >& me)
{
    clear(me.strings);
    me.id_pos_map.clear();
    resize(me.limits, 1, Exact());
    me.limitsValid = true;

    clear(me.ids);
    me.lastId = 0;
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename TString, typename TPos >
inline typename Reference<StringSet<TString, Dependent<Tight> > >::Type
value(StringSet<TString, Dependent<Tight> >& me, TPos pos)
{
    return *me.strings[pos];
}

template <typename TString, typename TPos >
inline typename Reference<StringSet<TString, Dependent<Tight> > const >::Type
value(StringSet<TString, Dependent<Tight> >const & me, TPos pos)
{
    return *me.strings[pos];
}

// --------------------------------------------------------------------------
// Function getValueById()
// --------------------------------------------------------------------------

template <typename TString, typename TId>
inline typename Reference<StringSet<TString, Dependent<Tight> > >::Type
getValueById(StringSet<TString, Dependent<Tight> > & me,
            TId const id)
{
    SEQAN_ASSERT_GT_MSG(me.id_pos_map.count(id), 0u, "String id must be known!");
    return (value(me, me.id_pos_map.find(id)->second));
}

// --------------------------------------------------------------------------
// Function assignValueById()
// --------------------------------------------------------------------------

template<typename TString, typename TString2>
inline typename Id<StringSet<TString, Dependent<Tight> > >::Type
assignValueById(StringSet<TString, Dependent<Tight> >& me,
                TString2& obj)
{
    appendValue(me, obj);
    SEQAN_ASSERT_EQ(length(me.limits), length(me) + 1);
    return positionToId(me, length(me.strings) - 1);
}


template<typename TString, typename TId1>
inline typename Id<StringSet<TString, Dependent<Tight> > >::Type
assignValueById(StringSet<TString, Dependent<Tight> >& me,
                TString& obj,
                TId1 id)
{
    typedef StringSet<TString, Dependent<Tight> > TStringSet;
    typedef typename TStringSet::TIdPosMap::const_iterator TIter;
    typedef typename Id<TStringSet>::Type TId;

    if (me.lastId < (TId) id) me.lastId = (TId) (id + 1);

    TIter pos = me.id_pos_map.find(id);
    if (pos != me.id_pos_map.end()) {
        me.strings[pos->second] = &obj;
        me.limitsValid = false;
        return id;
    }
    appendValue(me.strings, &obj);
    appendValue(me.ids, id);
    me.id_pos_map.insert(std::make_pair(id, length(me.strings) - 1));
    appendValue(me.limits, lengthSum(me) + length(obj));
    return id;
}

// --------------------------------------------------------------------------
// Function removeValueById()
// --------------------------------------------------------------------------

template<typename TString, typename TId>
inline void
removeValueById(StringSet<TString, Dependent<Tight> >& me, TId const id)
{
    typedef StringSet<TString, Dependent<Tight> > TStringSet;
    typedef typename Size<TStringSet>::Type TSize;
    typedef typename TStringSet::TIdPosMap::iterator TIter;

    SEQAN_ASSERT_EQ(length(me.limits), length(me) + 1);
    TIter pos = me.id_pos_map.find(id);
    if (pos != me.id_pos_map.end()) {
        TSize remPos = pos->second;
        erase(me.strings, remPos);
        erase(me.ids, remPos);
        me.id_pos_map.erase(pos);
        resize(me.limits, length(me.limits) - 1, Generous());

        for(TIter itChange = me.id_pos_map.begin(); itChange != me.id_pos_map.end(); ++itChange) {
            if (itChange->second > remPos) --(itChange->second);
        }
    }
    SEQAN_ASSERT_EQ(length(me.limits), length(me) + 1);
}

// --------------------------------------------------------------------------
// Function positionToId()
// --------------------------------------------------------------------------

template <typename TString, typename TPos>
inline typename Id<StringSet<TString, Dependent<Tight> > >::Type
positionToId(StringSet<TString, Dependent<Tight> > & me,
            TPos const pos)
{
    return me.ids[pos];
}

template <typename TString, typename TPos>
inline typename Id<StringSet<TString, Dependent<Tight> > >::Type
positionToId(StringSet<TString, Dependent<Tight> > const & me,
            TPos const pos)
{
    return me.ids[pos];
}

// --------------------------------------------------------------------------
// Function idToPosition()
// --------------------------------------------------------------------------

template <typename TString, typename TId>
inline typename Position<StringSet<TString, Dependent<Tight> > >::Type
idToPosition(StringSet<TString, Dependent<Tight> > const & me,
            TId const id)
{
    return me.id_pos_map.find(id)->second;
/*
    for(unsigned i = 0; i < length(me.ids); ++i)
        if ((TId) me.ids[i] == id)
            return i;
    return 0;
    */
}

}  // namespace seqan2

#endif  // #ifndef SEQAN_SEQUENCE_STRING_SET_DEPENDENT_TIGHT_H_
