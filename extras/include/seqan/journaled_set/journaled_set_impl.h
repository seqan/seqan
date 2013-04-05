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

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_IMPL_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_IMPL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class StringSet                                          [Owner<JournaledSet>]
// ----------------------------------------------------------------------------

template <typename TString>
class StringSet<TString, Owner<JournaledSet> >
{
public:
    typedef typename Position<StringSet>::Type TPosition_;
    typedef typename Host<TString>::Type THostString_;
    typedef String<TString>  TStrings_;
    typedef typename StringSetLimits<StringSet>::Type TLimits_;

    Holder<THostString_> _globalRefHolder;
    TStrings_   strings;
    TLimits_    limits;
    bool        limitsValid;

    StringSet() : limitsValid(true)
    {
        appendValue(limits,0);
    }

    ~StringSet() {}

    template <typename TPosition>
    inline typename Reference<StringSet>::Type
    operator[](TPosition const pos)
    {
        return value(*this, pos);
    }

    template <typename TPosition>
    inline typename Reference<StringSet const>::Type
    operator[](TPosition const pos) const
    {
        return value(*this, pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

template <typename TString>
struct Host<StringSet<TString, Owner<JournaledSet> > >
{
    typedef typename Host<TString>::Type Type;
};

template <typename TString>
struct Host<StringSet<TString, Owner<JournaledSet> > const >
{
    typedef typename Host<TString const>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TString, typename TPosition>
inline typename Reference<StringSet<TString, Owner<JournaledSet> > >::Type
value(StringSet<TString, Owner<JournaledSet> > & me,
      TPosition const pos)
{
    return me.strings[pos];
}

template <typename TString, typename TPosition>
inline typename Reference<StringSet<TString, Owner<JournaledSet> > const >::Type
value(StringSet<TString, Owner<JournaledSet> > const & me,
      TPosition const pos)
{
    return me.strings[pos];
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TExpand>
void appendValue(StringSet<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournaledSet> > & journalSet,
                 String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & newElement,
                 Tag<TExpand> tag)
{
    if (_validStringSetLimits(journalSet))
        appendValue(journalSet.limits, lengthSum(journalSet) + length(newElement), tag);
    appendValue(journalSet.strings, newElement, tag);
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TExpand>
void
appendValue(StringSet<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournaledSet> > & journalSet,
        String<TValue, THostSpec> & newElement,
        Tag<TExpand> tag)
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;

    TJournalString jrn(newElement);

    if (_validStringSetLimits(journalSet))
            appendValue(journalSet.limits, lengthSum(journalSet) + length(jrn), tag);
    appendValue(journalSet.strings, jrn, tag);
}

template <typename TString, typename TString2, typename TExpand>
void
appendValue(StringSet<TString, Owner<JournaledSet> > & journalSet,
            TString2 const & newElement,
            Tag<TExpand> tag)
{
    if (_validStringSetLimits(journalSet))
        appendValue(journalSet.limits, lengthSum(journalSet) + length(newElement), tag);
    appendValue(journalSet.strings, newElement, tag);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TString>
inline void
clear(StringSet<TString, Owner<JournaledSet> > & journalSet)
{
    clear(journalSet._globalRefHolder);
    clear(journalSet.strings);
    resize(journalSet.limits, 1, Exact());
    journalSet.limitsValid = true;
}

// ----------------------------------------------------------------------------
// Function erase()
// ----------------------------------------------------------------------------

template <typename TString, typename TPos>
inline typename Size<StringSet<TString, Owner<JournaledSet> > >::Type
erase(StringSet<TString, Owner<JournaledSet> > & journalSet, TPos pos)
{
    erase(journalSet.strings, pos);
    journalSet.limitsValid = false;
    return length(journalSet);
}

template <typename TString, typename TPos, typename TPosEnd>
inline typename Size<StringSet<TString, Owner<JournaledSet> > >::Type
erase(StringSet<TString, Owner<JournaledSet> > & journalSet, TPos pos, TPosEnd posEnd)
{
    erase(journalSet.strings, pos, posEnd);
    journalSet.limitsValid = false;
    return length(journalSet);
}

// --------------------------------------------------------------------------
// Function assignValue()
// --------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TPos>
inline void assignValue(
    StringSet<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournaledSet> > & journalSet,
    TPos pos,
    String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & newElement)
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
    typedef StringSet<TJournalString, Owner<JournaledSet> > TJournaledSet;
    typedef typename Size<TJournaledSet>::Type TSize;
    typedef typename StringSetLimits<TJournaledSet>::Type TLimits;
    typedef typename Value<TLimits>::Type TLimitValue;
    typedef typename MakeSigned<TLimitValue>::Type TSignedLimitValue;

    TSignedLimitValue oldSize = length(journalSet[pos]);
    set(journalSet[pos], newElement);
    if (_validStringSetLimits(journalSet))
    {
        TSignedLimitValue delta = (TSignedLimitValue)length(newElement) - oldSize;
        TSize size = length(journalSet);
        while (pos < size)
            journalSet.limits[++pos] += delta;
    }
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TPos>
inline void assignValue(
    StringSet<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournaledSet> > & journalSet,
    TPos pos,
    String<TValue, THostSpec> & newElement)
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
    typedef StringSet<TJournalString, Owner<JournaledSet> > TJournaledSet;
    typedef typename Size<TJournaledSet>::Type TSize;
    typedef typename StringSetLimits<TJournaledSet>::Type TLimits;
    typedef typename Value<TLimits>::Type TLimitValue;
    typedef typename MakeSigned<TLimitValue>::Type TSignedLimitValue;

    TSignedLimitValue oldSize = length(journalSet[pos]);
    TJournalString newJrn(newElement);
    set(journalSet[pos], newJrn);
    if (_validStringSetLimits(journalSet))
    {
        TSignedLimitValue delta = (TSignedLimitValue)length(newJrn) - oldSize;
        TSize size = length(journalSet);
        while (pos < size)
            journalSet.limits[++pos] += delta;
    }
}

template <typename TString, typename TPos,  typename TString2>
inline void assignValue(
    StringSet<TString, Owner<JournaledSet> > & journalSet,
    TPos pos,
    TString2 const & newElement)
{
    typedef StringSet<TString, Owner<JournaledSet> > TJournaledSet;
    typedef typename Size<TJournaledSet>::Type TSize;
    typedef typename StringSetLimits<TJournaledSet>::Type TLimits;
    typedef typename Value<TLimits>::Type TLimitValue;
    typedef typename MakeSigned<TLimitValue>::Type TSignedLimitValue;

    TSignedLimitValue oldSize = length(journalSet[pos]);
    assign(journalSet[pos], newElement);
    if (_validStringSetLimits(journalSet))
    {
        TSignedLimitValue delta = (TSignedLimitValue)length(newElement) - oldSize;
        TSize size = length(journalSet);
        while (pos < size)
            journalSet.limits[++pos] += delta;
    }
}

// ----------------------------------------------------------------------------
// Function globalReference()
// ----------------------------------------------------------------------------

/**
.Function.globalReference:
..summary:Returns the global reference sequence of a @Spec.Journaled Set@.
..class:Spec.Journaled Set
..cat:Sequences
..signature:globalReference(stringSet)
..param.stringSet: The string set that stores the sequences.
...type:Spec.Journaled Set
..returns:The global reference sequence of the @Spec.Journaled Set@.
...metafunction:Metafunction.Host
..see:Function.setGlobalReference
..see:Function.createGlobalReference
..include:seqan/journal_set.h
*/

template <typename TString>
inline typename Host<StringSet<TString, Owner<JournaledSet> > >::Type &
globalReference(StringSet<TString, Owner<JournaledSet> > const & journalSet)
{
    return value(journalSet._globalRefHolder);
}

// ----------------------------------------------------------------------------
// Function setGlobalReference()
// ----------------------------------------------------------------------------

/**
.Function.setGlobalReference:
..summary:Sets the global reference of a @Spec.Journaled Set@.
..class:Spec.Journaled Set
..cat:Sequences
..signature:setGlobalReference(stringSet, reference)
..param.stringSet: The string set that stores the sequences.
...type:Spec.Journaled Set
..param.reference: The new reference sequence of the Journaled Set.
...type:nolink:$Host<StringSet<TString, Owner<JournaledSet> > >::Type$
..remarks:Uses an @Class.Holder@ to store a reference to the new global reference sequence instead of copying it.
..see:Function.createGlobalReference
..see:Function.globalReference
..include:seqan/journal_set.h
*/
template <typename TString>
inline void
setGlobalReference(StringSet<TString, Owner<JournaledSet> > & journalSet,
                   typename Host<TString>::Type & newGlobalRef)
{
    setValue(journalSet._globalRefHolder, newGlobalRef);
}

// ----------------------------------------------------------------------------
// Function createGlobalReference()
// ----------------------------------------------------------------------------

/**
.Function.createGlobalReference:
..summary:Creates the global reference of a @Spec.Journaled Set@.
..class:Spec.Journaled Set
..cat:Sequences
..signature:createGlobalReference(stringSet, reference)
..param.stringSet: The string set that stores the sequences.
...type:Spec.Journaled Set
..param.reference: The new reference sequence of the Journaled Set.
...type:nolink:$Host<StringSet<TString, Owner<JournaledSet> > >::Type$
..remarks:Stores a copy of the passed global reference sequence.
..see:Function.setGlobalReference
..see:Function.globalReference
..include:seqan/journal_set.h
*/
template <typename TString>
inline void
createGlobalReference(StringSet<TString, Owner<JournaledSet> > & journalSet,
                   typename Host<TString>::Type const & newGlobalRef)
{
    create(journalSet._globalRefHolder, newGlobalRef);
}

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_IMPL_H_
