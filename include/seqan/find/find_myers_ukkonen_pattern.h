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

#ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_PATTERN_H
#define SEQAN_HEADER_FIND_MYERS_UKKONEN_PATTERN_H

#include <seqan/find/find_myers_ukkonen_base.h>

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// State Data
//////////////////////////////////////////////////////////////////////////////

// TODO: should go elsewhere
// template <typename TNeedle, typename TSpec>
// class Pattern{};
// TODO: should go elsewhere
template <typename TNeedle, typename TSpec>
class PatternState_{};

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
class PatternState_<TNeedle, Myers<TSpec, False, TFindBeginPatternSpec> > {};

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
class PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> >:
    public MyersSmallState_<TNeedle, TSpec>
{
public:
    typedef MyersSmallState_<TNeedle, TSpec>    TSmallState;
    typedef MyersLargeState_<TNeedle, TSpec>    TLargeState;
    typedef typename TSmallState::TWord         TWord;

    enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

    TLargeState *largeState;

//____________________________________________________________________________

    PatternState_():
        largeState(NULL) {}

    PatternState_(PatternState_ const & other):
        TSmallState(other),
        largeState(NULL)
    {
        if (other.largeState)
            largeState = new TLargeState(*other.largeState);
    }

    // Add move constructor.
    PatternState_(PatternState_ && other):
        TSmallState(std::forward<PatternState_>(other)),
        largeState(NULL)
    {
        largeState = other.largeState;
        other.largeState = NULL;
    }

    ~PatternState_()
    {
        delete largeState;
    }

    PatternState_ &
    operator = (PatternState_ const & other)
    {
        TSmallState::operator=(other);
        if (other.largeState)
        {
            if (largeState == NULL)
                largeState = new TLargeState;
            (*largeState) = *(other.largeState);
        } else {
            delete largeState;
            largeState = NULL;
        }
        return *this;
    }

    // Add move assignment.
    PatternState_&
    operator=(PatternState_ && other)
    {
        TSmallState::operator=(std::forward<PatternState_>(other));
        delete largeState;
        largeState = NULL; // *this is now in a valid state in case of self assignment.
        largeState = other.largeState;
        other.largeState = NULL;  // We moved the data.
        return *this;
    }

};

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
class Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >:
    public MyersSmallPattern_<TNeedle, TSpec>,
    public FindBegin_<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >,
    public PatternState_<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >
{

public:
    typedef MyersSmallPattern_<TNeedle, TSpec>      TSmallPattern;
    typedef MyersLargePattern_<TNeedle, TSpec>      TLargePattern;
    typedef typename TSmallPattern::TWord           TWord;

    enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

    typedef PatternState_<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPatternState;

    TLargePattern *largePattern;    // extra preprocessing info for large patterns

//____________________________________________________________________________

    Pattern():
        largePattern(NULL) {}

    Pattern(int _limit):
        largePattern(NULL)
    {
        setScoreLimit(*this, _limit);
    }

    Pattern(Pattern const & other) :
        TSmallPattern(other),
        TPatternState(other),
        largePattern(NULL)
    {
        if (other.largePattern)
            largePattern = new TLargePattern(*other.largePattern);
    }

    // Add move constructor.
    Pattern(Pattern && other) :
        TSmallPattern(std::forward<Pattern>(other)),
        TPatternState(std::forward<Pattern>(other)),
        largePattern(NULL)
    {
        largePattern = other.largePattern;
        other.largePattern = NULL;
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 && ndl, int _limit = -1) :
        largePattern(NULL)
    {
        setScoreLimit(*this, _limit);
        setHost(*this, std::forward<TNeedle2>(ndl));
    }

    ~Pattern()
    {
        delete largePattern;
    }

    Pattern &
    operator = (Pattern const & other)
    {
        TSmallPattern::operator=(other);
        TPatternState::operator=(other);
        if (other.largePattern)
        {
            if (largePattern == NULL)
                largePattern = new TLargePattern;
            (*largePattern) = *(other.largePattern);
        } else {
            delete largePattern;
            largePattern = NULL;
        }
        return *this;
    }

    // Add move assignment.
    Pattern &
    operator= (Pattern && other)
    {
        TSmallPattern::operator=(std::forward<Pattern>(other));
        TPatternState::operator=(std::forward<Pattern>(other));
        delete largePattern;  // delete old data.
        largePattern = NULL;  // *this is now in a valid state in case of self-assignment.
        largePattern = other.largePattern;  // Move data.
        other.largePattern = NULL;  // Reset other to default.
        return *this;
    }
};

//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct FindBeginPatternSpec <Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >
{
    typedef TFindBeginPatternSpec Type;
};
template <typename TNeedle, typename THasState, typename TFindBeginPatternSpec>
struct FindBeginPatternSpec <Pattern<TNeedle, Myers<FindPrefix, THasState, TFindBeginPatternSpec> > >
{// no find begin for FindPrefix
    typedef void Type;
};


template <typename TPattern>
struct PatternState {};

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct PatternState<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >
{
    typedef PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > Type;
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________
/*!
 * @fn MyersPattern#scoreLimit
 * @headerfile <seqan/find.h>
 * @brief The minimal score a match must reach in approximate searching.
 *
 * @signature TScoreValue scoreLimit(pattern);
 *
 * @param[in] pattern The pattern to query.
 *
 * @return TScoreValue The score limit value.
 */

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline int
scoreLimit(PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const & state)
{
    return - (int) state.maxErrors;
}

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline int
scoreLimit(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const & pattern)
{
    return - (int) pattern.maxErrors;
}

//____________________________________________________________________________
/*!
 * @fn MyersPattern#setSoreLimit
 * @headerfile <seqan/find.h>
 * @brief Set the minimal score a match must reach in approximate serach.
 *
 * @signature void setScoreLimit(pattern, limit);
 *
 * @param[in,out] pattern The pattern to set the limit for.
 * @param[in]     limit   The limit score value to set.
 *
 * @return TScoreValue The score limit value.
 */

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TScoreValue>
inline void
setScoreLimit(PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & state,
           TScoreValue minScore)
{
    // we need to convert the minimal score into a maximal penalty
    // that is why minScore is negated
    state.maxErrors = -minScore;
}

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TScoreValue>
inline void
setScoreLimit(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern,
           TScoreValue minScore)
{
    // we need to convert the minimal score into a maximal penalty
    // that is why minScore is negated
    pattern.maxErrors = -minScore;
}

//____________________________________________________________________________
/*!
 * @fn MyersPattern#getScore
 * @headerfile <seqan/find.h>
 * @brief Score of the last found match in approximate searching.
 *
 * @signature TScoreValue getScore(pattern);
 *
 * @param[in] pattern A myers pattern that can be used for approximate searching.
 *
 * @return TScoreValue The score of the last match found using <tt>pattern</tt>.  If no match was found then the value
 *                     is undefined.
 */

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline int
getScore(PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const & state)
{
    return -(int)state.errors;
}

template<typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline int
getScore(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const & state)
{
    return -(int)state.errors;
}

} // namespace seqan

#endif // #ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_PATTERN_H
