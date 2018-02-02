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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains the Verifier class.
// ==========================================================================

#ifndef SEQAN_APPS_SEARCHJOIN_VERIFIER_H_
#define SEQAN_APPS_SEARCHJOIN_VERIFIER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/find.h>

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Verifier
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec = void>
struct Verifier
{
    typedef Db<TText>           TDb;
    typedef Db<TText, TSpec>    TDbQuery;

    TDb /* const */         & db;
    Holder<TDbQuery>        query;
    unsigned long           lengthFilter;
    unsigned long           verifications;

    Verifier(TDb /* const */ & db) :
        db(db),
        query(),
        lengthFilter(0),
        verifications(0)
    {}

    template <typename TDelegate>
    inline bool
    operator() (typename Size<TDb>::Type dbId,
                typename Size<TDbQuery>::Type queryId,
                TDelegate & delegate)
    {
        return _verify(*this, dbId, queryId, delegate);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _verify()                                                [Verifier]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TDelegate>
inline bool
_verify(Verifier<TText, TSpec> & verifier,
        typename Size<Db<TText> >::Type dbId,
        typename Size<Db<TText, TSpec> >::Type queryId,
        TDelegate & delegate)
{
    //typedef Db<TText> /* const */                           TDb;
    //typedef typename Size<TDb>::Type                        TDbSize;
    typedef typename Value<TText>::Type                     TTextReference;
    typedef typename Size<TTextReference>::Type             TTextSize;
    typedef typename Size<TTextReference>::Type             TErrors;
    //typedef Finder<TTextReference>                          TFinder;
    //typedef MyersUkkonenGlobalBanded                        TAlgorithmSpec;
    //typedef PatternState_<TTextReference, TAlgorithmSpec>   TPatternState;

    // Get texts from database.
    TTextReference text1 = verifier.db.text[dbId];
    TTextReference text2 = value(verifier.query).text[queryId];
    TTextSize length1 = length(text1);
    TTextSize length2 = length(text2);
    TErrors maxErrors = getErrors(value(verifier.query), queryId);

    // Check that texts do not differ more than maxErrors in length.
    if (_max(length1, length2) - _min(length1, length2) > maxErrors)
    {
        SEQAN_OMP_PRAGMA(atomic)
        verifier.lengthFilter++;

        return false;
    }

    SEQAN_OMP_PRAGMA(atomic)
    verifier.verifications++;

    if (maxErrors == 0)
    {
        if (text1 != text2)
            return false;
    }
    else
    {
        TErrors errors = _computeEditDistanceBanded(text1, text2, maxErrors);

        if (errors > maxErrors)
            return false;
    }

    delegate(dbId, queryId);

    return true;
}

#endif  // #ifndef SEQAN_APPS_SEARCHJOIN_VERIFIER_H_
