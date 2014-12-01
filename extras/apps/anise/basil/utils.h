// ==========================================================================
//                                 BASIL
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_HOLTGREW_APPS_SV_SCANNER_UTILS_H_
#define SANDBOX_HOLTGREW_APPS_SV_SCANNER_UTILS_H_

#include <seqan/sequence.h>
#include <seqan/store.h>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function std::hash() for seqan::String.
// ----------------------------------------------------------------------------

namespace std
{

// ----------------------------------------------------------------------------
// Class std::hash<seqan::String<> >
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct hash<seqan::String<TValue, TSpec> >
{
    typedef size_t result_type;

    result_type operator()(seqan::String<TValue, TSpec> const & str) const
    {
        using seqan::ordValue;
        result_type res = 0;
        for (auto && c : str)
            res = 5 * res + ordValue(c);
        return res;
    }
};

}  // namespace std

// ----------------------------------------------------------------------------
// Class HashPair
// ----------------------------------------------------------------------------

// Replacement for std::hash<> for a pair of seqan::CharString and bool.

struct HashPair
{
    typedef size_t result_type;

    result_type operator()(std::pair<seqan::CharString, bool> const & p) const
    {
        std::hash<seqan::CharString> h;
        return h(p.first);
    }
};

// ----------------------------------------------------------------------------
// Function std::copy_if().
// ----------------------------------------------------------------------------

#if __cplusplus == 1

namespace std {

// The function copy_if was dropped from the standard library by accident.
//
// We add it back here.
template<typename In, typename Out, typename Pred>
Out copy_if(In first, In last, Out res, Pred Pr)
{
    while (first != last) {
        if (Pr(*first))
            *res++ = *first;
        ++first;
    }
    return res;
}

}  // namespace std

#endif  // #if __cplusplus == 1

// ----------------------------------------------------------------------------
// Function startsWithClipping().
// ----------------------------------------------------------------------------

// Returns whether the CIGAR string starts with a clipping character.

template <typename TOperation, typename TCount, typename TSpec>
inline bool startsWithClipping(seqan::String<seqan::CigarElement<TOperation, TCount>, TSpec> const & cigarString,
                               unsigned minCount = 0)
{
    if (empty(cigarString))
        return false;
    return ((front(cigarString).operation == 'S') || (front(cigarString).operation == 'H')) &&
            (front(cigarString).count >= minCount);
}

// ----------------------------------------------------------------------------
// Function startsWithClipping().
// ----------------------------------------------------------------------------

// Returns whether the CIGAR string ends with a clipping character.

template <typename TOperation, typename TCount, typename TSpec>
inline bool endsWithClipping(seqan::String<seqan::CigarElement<TOperation, TCount>, TSpec> const & cigarString,
                             unsigned minCount = 0)
{
    if (empty(cigarString))
        return false;
    return ((back(cigarString).operation == 'S') || (back(cigarString).operation == 'H')) &&
            (back(cigarString).count >= minCount);
}

// ----------------------------------------------------------------------------
// Function hasClipping()
// ----------------------------------------------------------------------------

inline bool hasClipping(seqan::BamAlignmentRecord const & record, unsigned minCount = 0)
{
    return startsWithClipping(record.cigar, minCount) || endsWithClipping(record.cigar, minCount);
}

// ----------------------------------------------------------------------------
// Function clippingLength()
// ----------------------------------------------------------------------------

// Returns length of clipping.

inline int clippingLength(seqan::BamAlignmentRecord const & record, unsigned minCount = 0)
{
    if (startsWithClipping(record.cigar, minCount))
        return front(record.cigar).count;
    else if (endsWithClipping(record.cigar, minCount))
        return back(record.cigar).count;
    else
        return 0;
}

// ----------------------------------------------------------------------------
// Function clippingPosition()
// ----------------------------------------------------------------------------

// Returns the position of the clipping in the reference.

inline int clippingPosition(seqan::BamAlignmentRecord const & record, unsigned minCount = 0)
{
    if (startsWithClipping(record.cigar, minCount))
        return record.beginPos;
    else if (endsWithClipping(record.cigar, minCount))
        return record.beginPos + getAlignmentLengthInRef(record);
    else
        SEQAN_FAIL("Record's CIGAR string has no clipping!");
    return -1;
}

#endif  // #ifndef SANDBOX_HOLTGREW_APPS_SV_SCANNER_UTILS_H_
