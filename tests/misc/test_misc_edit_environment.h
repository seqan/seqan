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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

#ifndef SEQAN_TESTS_MISC_TEST_MISG_EDIT_ENVIRONMENT_H_
#define SEQAN_TESTS_MISC_TEST_MISG_EDIT_ENVIRONMENT_H_

#include <sstream>

#include <seqan/misc/edit_environment.h>

// Test StringEnumerator for Hamming distance.
SEQAN_DEFINE_TEST(test_misc_edit_environment_string_enumerator_hamming)
{
    using namespace seqan;

    DnaString original = "CGAT";

    typedef StringEnumerator<DnaString, EditEnvironment<HammingDistance, 1> > THammingEnumerator;
    typedef Iterator<THammingEnumerator>::Type THammingIterator;

    std::stringstream result;

    THammingEnumerator hammingEnumerator(original);
    for (THammingIterator itH = begin(hammingEnumerator); !atEnd(itH); goNext(itH))
        result << *itH << '$';

    CharString expected = "AGAT$CGAT$GGAT$TGAT$CAAT$CCAT$CTAT$CGCT$CGGT$CGTT$CGAA$CGAC$CGAG$";
    SEQAN_ASSERT_EQ(expected, result.str());

    SEQAN_ASSERT_EQ(length(hammingEnumerator), 13u);
}

// More comprehensive tests for Hamming StringEnumerator Iterator.
SEQAN_DEFINE_TEST(test_misc_edit_environment_string_enumerator_iterator_hamming)
{
    using namespace seqan;

    DnaString original = "CGAT";

    typedef StringEnumerator<DnaString, EditEnvironment<HammingDistance, 1> > THammingEnumerator;
    typedef Iterator<THammingEnumerator, Standard>::Type THammingIterator;

    THammingEnumerator hammingEnumerator(original);

    THammingIterator it1 = begin(hammingEnumerator);
    goEnd(it1);
    SEQAN_ASSERT(atEnd(it1));
    goBegin(it1);
    SEQAN_ASSERT(it1 == begin(hammingEnumerator, Standard()));

    THammingIterator it2 = end(hammingEnumerator);
    SEQAN_ASSERT(atEnd(it2));
    goBegin(it2);
    SEQAN_ASSERT(it2 == begin(hammingEnumerator, Standard()));
}

// Test StringEnumerator for Edit distance.
SEQAN_DEFINE_TEST(test_misc_edit_environment_string_enumerator_edit)
{
    using namespace seqan;

    DnaString original = "CGAT";

    typedef StringEnumerator<DnaString, EditEnvironment<LevenshteinDistance, 1> > TEditEnumerator;
    typedef Iterator<TEditEnumerator>::Type TEditIterator;

    std::stringstream result;

    TEditEnumerator editEnumerator(original);
    for (TEditIterator itH = begin(editEnumerator); !atEnd(itH); goNext(itH))
        result << *itH << '$';

    CharString expected = "CGAT$CAAT$CCAT$CTAT$CGCT$CGGT$CGTT$GAT$CAT$CGT$CGA$CAGAT$CCGAT$CGGAT$CTGAT$CGAAT$CGCAT$CGGAT$CGTAT$CGAAT$CGACT$CGAGT$CGATT$";
    SEQAN_ASSERT_EQ(expected, result.str());

    SEQAN_ASSERT_EQ(length(editEnumerator), 23u);
}

// More comprehensive tests for Edit StringEnumerator Iterator.
SEQAN_DEFINE_TEST(test_misc_edit_environment_string_enumerator_iterator_edit)
{
    using namespace seqan;

    DnaString original = "CGAT";

    typedef StringEnumerator<DnaString, EditEnvironment<LevenshteinDistance, 1> > TEditEnumerator;
    typedef Iterator<TEditEnumerator, Standard>::Type TEditIterator;

    TEditEnumerator editEnumerator(original);

    TEditIterator it1 = begin(editEnumerator);
    goEnd(it1);
    SEQAN_ASSERT(atEnd(it1));
    goBegin(it1);
    SEQAN_ASSERT(it1 == begin(editEnumerator, Standard()));

    TEditIterator it2 = end(editEnumerator);
    SEQAN_ASSERT(atEnd(it2));
    goBegin(it2);
    SEQAN_ASSERT(it2 == begin(editEnumerator, Standard()));
}

#endif  // #ifndef SEQAN_TESTS_MISC_TEST_MISG_EDIT_ENVIRONMENT_H_
