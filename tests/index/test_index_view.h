// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
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

#ifndef TESTS_INDEX_VIEW_H_
#define TESTS_INDEX_VIEW_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function compareTreeIterators()
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Copied from test_stree_iterators.h

template <typename TIndex, typename TIndexView>
void compareTreeIterators(TIndex &index, TIndexView & /*indexView */)
{
    typedef typename Iterator<TIndex, TopDown<ParentLinks<> > >::Type       TIter;
    typedef typename Iterator<TIndexView, TopDown<ParentLinks<> > >::Type   TIterView;
    typedef Factory<TIter>                                                  TFactory;
    typedef typename View<TFactory>::Type                                   TTFactoryView;

    TFactory factory(index, 1u, length(index));
    TTFactoryView factoryView = view(factory);
    TIterView it2 = getObject(factoryView, 0u);
    TIter it1(index);

    while (!atEnd(it1) && !atEnd(it2))
    {
        SEQAN_ASSERT_EQ(representative(it1), representative(it2));
        SEQAN_ASSERT_EQ(parentEdgeLabel(it1), parentEdgeLabel(it2));
        SEQAN_ASSERT_EQ(countOccurrences(it1), countOccurrences(it2));
        SEQAN_ASSERT_EQ(isRoot(it1), isRoot(it2));
        SEQAN_ASSERT_EQ(isLeaf(it1), isLeaf(it2));
        goNext(it1);
        goNext(it2);
    }

    SEQAN_ASSERT_EQ(atEnd(it1), atEnd(it1));
}

// ----------------------------------------------------------------------------
// Function compareIndexFibres()
// ----------------------------------------------------------------------------

template <typename TIndex1, typename TIndex2, typename TFibre>
void compareIndexFibres(TIndex1 &index1, TIndex2 &index2, TFibre const fibre)
{
    typedef typename Fibre<TIndex1, TFibre const>::Type   TFibre1;
    typedef typename Fibre<TIndex2, TFibre const>::Type   TFibre2;

    TFibre1 & fibre1 = getFibre(index1, fibre);
    TFibre2 & fibre2 = getFibre(index2, fibre);

    SEQAN_ASSERT_EQ(length(fibre1), length(fibre2));
    for (unsigned i = 0; i < length(fibre1); ++i)
        SEQAN_ASSERT_EQ(fibre1[i], fibre2[i]);
}

// ============================================================================
// Tests
// ============================================================================

// ----------------------------------------------------------------------------
// Test test_index_view_basic
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_index_view_basic)
{
    typedef CharString                                      TText;
    typedef Index<TText, IndexSa<> >                        TIndex;
    typedef typename View<CharString>::Type                 TTextView;
    typedef typename View<Index<TText, IndexSa<> > >::Type  TIndexView;

    // ------------------------------------------------------------------------

    // Construct a text and its view.
    TText text("text");
    TTextView textView(text);

    // ------------------------------------------------------------------------

    // Construct an index and its view.
    TIndex index(text);
    indexCreate(index, FibreSA());
    TIndexView indexView = view(index);

    // ------------------------------------------------------------------------

    // NOTE(esiragusa): The index view constructor must not be called with the text view.
//    TIndexView indexView2(textView);

    // NOTE(esiragusa): The index view constructor interprets the original index as the text.
//    TIndexView indexView3(index);

    // ------------------------------------------------------------------------

    // Assert that the Text fibre and its view behave equally.
    compareIndexFibres(index, indexView, FibreText());

    // ------------------------------------------------------------------------

    // Assert that the SA fibre and its view behave equally.
    compareIndexFibres(index, indexView, FibreSA());

    // ------------------------------------------------------------------------

    // Assert that the index iterator and its view behave equally.
    compareTreeIterators(index, indexView);
}

#endif  // TESTS_INDEX_VIEW_H_
