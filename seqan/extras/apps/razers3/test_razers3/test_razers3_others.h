/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de
  ===========================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  ===========================================================================
  Author: @@Your Name@@ <@@Your Email@@>
  ===========================================================================
  @@By its name, this would test the header template/others.h.@@
  ===========================================================================
*/

// @@Lines and sections with two at-chars are comments that should either
// be removed or replaced by tests defined by you.@@

// @@Replace symbol accordingly to your file name!@@
#ifndef TEST_RAZERS3_TEST_PARALLEL_OTHERS_H_
#define TEST_RAZERS3_TEST_PARALLEL_OTHERS_H_

#include <seqan/basic.h>     // @@For the testing infrastructure.@@

#define RAZERS_PARALLEL_READS_INDEPENDENT
#include "../razers_parallel_reads.h"  // @@Replace with your header under test.@@
#include <seqan/sequence.h>


struct TestHit
{
    int ndlSeqNo;
};

// @@This is an empty example test.  Replace it with your own test.@@
SEQAN_DEFINE_TEST(test_split_algorithm)
{
    using namespace seqan;

    String<int> tav;
    resize(tav, 4, 1, Exact());

    String<int> expTav;
    resize(expTav, 4, 1, Exact());

    expTav[0] = -1;
    expTav[1] = 2;
    expTav[2] = 1;
    expTav[3] = 1;
    updateTav(tav, 0);
    SEQAN_ASSERT_EQ_MSG(tav, expTav, "first operation failed");

    expTav[0] = -1;
    expTav[1] = -1;
    expTav[2] = 3;
    expTav[3] = 1;
    updateTav(tav, 1);
    SEQAN_ASSERT_EQ_MSG(tav, expTav, "second operation failed");

    expTav[0] = -1;
    expTav[1] = -1;
    expTav[2] = 4;
    expTav[3] = -1;
    updateTav(tav, 3);
    SEQAN_ASSERT_EQ_MSG(tav, expTav, "third operation failed");

    expTav[0] = -1;
    expTav[1] = -1;
    expTav[2] = 4;
    expTav[3] = -1;
    updateTav(tav, 2);
    SEQAN_ASSERT_EQ_MSG(tav, expTav, "fourth operation failed");

}

SEQAN_DEFINE_TEST(test_to_bucket)
{
    using namespace seqan;

    #pragma omp parallel
    {
    #pragma omp single
    {

        int id = 0;
        int result = toBucket<0>(id);
        SEQAN_ASSERT_EQ_MSG(result, 0, "id = 0, stage = 0");

        id = 1;
        result = toBucket<0>(id);
        SEQAN_ASSERT_EQ_MSG(result, 1, "id = 1, stage = 0");

        id = 256;
        result = toBucket<0>(id);
        SEQAN_ASSERT_EQ_MSG(result, 0, "id = 256, stage = 0");

        id = 259;
        result = toBucket<0>(id);
        SEQAN_ASSERT_EQ_MSG(result, 3, "id = 259, stage = 0");

        id = 256;
        result = toBucket<1>(id);
        SEQAN_ASSERT_EQ_MSG(result, 1, "id = 256, stage = 1");

        id = 259;
        result = toBucket<1>(id);
        SEQAN_ASSERT_EQ_MSG(result, 1, "id = 259, stage = 1");

        id = 1333;
        result = toBucket<0>(id);
        SEQAN_ASSERT_EQ_MSG(result, 53, "id = 1333, stage = 0");

        id = 1333;
        result = toBucket<1>(id);
        SEQAN_ASSERT_EQ_MSG(result, 5, "id = 1333, stage = 1");

    }
    }
}

SEQAN_DEFINE_TEST(test_radix_pass)
{
    using namespace seqan;

    #pragma omp parallel
    {
    #pragma omp single
    {

        String<TestHit> hits;
        resize(hits, 4, Exact());
        hits[0].ndlSeqNo = 1;
        hits[1].ndlSeqNo = 12;
        hits[2].ndlSeqNo = 515;
        hits[3].ndlSeqNo = 7;

        String<TestHit> sorted;
        resize(sorted, 4, Exact());

        myRadixPass<0>(sorted, hits);

        SEQAN_ASSERT_EQ_MSG(sorted[0].ndlSeqNo, 1, "radix pass 1, elem 0");
        SEQAN_ASSERT_EQ_MSG(sorted[1].ndlSeqNo, 515, "radix pass 1, elem 1");
        SEQAN_ASSERT_EQ_MSG(sorted[2].ndlSeqNo, 7, "radix pass 1, elem 2");
        SEQAN_ASSERT_EQ_MSG(sorted[3].ndlSeqNo, 12, "radix pass 1, elem 3");

        myRadixPass<1>(hits, sorted);

        SEQAN_ASSERT_EQ_MSG(hits[0].ndlSeqNo, 1, "radix pass 2, elem 0");
        SEQAN_ASSERT_EQ_MSG(hits[1].ndlSeqNo, 7, "radix pass 2, elem 1");
        SEQAN_ASSERT_EQ_MSG(hits[2].ndlSeqNo, 12, "radix pass 2, elem 2");
        SEQAN_ASSERT_EQ_MSG(hits[3].ndlSeqNo, 515, "radix pass 2, elem 3");

    }
    }
}

SEQAN_DEFINE_TEST(test_radix_sort)
{
    using namespace seqan;

    #pragma omp parallel
    {
    #pragma omp single
    {

        String<TestHit> hits;
        resize(hits, 4, Exact());
        hits[0].ndlSeqNo = 1;
        hits[1].ndlSeqNo = 12;
        hits[2].ndlSeqNo = 515;
        hits[3].ndlSeqNo = 7;

        myRadixSort(hits);

        SEQAN_ASSERT_EQ_MSG(hits[0].ndlSeqNo, 1, "radix sort, elem 0");
        SEQAN_ASSERT_EQ_MSG(hits[1].ndlSeqNo, 7, "radix sort, elem 1");
        SEQAN_ASSERT_EQ_MSG(hits[2].ndlSeqNo, 12, "radix sort, elem 2");
        SEQAN_ASSERT_EQ_MSG(hits[3].ndlSeqNo, 515, "radix sort, elem 3");

        resize(hits, 400, Exact());
        for (int i = 0; i < (int)length(hits); ++i)
            hits[i].ndlSeqNo = 500 - i;

        myRadixSort(hits);

        for (int i = 0; i < (int)length(hits); ++i)
            SEQAN_ASSERT_EQ_MSG(hits[i].ndlSeqNo, (101 + i), "radix sort, loop");


    }
    }
}

// @@Replace symbol accordingly to your file name!@@
#endif  // TEST_RAZERS3_TEST_PARALLEL_OTHERS_H_
