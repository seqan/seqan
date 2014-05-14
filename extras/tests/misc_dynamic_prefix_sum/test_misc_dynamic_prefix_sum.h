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
// Unit tests for dynamic prefix sums.
// ==========================================================================

#ifndef EXTRAS_TESTS_MISC_DYNAMIC_PREFIX_SUM_TEST_MISC_DYNAMIC_PREFIX_SUM_H_
#define EXTRAS_TESTS_MISC_DYNAMIC_PREFIX_SUM_TEST_MISC_DYNAMIC_PREFIX_SUM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/random.h>

#include <seqan/misc_dynamic_prefix_sum.h>

using namespace seqan;

// Test constructor.
SEQAN_DEFINE_TEST(test_misc_dynamic_prefix_sum_constructor)
{
    typedef DynamicPrefixSumTree<unsigned, int, 32> TPrefixSumTree;
    typedef typename Value<TPrefixSumTree>::Type TNode;

    TPrefixSumTree prefixSumTree;

    SEQAN_ASSERT_NEQ(prefixSumTree.rootPtr, (TNode*) 0);
    SEQAN_ASSERT_EQ(prefixSumTree.rootPtr->_isLeaf, true);

    SEQAN_ASSERT_EQ(capacity(prefixSumTree.rootPtr->_keyTab), DpsHelper_<32>::MAX_KEY_NUMBER);
    SEQAN_ASSERT_EQ(capacity(prefixSumTree.rootPtr->_cargoTab), DpsHelper_<32>::MAX_KEY_NUMBER);
    SEQAN_ASSERT_EQ(capacity(prefixSumTree.rootPtr->_childTab), DpsHelper_<32>::MAX_CHILD_NUMBER);
}

SEQAN_DEFINE_TEST(test_misc_dynamic_prefix_sum_empty)
{
    typedef DynamicPrefixSumTree<unsigned, int> TPrefixSumTree;
    typedef typename Value<TPrefixSumTree>::Type TNode;

    TPrefixSumTree tree;

    SEQAN_ASSERT_EQ(empty(tree), true);

    appendValue(tree.rootPtr->_keyTab, 0);
    appendValue(tree.rootPtr->_cargoTab, 1);
    SEQAN_ASSERT_EQ(empty(tree), false);

    appendValue(tree.rootPtr->_childTab, (TNode*) 0);
    SEQAN_ASSERT_EQ(empty(tree), false);

    clear(tree.rootPtr->_keyTab);
    SEQAN_ASSERT_EQ(empty(tree), false);

    clear(tree.rootPtr->_childTab);
    SEQAN_ASSERT_EQ(empty(tree), true);
}

SEQAN_DEFINE_TEST(test_misc_dynamic_prefix_sum_insert)
{
    typedef DynamicPrefixSumTree<unsigned, int, 3> TPrefixSumTree;
    typedef typename Value<TPrefixSumTree>::Type TNode;
    typedef typename TNode::TKeyTab   TKeyTable;
    typedef typename TNode::TCargoTab TCargoTable;
    typedef typename TNode::TChildTab TChildTable;

    TPrefixSumTree tree;

    SEQAN_ASSERT_EQ(empty(tree), true);

    insert(tree, 10, -2);

    SEQAN_ASSERT_EQ(empty(tree), false);

    {
        TKeyTable & keyTab = tree.rootPtr->_keyTab;
        TCargoTable & cargoTab = tree.rootPtr->_cargoTab;
        TChildTable & childTab = tree.rootPtr->_childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 1u);
        SEQAN_ASSERT_EQ(keyTab[0], 10u);

        SEQAN_ASSERT_EQ(length(cargoTab), 1u);
        SEQAN_ASSERT_EQ(cargoTab[0], -2);

        SEQAN_ASSERT_EQ(length(childTab), 0u);

        // Full Node.
        insert(tree, 20,  7);
        insert(tree,  5,  4);
        insert(tree, 25, -6);
        insert(tree,  0,  1);

        SEQAN_ASSERT_EQ(length(keyTab), 5u);
        SEQAN_ASSERT_EQ(keyTab[0], 0u);
        SEQAN_ASSERT_EQ(keyTab[1], 5u);
        SEQAN_ASSERT_EQ(keyTab[2], 10u);
        SEQAN_ASSERT_EQ(keyTab[3], 20u);
        SEQAN_ASSERT_EQ(keyTab[4], 25u);

        SEQAN_ASSERT_EQ(length(cargoTab), 5u);
        SEQAN_ASSERT_EQ(cargoTab[0], 1);    // 1
        SEQAN_ASSERT_EQ(cargoTab[1], 5);    // 4
        SEQAN_ASSERT_EQ(cargoTab[2], 3);    // -2
        SEQAN_ASSERT_EQ(cargoTab[3], 10);    // 7
        SEQAN_ASSERT_EQ(cargoTab[4], 4);    // -6

        SEQAN_ASSERT_EQ(length(childTab), 0u);
    }

    // Split root and insert in right child.
    insert(tree, 15, -3);

    {
        TKeyTable & keyTab = tree.rootPtr->_keyTab;
        TCargoTable & cargoTab = tree.rootPtr->_cargoTab;
        TChildTable & childTab = tree.rootPtr->_childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 1u);
        SEQAN_ASSERT_EQ(keyTab[0], 10u);

        SEQAN_ASSERT_EQ(length(cargoTab), 1u);
        SEQAN_ASSERT_EQ(cargoTab[0], 3);

        SEQAN_ASSERT_EQ(length(childTab), 2u);
    }

    // Check left child.
    {
        TNode& leftChild = getChild(*tree.rootPtr, 0);
        TKeyTable & keyTab = leftChild._keyTab;
        TCargoTable & cargoTab = leftChild._cargoTab;
        TChildTable & childTab = leftChild._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 2u);
        SEQAN_ASSERT_EQ(keyTab[0], 0u);
        SEQAN_ASSERT_EQ(keyTab[1], 5u);

        SEQAN_ASSERT_EQ(length(cargoTab), 2u);
        SEQAN_ASSERT_EQ(cargoTab[0], 1);
        SEQAN_ASSERT_EQ(cargoTab[1], 5);

        SEQAN_ASSERT_EQ(length(childTab), 0u);
    }

    // Check right child.
    {
        TNode& rightChild = getChild(*tree.rootPtr, 1);
        TKeyTable & keyTab = rightChild._keyTab;
        TCargoTable & cargoTab = rightChild._cargoTab;
        TChildTable & childTab = rightChild._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 3u);
        SEQAN_ASSERT_EQ(keyTab[0], 15u);
        SEQAN_ASSERT_EQ(keyTab[1], 20u);
        SEQAN_ASSERT_EQ(keyTab[2], 25u);

        SEQAN_ASSERT_EQ(length(cargoTab), 3u);
        SEQAN_ASSERT_EQ(cargoTab[0], -3);
        SEQAN_ASSERT_EQ(cargoTab[1], 4);
        SEQAN_ASSERT_EQ(cargoTab[2], -2);

        SEQAN_ASSERT_EQ(length(childTab), 0u);
    }

    // Fill left child.
    {
        insert(tree, 2, 2);
        insert(tree, 8, 5);
        insert(tree, 3, -8);

        TNode leftChild = getChild(*tree.rootPtr, 0);
        TKeyTable keyTab = leftChild._keyTab;
        TCargoTable cargoTab = leftChild._cargoTab;
        TChildTable childTab = leftChild._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 5u);
        SEQAN_ASSERT_EQ(keyTab[0], 0u);
        SEQAN_ASSERT_EQ(keyTab[1], 2u);
        SEQAN_ASSERT_EQ(keyTab[2], 3u);
        SEQAN_ASSERT_EQ(keyTab[3], 5u);
        SEQAN_ASSERT_EQ(keyTab[4], 8u);

        SEQAN_ASSERT_EQ(length(cargoTab), 5u);
        SEQAN_ASSERT_EQ(cargoTab[0], 1);  // 1
        SEQAN_ASSERT_EQ(cargoTab[1], 3);  // 2
        SEQAN_ASSERT_EQ(cargoTab[2], -5); // -8
        SEQAN_ASSERT_EQ(cargoTab[3], -1); // 4
        SEQAN_ASSERT_EQ(cargoTab[4], 4);  // 5

        SEQAN_ASSERT_EQ(length(childTab), 0u);

        TNode root = _getRoot(tree);

        keyTab = root._keyTab;
        cargoTab = root._cargoTab;
        childTab = root._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 1u);
        SEQAN_ASSERT_EQ(keyTab[0], 10u);

        SEQAN_ASSERT_EQ(length(cargoTab), 1u);
        SEQAN_ASSERT_EQ(cargoTab[0], 2);

        SEQAN_ASSERT_EQ(length(childTab), 2u);
    }

    // Split left child add to left half of split node.
    {
        insert(tree, 1, -100);

        TNode root = _getRoot(tree);

        TKeyTable   keyTab   = root._keyTab;
        TCargoTable cargoTab = root._cargoTab;
        TChildTable childTab = root._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 2u);
        SEQAN_ASSERT_EQ(keyTab[0], 3u);
        SEQAN_ASSERT_EQ(keyTab[1], 10u);

        SEQAN_ASSERT_EQ(length(cargoTab), 2u);
        SEQAN_ASSERT_EQ(cargoTab[0], -105);  // -8
        SEQAN_ASSERT_EQ(cargoTab[1], -98);   // -2

        SEQAN_ASSERT_EQ(length(childTab), 3u);

        TNode& leftChild = getChild(*tree.rootPtr, 0);
        keyTab = leftChild._keyTab;
        cargoTab = leftChild._cargoTab;
        childTab = leftChild._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 3u);
        SEQAN_ASSERT_EQ(keyTab[0], 0u);  // 1
        SEQAN_ASSERT_EQ(keyTab[1], 1u);  // -100
        SEQAN_ASSERT_EQ(keyTab[2], 2u);  // 2

        SEQAN_ASSERT_EQ(length(cargoTab), 3u);
        SEQAN_ASSERT_EQ(cargoTab[0], 1);    // 1
        SEQAN_ASSERT_EQ(cargoTab[1], -99);  // -100
        SEQAN_ASSERT_EQ(cargoTab[2], -97);  // 2

        SEQAN_ASSERT_EQ(length(childTab), 0u);

        TNode splitChild = getChild(*tree.rootPtr, 1);
        keyTab = splitChild._keyTab;
        cargoTab = splitChild._cargoTab;
        childTab = splitChild._childTab;

        SEQAN_ASSERT_EQ(keyTab[0], 5u);
        SEQAN_ASSERT_EQ(keyTab[1], 8u);

        SEQAN_ASSERT_EQ(cargoTab[0], 4);  // 4
        SEQAN_ASSERT_EQ(cargoTab[1], 9);  // 5

        SEQAN_ASSERT_EQ(length(childTab), 0u);
    }

    // Fill right child.
    {
        insert(tree, 50, -2);
        insert(tree, 11, -2);

        TNode rightChild = getChild(*tree.rootPtr, 2);
        TKeyTable keyTab = rightChild._keyTab;
        TCargoTable cargoTab = rightChild._cargoTab;
        TChildTable childTab = rightChild._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 5u);
        SEQAN_ASSERT_EQ(keyTab[0], 11u);  // -2
        SEQAN_ASSERT_EQ(keyTab[1], 15u);  // -3
        SEQAN_ASSERT_EQ(keyTab[2], 20u);  // 7
        SEQAN_ASSERT_EQ(keyTab[3], 25u);  // -6
        SEQAN_ASSERT_EQ(keyTab[4], 50u);  // -2

        SEQAN_ASSERT_EQ(length(cargoTab), 5u);
        SEQAN_ASSERT_EQ(cargoTab[0], -2);  // -2
        SEQAN_ASSERT_EQ(cargoTab[1], -5);  // -3
        SEQAN_ASSERT_EQ(cargoTab[2], 2); // 7
        SEQAN_ASSERT_EQ(cargoTab[3], -4); // -6
        SEQAN_ASSERT_EQ(cargoTab[4], -6);  // -2

        SEQAN_ASSERT_EQ(length(childTab), 0u);

        TNode root = _getRoot(tree);

        keyTab = root._keyTab;
        cargoTab = root._cargoTab;
        childTab = root._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 2u);
        SEQAN_ASSERT_EQ(keyTab[0], 3u);
        SEQAN_ASSERT_EQ(keyTab[1], 10u);

        SEQAN_ASSERT_EQ(length(cargoTab), 2u);
        SEQAN_ASSERT_EQ(cargoTab[0], -105);
        SEQAN_ASSERT_EQ(cargoTab[1], -98);

        SEQAN_ASSERT_EQ(length(childTab), 3u);
    }

    // Split right child and add to left split node.
    {
        insert(tree, 19, -2);

        TNode root = _getRoot(tree);

        TKeyTable keyTab = root._keyTab;
        TCargoTable cargoTab = root._cargoTab;
        TChildTable childTab = root._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 3u);
        SEQAN_ASSERT_EQ(keyTab[0], 3u);
        SEQAN_ASSERT_EQ(keyTab[1], 10u);
        SEQAN_ASSERT_EQ(keyTab[2], 20u);

        SEQAN_ASSERT_EQ(length(cargoTab), 3u);
        SEQAN_ASSERT_EQ(cargoTab[0], -105);
        SEQAN_ASSERT_EQ(cargoTab[1], -98);
        SEQAN_ASSERT_EQ(cargoTab[2], -98);

        SEQAN_ASSERT_EQ(length(childTab), 4u);

        TNode rightChild = getChild(*tree.rootPtr, 2);
        keyTab = rightChild._keyTab;
        cargoTab = rightChild._cargoTab;
        childTab = rightChild._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 3u);
        SEQAN_ASSERT_EQ(keyTab[0], 11u);  // -2
        SEQAN_ASSERT_EQ(keyTab[1], 15u);  // -3
        SEQAN_ASSERT_EQ(keyTab[2], 19u);  // -3

        SEQAN_ASSERT_EQ(length(cargoTab), 3u);
        SEQAN_ASSERT_EQ(cargoTab[0], -2);  // -2
        SEQAN_ASSERT_EQ(cargoTab[1], -5);  // -3
        SEQAN_ASSERT_EQ(cargoTab[2], -7);  // -2

        SEQAN_ASSERT_EQ(length(childTab), 0u);

        TNode splitChild = getChild(*tree.rootPtr, 3);
        keyTab = splitChild._keyTab;
        cargoTab = splitChild._cargoTab;
        childTab = splitChild._childTab;

        SEQAN_ASSERT_EQ(keyTab[0], 25u);  // -6
        SEQAN_ASSERT_EQ(keyTab[1], 50u);  // -2

        SEQAN_ASSERT_EQ(cargoTab[0], -6);  // -6
        SEQAN_ASSERT_EQ(cargoTab[1], -8);  // -2

        SEQAN_ASSERT_EQ(length(childTab), 0u);
    }

    // Split root.
    {
        insert(tree, 100, 10);
        insert(tree, 12, 10);
        insert(tree, 150, 10);
        insert(tree, 16, 10);
        insert(tree, 200, 10);
        insert(tree, 17, 10);
        insert(tree, 250, 10);

        TNode root = _getRoot(tree);
        TKeyTable keyTab = root._keyTab;
        TCargoTable cargoTab = root._cargoTab;
        TChildTable childTab = root._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 5u);
        SEQAN_ASSERT_EQ(keyTab[0], 3u);
        SEQAN_ASSERT_EQ(keyTab[1], 10u);
        SEQAN_ASSERT_EQ(keyTab[2], 15u);
        SEQAN_ASSERT_EQ(keyTab[3], 20u);
        SEQAN_ASSERT_EQ(keyTab[4], 100u);

        SEQAN_ASSERT_EQ(length(cargoTab), 5u);
        SEQAN_ASSERT_EQ(cargoTab[0], -105);
        SEQAN_ASSERT_EQ(cargoTab[1], -98);
        SEQAN_ASSERT_EQ(cargoTab[2], -93);
        SEQAN_ASSERT_EQ(cargoTab[3], -68);
        SEQAN_ASSERT_EQ(cargoTab[4], -66);

        SEQAN_ASSERT_EQ(length(childTab), 6u);

        // Actually split the root.
        insert(tree, 7, 10);
        insert(tree, 9, 10);
        insert(tree, 4, 10);
        insert(tree, 6, 10);

        root = _getRoot(tree);
        keyTab = root._keyTab;
        cargoTab = root._cargoTab;
        childTab = root._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 1u);
        SEQAN_ASSERT_EQ(keyTab[0], 15u);

        SEQAN_ASSERT_EQ(length(cargoTab), 1u);
        SEQAN_ASSERT_EQ(cargoTab[0], -53);

        SEQAN_ASSERT_EQ(length(childTab), 2u);

        // New left child < 15

        TNode leftChildL1 = getChild(_getRoot(tree), 0);
        keyTab = leftChildL1._keyTab;
        cargoTab = leftChildL1._cargoTab;
        childTab = leftChildL1._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 3u);
        SEQAN_ASSERT_EQ(keyTab[0], 3u);
        SEQAN_ASSERT_EQ(keyTab[1], 7u);
        SEQAN_ASSERT_EQ(keyTab[2], 10u);
        SEQAN_ASSERT_EQ(length(cargoTab), 3u);
        SEQAN_ASSERT_EQ(cargoTab[0], -105);
        SEQAN_ASSERT_EQ(cargoTab[1], -71);
        SEQAN_ASSERT_EQ(cargoTab[2], -58);
        SEQAN_ASSERT_EQ(length(childTab), 4u);

        // Children of leftChild1
        {  // Leaf 0
            TNode l1Leaf0 = getChild(leftChildL1, 0);
            keyTab = l1Leaf0._keyTab;
            cargoTab = l1Leaf0._cargoTab;
            childTab = l1Leaf0._childTab;

            SEQAN_ASSERT_EQ(length(keyTab), 3u);
            SEQAN_ASSERT_EQ(keyTab[0], 0u);
            SEQAN_ASSERT_EQ(keyTab[1], 1u);
            SEQAN_ASSERT_EQ(keyTab[2], 2u);
            SEQAN_ASSERT_EQ(length(cargoTab), 3u);
            SEQAN_ASSERT_EQ(cargoTab[0], 1);
            SEQAN_ASSERT_EQ(cargoTab[1], -99);
            SEQAN_ASSERT_EQ(cargoTab[2], -97);
            SEQAN_ASSERT_EQ(length(childTab), 0u);

            // Leaf 1
            TNode l1Leaf1 = getChild(leftChildL1, 1);
            keyTab = l1Leaf1._keyTab;
            cargoTab = l1Leaf1._cargoTab;
            childTab = l1Leaf1._childTab;

            SEQAN_ASSERT_EQ(length(keyTab), 3u);
            SEQAN_ASSERT_EQ(keyTab[0], 4u);
            SEQAN_ASSERT_EQ(keyTab[1], 5u);
            SEQAN_ASSERT_EQ(keyTab[2], 6u);
            SEQAN_ASSERT_EQ(length(cargoTab), 3u);
            SEQAN_ASSERT_EQ(cargoTab[0], 10);
            SEQAN_ASSERT_EQ(cargoTab[1], 14);
            SEQAN_ASSERT_EQ(cargoTab[2], 24);
            SEQAN_ASSERT_EQ(length(childTab), 0u);

            //Leaf 2
            TNode l1Leaf2 = getChild(leftChildL1, 2);
            keyTab = l1Leaf2._keyTab;
            cargoTab = l1Leaf2._cargoTab;
            childTab = l1Leaf2._childTab;

            SEQAN_ASSERT_EQ(length(keyTab), 2u);
            SEQAN_ASSERT_EQ(keyTab[0], 8u);
            SEQAN_ASSERT_EQ(keyTab[1], 9u);
            SEQAN_ASSERT_EQ(length(cargoTab), 2u);
            SEQAN_ASSERT_EQ(cargoTab[0], 5);
            SEQAN_ASSERT_EQ(cargoTab[1], 15);
            SEQAN_ASSERT_EQ(length(childTab), 0u);

            //Leaf 3
            TNode l1Leaf3 = getChild(leftChildL1, 3);
            keyTab = l1Leaf3._keyTab;
            cargoTab = l1Leaf3._cargoTab;
            childTab = l1Leaf3._childTab;

            SEQAN_ASSERT_EQ(length(keyTab), 2u);
            SEQAN_ASSERT_EQ(keyTab[0], 11u);
            SEQAN_ASSERT_EQ(keyTab[1], 12u);
            SEQAN_ASSERT_EQ(length(cargoTab), 2u);
            SEQAN_ASSERT_EQ(cargoTab[0], -2);
            SEQAN_ASSERT_EQ(cargoTab[1], 8);
            SEQAN_ASSERT_EQ(length(childTab), 0u);
        }

        // New right child > 15
        TNode rightChildL1 = getChild(_getRoot(tree), 1);
        keyTab = rightChildL1._keyTab;
        cargoTab = rightChildL1._cargoTab;
        childTab = rightChildL1._childTab;

        SEQAN_ASSERT_EQ(length(keyTab), 2u);
        SEQAN_ASSERT_EQ(keyTab[0], 20u);
        SEQAN_ASSERT_EQ(keyTab[1], 100u);
        SEQAN_ASSERT_EQ(length(cargoTab), 2u);
        SEQAN_ASSERT_EQ(cargoTab[0], 25);
        SEQAN_ASSERT_EQ(cargoTab[1], 27);
        SEQAN_ASSERT_EQ(length(childTab), 3u);

        // Children of rightChildL1
        {  // Leaf 0
            TNode r1Leaf0 = getChild(rightChildL1, 0);
            keyTab = r1Leaf0._keyTab;
            cargoTab = r1Leaf0._cargoTab;
            childTab = r1Leaf0._childTab;

            SEQAN_ASSERT_EQ(length(keyTab), 3u);
            SEQAN_ASSERT_EQ(keyTab[0], 16u);
            SEQAN_ASSERT_EQ(keyTab[1], 17u);
            SEQAN_ASSERT_EQ(keyTab[2], 19u);
            SEQAN_ASSERT_EQ(length(cargoTab), 3u);
            SEQAN_ASSERT_EQ(cargoTab[0], 10);
            SEQAN_ASSERT_EQ(cargoTab[1], 20);
            SEQAN_ASSERT_EQ(cargoTab[2], 18);
            SEQAN_ASSERT_EQ(length(childTab), 0u);

            // Leaf 1
            TNode r1Leaf1 = getChild(rightChildL1, 1);
            keyTab = r1Leaf1._keyTab;
            cargoTab = r1Leaf1._cargoTab;
            childTab = r1Leaf1._childTab;

            SEQAN_ASSERT_EQ(length(keyTab), 2u);
            SEQAN_ASSERT_EQ(keyTab[0], 25u);
            SEQAN_ASSERT_EQ(keyTab[1], 50u);
            SEQAN_ASSERT_EQ(length(cargoTab), 2u);
            SEQAN_ASSERT_EQ(cargoTab[0], -6);
            SEQAN_ASSERT_EQ(cargoTab[1], -8);
            SEQAN_ASSERT_EQ(length(childTab), 0u);

            //Leaf 2
            TNode r1Leaf2 = getChild(rightChildL1, 2);
            keyTab = r1Leaf2._keyTab;
            cargoTab = r1Leaf2._cargoTab;
            childTab = r1Leaf2._childTab;

            SEQAN_ASSERT_EQ(length(keyTab), 3u);
            SEQAN_ASSERT_EQ(keyTab[0], 150u);
            SEQAN_ASSERT_EQ(keyTab[1], 200u);
            SEQAN_ASSERT_EQ(keyTab[2], 250u);
            SEQAN_ASSERT_EQ(length(cargoTab), 3u);
            SEQAN_ASSERT_EQ(cargoTab[0], 10);
            SEQAN_ASSERT_EQ(cargoTab[1], 20);
            SEQAN_ASSERT_EQ(cargoTab[2], 30);
            SEQAN_ASSERT_EQ(length(childTab), 0u);
        }
    }
}


struct TestCompareLess_
{
    template <typename T>
    inline bool
    operator()(T const & i, T const & j)
    {
        return i.i1 < j.i1;
    }
};

template <typename TFreq, typename TPos>
inline int
_testComputePrefixSum(TFreq const & data, TPos const & pos)
{
    typedef typename Iterator<TFreq const, Standard>::Type TIterator;

    int sum = 0;
    Pair<unsigned, int> key(pos, 0);
    TIterator itEnd = std::upper_bound(begin(data, Standard()), end(data, Standard()), key, TestCompareLess_());
    for (TIterator it = begin(data, Standard()); it != itEnd; ++it)
        sum += it->i2;
    return sum;
}

SEQAN_DEFINE_TEST(test_misc_dynamic_prefix_sum_prefix_sum)
{
    typedef DynamicPrefixSumTree<unsigned, int, 64> TPrefixSumTree;

    String<Pair<unsigned, int> > testData;
    TPrefixSumTree tree;

    Rng<MersenneTwister> rng(42);
    Pdf<Uniform<int> > freqPdf(-10000, 10000);
    Pdf<Uniform<unsigned> > posPdf(0, 1000000);
    String<bool, Packed<> > bitVec;
    resize(bitVec, 1000000, false, Exact());

    const unsigned LENGTH = 100000;
    resize(testData, LENGTH, Exact());
    for (unsigned i = 0; i < LENGTH; ++i)
    {
        testData[i].i1 = pickRandomNumber(rng, posPdf);
        if (bitVec[testData[i].i1])  // Do not pick the same number twice.
        {
            --i;
            continue;
        }
        bitVec[testData[i].i1] = true;
        testData[i].i2 = pickRandomNumber(rng, freqPdf);
        insert(tree, testData[i].i1, testData[i].i2);
    }

    std::sort(begin(testData, Standard()), end(testData, Standard()));
    int testSum = 0;
    for (unsigned i = 0; i < LENGTH; ++i)
    {
        testSum += testData[i].i2;
        SEQAN_ASSERT_EQ(prefixSum(tree, testData[i].i1), testSum);
    }

    // Find random positions.
    for (unsigned i = 0; i < LENGTH; ++i)
    {
        unsigned pos = pickRandomNumber(rng, posPdf);
        SEQAN_ASSERT_EQ(prefixSum(tree, pos), _testComputePrefixSum(testData, pos));
    }
}


#endif  // EXTRAS_TESTS_MISC_DYNAMIC_PREFIX_SUM_TEST_MISC_DYNAMIC_PREFIX_SUM_H_
