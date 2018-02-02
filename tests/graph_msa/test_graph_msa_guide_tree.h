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
// Author: Tobias Rausch <rausch@embl.de>
// ==========================================================================
// Tests for the Graph MSA module.
// ==========================================================================

#ifndef TESTS_TEST_GRAPH_MSA_GUIDE_TREE_H_
#define TESTS_TEST_GRAPH_MSA_GUIDE_TREE_H_

void Test_GuideTree_NeighbourJoining()
{
    using namespace seqan;
//____________________________________________________________________________
// Neighbor Joining

    // Create a distance matrix
    String<double> mat;
    resize(mat, 8*8, 0);
    assignValue(mat, 0*8+1, 7);assignValue(mat, 0*8+2, 8);assignValue(mat, 0*8+3, 11);assignValue(mat, 0*8+4, 13);assignValue(mat, 0*8+5, 16);assignValue(mat, 0*8+6, 13);assignValue(mat, 0*8+7, 17);
    assignValue(mat, 1*8+2, 5);assignValue(mat, 1*8+3, 8);assignValue(mat, 1*8+4, 10);assignValue(mat, 1*8+5, 13);assignValue(mat, 1*8+6, 10);assignValue(mat, 1*8+7, 14);
    assignValue(mat, 2*8+3, 5);assignValue(mat, 2*8+4, 7);assignValue(mat, 2*8+5, 10);assignValue(mat, 2*8+6, 7);assignValue(mat, 2*8+7, 11);
    assignValue(mat, 3*8+4, 8);assignValue(mat, 3*8+5, 11);assignValue(mat, 3*8+6, 8);assignValue(mat, 3*8+7, 12);
    assignValue(mat, 4*8+5, 5);assignValue(mat, 4*8+6, 6);assignValue(mat, 4*8+7, 10);
    assignValue(mat, 5*8+6, 9);assignValue(mat, 5*8+7, 13);
    assignValue(mat, 6*8+7, 8);

    typedef Graph<Tree<double> > TGraph;
    TGraph guideTreeOut;
    njTree(mat, guideTreeOut);
    //std::cout << guideTreeOut << std::endl;

    SEQAN_ASSERT(numVertices(guideTreeOut) == 15);
    SEQAN_ASSERT(findEdge(guideTreeOut, 8, 1) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 8, 0) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 9, 5) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 9, 4) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 10, 2) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 10, 8) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 10, 2) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 10, 8) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 11, 3) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 11, 10) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 12, 9) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 12, 11) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 13, 12) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 13, 6) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 14, 13) != 0);
    SEQAN_ASSERT(findEdge(guideTreeOut, 14, 7) != 0);
    SEQAN_ASSERT(getRoot(guideTreeOut) == 14);
}

template<typename TTag>
void
Test_UpgmaGuideTree(int seed) {
    using namespace seqan;

    typedef unsigned int TSize;

    std::mt19937 rng(seed);

    for(TSize i = 0; i < 10; ++i) {
        // Set-up a sparse distance matrix
        std::uniform_int_distribution<TSize> pdfN(2, 11);
        TSize n = pdfN(rng);
        Graph<Undirected<double> > distGraph;
        String<double> distMatrix;
        resize(distMatrix, n * n, 0);
        TSize all = (n * (n - 1)) / 2;
        typedef std::set<double> TDistanceSet;
        //typedef TDistanceSet::iterator TSetIt;
        TDistanceSet distances;
        String<double> myDist;
        typedef Iterator<String<double> >::Type TStringIter;
        while (distances.size() < all) {
            std::uniform_real_distribution<double> pdf(0, 1000000);
            double newVal = pdf(rng);
            if (distances.insert(newVal).second) {
                appendValue(myDist, newVal);
            }
        }
        double infCargo = _getInfinity<double>();
        //clear(myDist); appendValue(myDist, infCargo); appendValue(myDist, infCargo); appendValue(myDist, 84);
        TStringIter strIt = begin(myDist);
        for(TSize row = 0; row < n; ++row)
            addVertex(distGraph);
        for(TSize row = 0; row < n; ++row) {
            for(TSize col = n - 1; col > row; --col) {
                addEdge(distGraph, row, col, value(strIt));
                value(distMatrix, row * n + col) = value(strIt);
                goNext(strIt);
            }
        }
        //removeEdge(distGraph, 0, 1);removeEdge(distGraph, 0, 2);
        Graph<Undirected<double> > distGraphCopy;
        String<double> distMatrixCopy;
        for(TSize row = 0; row < n; ++row) {
            for(TSize col = n - 1; col > row; --col) {
                distGraphCopy = distGraph;
                distMatrixCopy = distMatrix;

                std::uniform_real_distribution<double> pdf(0, 1.0);
                if (pdf(rng) < 0.5) {
                    value(distMatrix, row * n + col) = infCargo;
                    removeEdge(distGraph, row, col);
                }

                String<size_t> _;
                if (connectedComponents(_, distGraph)) {
                    move(distGraph, distGraphCopy);
                    move(distMatrix, distMatrixCopy);
                }
            }
        }
        // Guide Tree
        Graph<Tree<double> > guideTreeGraph;
        upgmaTree(distGraph, guideTreeGraph, TTag());
        Graph<Tree<double> > guideTreeMat;
        upgmaTree(distMatrix, guideTreeMat, TTag());
        typedef Iterator<Graph<Tree<double> >, BfsIterator>::Type TBfsIterator;
        String<TSize> set1;
        TBfsIterator itBfs(guideTreeGraph, getRoot(guideTreeGraph));
        for(;!atEnd(itBfs);goNext(itBfs)) appendValue(set1, value(itBfs));
        String<TSize> set2;
        TBfsIterator itBfs2(guideTreeMat, getRoot(guideTreeMat));
        for(;!atEnd(itBfs2);goNext(itBfs2)) appendValue(set2, value(itBfs2));
        SEQAN_ASSERT(set1 == set1);
        /*
        if (set1 != set2) {
            std::cout << "Randomized test failed:" << std::endl;
            std::cout << "Upgma Guide Trees:" << std::endl;
            std::cout << guideTreeMat << std::endl;
            std::cout << guideTreeGraph << std::endl;
            for(TSize i=0;i<n;++i) {
                for(TSize j=i+1;j<n;++j) {
                    std::cout << value(distMatrix, i*n+j) << ",";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
            exit(0);
        }
        */
    }
}

SEQAN_DEFINE_TEST(test_graph_msa_guide_tree_neighbour_joining)
{
    Test_GuideTree_NeighbourJoining();
}

SEQAN_DEFINE_TEST(test_graph_msa_guide_tree_upgma_weight_avg)
{
    for (int i = 0; i < 10; ++i)
        Test_UpgmaGuideTree<seqan::UpgmaWeightAvg>(i);
}

SEQAN_DEFINE_TEST(test_graph_msa_guide_tree_upgma_avg)
{
    for (int i = 0; i < 10; ++i)
        Test_UpgmaGuideTree<seqan::UpgmaAvg>(i);
}

SEQAN_DEFINE_TEST(test_graph_msa_guide_tree_upgma_min)
{
    for (int i = 0; i < 10; ++i)
        Test_UpgmaGuideTree<seqan::UpgmaMin>(i);
}

SEQAN_DEFINE_TEST(test_graph_msa_guide_tree_upgma_max)
{
    for (int i = 0; i < 10; ++i)
        Test_UpgmaGuideTree<seqan::UpgmaMax>(i);
}

#endif  // #ifndef TESTS_TEST_GRAPH_MSA_GUIDE_TREE_H_
