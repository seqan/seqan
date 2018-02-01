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

#include <seqan/basic.h>
#include <seqan/graph_types.h>

// Test external property maps (using String) with standard interface
SEQAN_DEFINE_TEST(test_graph_types_property_map_external)
{
    using namespace seqan;

    // Create graph typedefs.
    typedef Graph<Directed<void> >                  TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type   TEdgeDescriptor;

    // Initialize graph.
    TGraph g;
    TVertexDescriptor v0 = addVertex(g);
    TVertexDescriptor v1 = addVertex(g);
    addVertex(g);
    TEdgeDescriptor e2 = addEdge(g, 0, 2);
    TEdgeDescriptor e1 = addEdge(g, v0, v1);

    // Test generic/external property map interface.
    String<int> dMap;
    resizeVertexMap(dMap, g);

    String<char> eMap;
    resizeEdgeMap(eMap, g);

    assignProperty(dMap, v0, 3);
    assignProperty(dMap, v1, 1);
    assignProperty(eMap, e1, 'a');
    assignProperty(eMap, e2, 'b');
    SEQAN_ASSERT(getProperty(dMap, v0) == 3);
    SEQAN_ASSERT(getProperty(dMap, v1) == 1);
    SEQAN_ASSERT(getProperty(eMap, e2) == 'b');
    SEQAN_ASSERT(getProperty(eMap, e1) == 'a');

    property(dMap, v1) = 2;
    property(eMap, e2) = 'c';
    SEQAN_ASSERT(getProperty(dMap, v1) == 2);
    SEQAN_ASSERT(getProperty(eMap, e2) == 'c');

    String<int> const dMap2(dMap);
    SEQAN_ASSERT(getProperty(dMap2, v0) == 3);
    SEQAN_ASSERT(getProperty(dMap2, v1) == 2);
    SEQAN_ASSERT(property(dMap2, v1) == 2);
}

// Test external property maps (using String) with assignVertexMap().
SEQAN_DEFINE_TEST(test_graph_types_property_map_external_assign_vertex_map)
{
    using namespace seqan;

    // Create graph typedefs.
    typedef Graph<Directed<void> >                  TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    // Initialize graph.
    TGraph g;
    TVertexDescriptor v0 = addVertex(g);
    TVertexDescriptor v1 = addVertex(g);
    TVertexDescriptor v2 = addVertex(g);

    char names[] = {
    'r', 's', 't'
    };
    String<char> nameMap;
    assignVertexMap(nameMap, g, names);
    SEQAN_ASSERT(getProperty(nameMap, v0) == 'r');
    SEQAN_ASSERT(getProperty(nameMap, v1) == 's');
    SEQAN_ASSERT(getProperty(nameMap, v2) == 't');
}

// Test external property maps (using String) with assignEdgeMap().
SEQAN_DEFINE_TEST(test_graph_types_property_map_external_assign_edge_map)
{
    using namespace seqan;

    // Create graph typedefs.
    typedef Graph<Directed<void> >                  TGraph;

    // Initialize graph.
    TGraph g10;

    unsigned int weights[] = {
    4, 8
    };
    addVertex(g10);
    addVertex(g10);
    addVertex(g10);
    addEdge(g10, 0, 1);
    addEdge(g10, 0, 2);

    String<int> weightMap;
    assignEdgeMap(weightMap, g10, weights);
    SEQAN_ASSERT(getProperty(weightMap, findEdge(g10, 0, 1)) == 4);
    SEQAN_ASSERT(getProperty(weightMap, findEdge(g10, 0, 2)) == 8);
}

// Test internal pointer property maps.
SEQAN_DEFINE_TEST(test_graph_types_property_internal_pointer_map)
{
    using namespace seqan;

    // Graph typedefs.
    typedef Pair<char, int>                TPair;
    typedef Directed<TPair>                TEdges;
    typedef Graph<TEdges>                  TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type   TEdgeDescriptor;

    // Build graph.
    TGraph g;
    TVertexDescriptor v0 = addVertex(g);
    TVertexDescriptor v1 = addVertex(g);
    addVertex(g);
    TEdgeDescriptor e1 = addEdge(g, 0, 2);
    TEdgeDescriptor e2 = addEdge(g, v0, v1);

    // Second Variant: Pointer to member using a class
    InternalPointerPropertyMap<char TPair::*, & TPair::i1> eMap4;
    InternalPointerPropertyMap<int TPair::*, & TPair::i2> eMap5;
    resizeEdgeMap(eMap4, g);
    resizeEdgeMap(eMap5, g);
    assignProperty(eMap4, e1, 'c');
    assignProperty(eMap5, e1, 10);
    assignProperty(eMap4, e2, 'd');
    assignProperty(eMap5, e2, 30);
    SEQAN_ASSERT(getProperty(eMap4, e1) == 'c');
    SEQAN_ASSERT(getProperty(eMap5, e1) == 10);
    SEQAN_ASSERT(getProperty(eMap4, e2) == 'd');
    SEQAN_ASSERT(getProperty(eMap5, e2) == 30);
    property(eMap4, e1) = 'z';
    property(eMap5, e1) = 100;
    SEQAN_ASSERT(getProperty(eMap4, e1) == 'z');
    SEQAN_ASSERT(getProperty(eMap5, e1) == 100);
    InternalPointerPropertyMap<char TPair::*, & TPair::i1> const eMap6(eMap4);
    SEQAN_ASSERT(getProperty(eMap6, e1) == 'z');
    SEQAN_ASSERT(getProperty(eMap6, e2) == 'd');
    SEQAN_ASSERT(property(eMap6, e2) == 'd');
}

SEQAN_DEFINE_TEST(test_graph_types_property_internal_pointer_map_assign_edge_map_member)
{
    using namespace seqan;

    // Create graph typedefs.
    typedef Pair<unsigned, unsigned>                TPair;
    typedef Graph<Directed<TPair> >                 TGraph;

    // Initialize graph.
    TGraph g10;
    unsigned int weights1[] = {4, 8};
    unsigned int weights2[] = {5, 9};

    addVertex(g10);
    addVertex(g10);
    addVertex(g10);

    addEdge(g10, 0, 1);
    addEdge(g10, 0, 2);

    InternalPointerPropertyMap<unsigned TPair::*, & TPair::i1> weightMap1;
    InternalPointerPropertyMap<unsigned TPair::*, & TPair::i2> weightMap2;
    assignEdgeMap(weightMap1, g10, weights1);
    assignEdgeMap(weightMap2, g10, weights2);
    SEQAN_ASSERT_EQ(getProperty(weightMap1, findEdge(g10, 0, 1)), 4u);
    SEQAN_ASSERT_EQ(getProperty(weightMap1, findEdge(g10, 0, 2)), 8u);
    SEQAN_ASSERT_EQ(getProperty(weightMap2, findEdge(g10, 0, 1)), 5u);
    SEQAN_ASSERT_EQ(getProperty(weightMap2, findEdge(g10, 0, 2)), 9u);
}

// Test internal property maps.
SEQAN_DEFINE_TEST(test_graph_types_property_internal_map)
{
    using namespace seqan;

    // Graph typedefs.
    typedef Directed<char>                 TEdges;
    typedef Graph<TEdges>                  TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type   TEdgeDescriptor;

    // Build graph.
    TGraph g;
    TVertexDescriptor v0 = addVertex(g);
    TVertexDescriptor v1 = addVertex(g);
    addVertex(g);
    TEdgeDescriptor e1 = addEdge(g, 0, 2);
    TEdgeDescriptor e2 = addEdge(g, v0, v1);

    // Second Variant: Pointer to member using a class
    InternalPropertyMap<char> eMap4;
    resizeEdgeMap(eMap4, g);
    assignProperty(eMap4, e1, 'c');
    assignProperty(eMap4, e2, 'd');
    SEQAN_ASSERT(getProperty(eMap4, e1) == 'c');
    SEQAN_ASSERT(getProperty(eMap4, e2) == 'd');
    property(eMap4, e1) = 'z';
    SEQAN_ASSERT(getProperty(eMap4, e1) == 'z');
    InternalPropertyMap<char> const eMap6(eMap4);
    SEQAN_ASSERT(getProperty(eMap6, e1) == 'z');
    SEQAN_ASSERT(getProperty(eMap6, e2) == 'd');
    SEQAN_ASSERT(property(eMap6, e2) == 'd');
}

SEQAN_DEFINE_TEST(test_graph_types_property_internal_map_assign_edge_map)
{
    using namespace seqan;

    // Create graph typedefs.
    typedef Graph<Directed<unsigned> >                 TGraph;

    // Initialize graph.
    TGraph g10;
    unsigned int weights1[] = {4, 8};

    addVertex(g10);
    addVertex(g10);
    addVertex(g10);

    addEdge(g10, 0, 1);
    addEdge(g10, 0, 2);

    InternalPropertyMap<unsigned> weightMap1;
    assignEdgeMap(weightMap1, g10, weights1);
    SEQAN_ASSERT_EQ(getProperty(weightMap1, findEdge(g10, 0, 1)), 4u);
    SEQAN_ASSERT_EQ(getProperty(weightMap1, findEdge(g10, 0, 2)), 8u);
}

SEQAN_BEGIN_TESTSUITE(test_graph_types_property_map)
{
    SEQAN_CALL_TEST(test_graph_types_property_map_external);
    SEQAN_CALL_TEST(test_graph_types_property_map_external_assign_vertex_map);
    SEQAN_CALL_TEST(test_graph_types_property_map_external_assign_edge_map);

    SEQAN_CALL_TEST(test_graph_types_property_internal_pointer_map);
    SEQAN_CALL_TEST(test_graph_types_property_internal_pointer_map_assign_edge_map_member);

    SEQAN_CALL_TEST(test_graph_types_property_internal_map);
    SEQAN_CALL_TEST(test_graph_types_property_internal_map_assign_edge_map);
}
SEQAN_END_TESTSUITE
