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

using namespace seqan;

SEQAN_DEFINE_TEST(test_graph_types_derived_oracle)
{
    Graph<Automaton<char> > g;
    createOracleOnReverse(g, "announce");
    SEQAN_ASSERT(parseString(g, 0, "e") == 1);
    SEQAN_ASSERT(parseString(g, 0, "ec") == 2);
    SEQAN_ASSERT(parseString(g, 0, "n") == 3);
    SEQAN_ASSERT(parseString(g, 0, "a") == 8);
    SEQAN_ASSERT(parseString(g, 0, "nn") == 7);

    Graph<Automaton<Dna> > g2;
    createOracle(g2, "ATATA");
    SEQAN_ASSERT(parseString(g2, 0, "A") == 1);
    SEQAN_ASSERT(parseString(g2, 0, "T") == 2);
    SEQAN_ASSERT(parseString(g2, 0, "AT") == 2);
}

SEQAN_DEFINE_TEST(test_graph_types_derived_trie)
{
    typedef Position<String<char> >::Type TPosition;

    Graph<Automaton<char> > g;
    String<String<TPosition> > pos;
    String<String<char> > keywords;
    appendValue(keywords, String<char>("announce"));
    appendValue(keywords, String<char>("annual"));
    appendValue(keywords, String<char>("annually"));
    createTrie(g, pos, keywords);

    SEQAN_ASSERT(parseString(g, 0, "a") == 1);
    SEQAN_ASSERT(parseString(g, 0, "an") == 2);
    SEQAN_ASSERT(parseString(g, 0, "ann") == 3);
    SEQAN_ASSERT(parseString(g, 0, "anno") == 4);
    SEQAN_ASSERT(parseString(g, 0, "annu") == 9);
    SEQAN_ASSERT(getProperty(pos, 11) == (unsigned int) 1); // In vertex 11 keyword 1 ends
    SEQAN_ASSERT(getProperty(pos, 13) == (unsigned int) 2);
    SEQAN_ASSERT(getProperty(pos, 8) == (unsigned int) 0);

    clear(g);
    clear(pos);
    createTrieOnReverse(g, pos, keywords);

    SEQAN_ASSERT(parseString(g, 0, "e") == 1);
    SEQAN_ASSERT(parseString(g, 0, "l") == 9);
    SEQAN_ASSERT(parseString(g, 0, "y") == 15);
    SEQAN_ASSERT(parseString(g, 0, "ec") == 2);
    SEQAN_ASSERT(getProperty(pos, 8) == (unsigned int) 0); // In vertex 8 keyword 0 ends
    SEQAN_ASSERT(getProperty(pos, 14) == (unsigned int) 1);
    SEQAN_ASSERT(getProperty(pos, 22) == (unsigned int) 2);

    Graph<Automaton<Dna> > gDna;
    clear(pos);
    String<String<Dna> > keyw;
    appendValue(keyw, String<Dna>("ATATATA"));
    appendValue(keyw, String<Dna>("TATAT"));
    appendValue(keyw, String<Dna>("ACGATAT"));
    createTrie(gDna, pos, keyw);

    SEQAN_ASSERT(parseString(gDna, 0, "A") == 1);
    SEQAN_ASSERT(parseString(gDna, 0, "T") == 8);
    SEQAN_ASSERT(parseString(gDna, 0, "AT") == 2);
    SEQAN_ASSERT(parseString(gDna, 0, "AC") == 13);
    SEQAN_ASSERT(getProperty(pos, 7) == (unsigned int) 0); // In vertex 7 keyword 0 ends
    SEQAN_ASSERT(getProperty(pos, 18) == (unsigned int) 2);
    SEQAN_ASSERT(getProperty(pos, 12) == (unsigned int) 1);


    //createSuffixTrie
    clear(g);
    clear(pos);
    char * str = (char *) "ABABBA";
    char * strend = end(str);
    char * it;
    typedef VertexDescriptor<Graph<Automaton<char> > >::Type TVertexDescriptor;
    TVertexDescriptor v;

    createSuffixTrie(g, pos, str);
    for (unsigned int i = 0; i < length(str); ++i)
    {
        v = parseString(g, 0, it = str + i, strend);
        SEQAN_ASSERT(it == strend);
        SEQAN_ASSERT(getProperty(pos, v) == i);
    }

    SEQAN_ASSERT(canParseString(g, "ABBA"));
    SEQAN_ASSERT(canParseString(g, "BAB"));

    SEQAN_ASSERT(!canParseString(g, "AA"));
    SEQAN_ASSERT(!canParseString(g, "BAC"));
    SEQAN_ASSERT(!canParseString(g, "C"));
    SEQAN_ASSERT(canParseString(g, ""));
}

SEQAN_DEFINE_TEST(test_graph_types_derived_set_oracle)
{
    typedef Position<String<char> >::Type TPosition;

    Graph<Automaton<char> > g;
    String<String<TPosition> > pos;
    String<String<char> > keywords;
    appendValue(keywords, String<char>("announce"));
    appendValue(keywords, String<char>("annual"));
    appendValue(keywords, String<char>("annually"));
    createSetOracle(g, pos, keywords);

    for (unsigned int i = 0; i < length(keywords); ++i)
    {
        String<char> & str = keywords[i];
        for (unsigned int j = 0; j < length(str); ++j)
        {
            SEQAN_ASSERT(canParseString(g, prefix(str, i)));
        }
    }

    SEQAN_ASSERT(!canParseString(g, "d"));
    SEQAN_ASSERT(!canParseString(g, "annly"));
    SEQAN_ASSERT(canParseString(g, ""));
}

SEQAN_BEGIN_TESTSUITE(test_graph_types_derived)
{
    SEQAN_CALL_TEST(test_graph_types_derived_set_oracle);
    SEQAN_CALL_TEST(test_graph_types_derived_trie);
    SEQAN_CALL_TEST(test_graph_types_derived_oracle);
}
SEQAN_END_TESTSUITE
