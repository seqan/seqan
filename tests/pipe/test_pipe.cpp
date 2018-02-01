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

#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>

#include <seqan/stream.h>
#include <seqan/pipe.h>
#include <seqan/system.h>
#include <seqan/index.h>

#include "test_pipe.h"

using namespace std;
using namespace seqan;


// Maximum size of data to use for the external tests.
const size_t MAX_SIZE = 1u << 20;


// TODO(holtgrew): The following test* functions should actually be defined with SEQAN_DEFINE_TEST().


template <typename TStringSpec>
void testExternalString(unsigned maxSize = 32*1024*1024)
{
	typedef String<unsigned, TStringSpec> TExtString;
    typedef typename Iterator<TExtString const, Standard>::Type TIter;

    Buffer<unsigned> buf;
    allocPage(buf, maxSize, buf);

    TExtString extString;
    for(unsigned i = 1; i <= maxSize; i = i << 1) 
	{
        // std::cout << i << " "; std::cout.flush();
        resize(buf, i);
        randomize(buf);

        Pipe<Buffer<unsigned>, Source<> > src(buf);
        extString << src;

        TIter I = begin(extString);
        for(unsigned *cur = buf.begin; cur != buf.end; ++cur) {
            if (*cur != *I) {
                SEQAN_ASSERT_FAIL("testExternalString failed at position %u", cur - buf.begin);
                // not reached
            }
            ++I;
        }
    }
    freePage(buf, buf);
}



void testPool(unsigned maxSize = 16*1024*1024) {
    Buffer<unsigned> buf;
    allocPage(buf, maxSize, buf);

    Pool<unsigned,PoolSpec<> > pool;
    for(unsigned i = 1; i <= maxSize; i = i << 1) {
        /*
        std::cout << i;
        std::cout.flush();
        */

        resize(buf, i);
        randomize(buf);

        Pipe<Buffer<unsigned>, Source<> > src(buf);
        pool << src;

        /*
        if (pool.memBuffer.begin)
            std::cout << "* ";
        else 
            std::cout << " ";
        std::cout.flush();
        */

        beginRead(pool);
        for(unsigned *cur = buf.begin; cur != buf.end; cur++) {
            if (*cur != *pool) {
                endRead(pool);
                SEQAN_ASSERT_FAIL("testPool failed at position %u", cur - buf.begin);
                // not reached
            }
            ++pool;
        }
        endRead(pool);
    }
    freePage(buf, buf);
}



void testMapper(unsigned maxSize = 16*1024*1024) {
    Buffer<unsigned> buf;
    allocPage(buf, maxSize, buf);

    Pool<unsigned,MapperSpec<MapperConfig<IdentityMap<unsigned> > > > mapper;
    for(unsigned i = 1; i <= maxSize; i = i << 1) {
        /*
        std::cout << i;
        std::cout.flush();
        */

        resize(buf, i);
        permute(buf);

        Pipe<Buffer<unsigned>, Source<> > src(buf);
        mapper << src;

        /*(
        if (mapper.memBuffer.begin)
            std::cout << "* ";
        else 
            std::cout << " ";
        std::cout.flush();
        */

        beginRead(mapper);
        for(unsigned j = 0; j < i; ++j) {
            if (*mapper != j) {
                freePage(buf, buf);
                SEQAN_ASSERT_FAIL("testMapper failed at position %u", j);
                // not reached
            }	
            ++mapper;
        }
        endRead(mapper);
    }
    freePage(buf, buf);
}



void testPartiallyFilledMapper(unsigned maxSize = 16*1024*1024) {
    Buffer<unsigned> buf;
    allocPage(buf, maxSize, buf);

    Pool<unsigned,MapperSpec<MapperConfig<IdentityMap<unsigned> > > > mapper;
    for(unsigned i = 1; i <= maxSize; i = i << 1) {
        /*
        std::cout << i;
        std::cout.flush();
        */

        resize(buf, i);
        permute(buf);

        // partially fill the mapper 
        mapper.undefinedValue = i;	// select i as an undefined value (all defined values are less than i)
        resize(mapper, i);
        resize(buf, i - i/3);
        Pipe<Buffer<unsigned>, Source<> > src(buf);
        beginWrite(mapper) && append(mapper, src) && endWrite(mapper);

        /*
        if (mapper.memBuffer.begin)
            std::cout << "* ";
        else 
            std::cout << " ";
        std::cout.flush();
        */

        unsigned undefCounter = 0, missCounter = 0;
        beginRead(mapper);
        for(unsigned j = 0; j < i; ++j) {
            if (*mapper == i) 
                ++undefCounter;
            else
                if (*mapper != j) {
                    ++missCounter;
                    if (!mapper.memBuffer.begin) { // external mapping -> no misses allowed
                        freePage(buf, buf);
                        SEQAN_ASSERT_FAIL("testPartiallyFilledMapper failed at position %u [ = %u ]", j, *mapper);
                        // not reached
                    }
                }
            ++mapper;
        }
        endRead(mapper);
        if (mapper.memBuffer.begin) {
            if (undefCounter + missCounter > i/3) {
                SEQAN_ASSERT_FAIL("testPartiallyFilledMapper failed [only %u of %u undefind", undefCounter + missCounter, i / 3);
                // not reached
            }
        } else
            if (undefCounter != i/3) {
                SEQAN_ASSERT_FAIL("testPartiallyFilledMapper failed [only %u of %u undefined]", undefCounter + missCounter, i / 3);
                // not reached
            }
    }
    freePage(buf, buf);
}



void testSorter(unsigned maxSize = 16*1024*1024) {
    Buffer<unsigned> buf;
    allocPage(buf, maxSize, buf);

    Pool<unsigned,SorterSpec<SorterConfig<SimpleCompare<unsigned> > > > sorter;
    for(unsigned i = 1; i <= maxSize; i = i << 1) {
        /*
        std::cout << i;
        std::cout.flush();
        */

        resize(buf, i);
        permute(buf);

        Pipe<Buffer<unsigned>, Source<> > src(buf);
        sorter << src;

        /*
        if (sorter.memBuffer.begin)
            std::cout << "* ";
        else 
            std::cout << " ";
        std::cout.flush();
        */

        beginRead(sorter);
        unsigned j = *sorter, pos = 0;
        while (!eof(sorter)) {
            if (*sorter < j || /* *sorter < 0 || */ *sorter >= i) {
                freePage(buf, buf);
                SEQAN_ASSERT_FAIL("testSorter failed at position %u", pos);
                // not reached
            }
            j = *sorter;
            ++sorter; ++pos;
        }
        endRead(sorter);
    }
    freePage(buf, buf);
}


SEQAN_DEFINE_TEST(test_pipe_test_external_string) {
    testExternalString<MMap<> >(MAX_SIZE);
    testExternalString<External<> >(MAX_SIZE);
}


SEQAN_DEFINE_TEST(test_pipe_test_simple_pool) {
    testPool(MAX_SIZE);
}


SEQAN_DEFINE_TEST(test_pipe_test_mapper) {
    testMapper(MAX_SIZE);
}


SEQAN_DEFINE_TEST(test_pipe_test_mapper_partially_filled) {
    testPartiallyFilledMapper(MAX_SIZE);
}


SEQAN_DEFINE_TEST(test_pipe_test_sorter) {
    testSorter(MAX_SIZE);
}

template <typename TStringSet, typename T>
inline void appendValues(TStringSet &stringSet, T t)
{
    appendValue(stringSet, t);
}

template <typename TStringSet, typename T, typename... Args>
inline void appendValues(TStringSet &stringSet, T t, Args... args)
{
    // use recursive variadic function, instead of va_list
    appendValue(stringSet, t);
    appendValues(stringSet, args...) ;
}

template <typename TPipe, typename TStrings>
inline void comparePipeStream(TPipe &pipe, TStrings const &strings)
{
    typename Size<TStrings>::Type numStrings = length(strings);
    SEQAN_ASSERT_EQ(length(pipe), numStrings);
    beginRead(pipe);
    for (unsigned i = 0; i < numStrings; ++i)
    {
        SEQAN_ASSERT_NOT(eof(pipe));
        std::stringstream pipeRes;
        std::stringstream expRes;
        pipeRes << *pipe;
        expRes << strings[i];   // here we have to use the stringstreams as well for the same tokenization (of \0)
        SEQAN_ASSERT_EQ_MSG(CharString(pipeRes.str()), CharString(expRes.str()), "Pipe output differs at position %d", i);
        ++pipe;
    }
    SEQAN_ASSERT(eof(pipe));
    endRead(pipe);
}

template <typename TPack>
void testPipeSampler()
{
    {
        typedef DnaString TString;
        typedef Pipe<TString, Source<> > TSource;
        typedef Pipe<TSource, Sampler<3, TPack> > TSamplerDC3;

        TString string = "";
        TSource source(string);
        TSamplerDC3 sampler(source);


        string = "";
        SEQAN_ASSERT_EQ(length(source), 0u);
        StringSet<CharString> expectedOutput;
        comparePipeStream(sampler, expectedOutput);


        string = "T";
        SEQAN_ASSERT_EQ(length(source), 1u);
        appendValue(expectedOutput, "< 1 , [T A A] >");
        comparePipeStream(sampler, expectedOutput);


        string = "TA";
        SEQAN_ASSERT_EQ(length(source), 2u);
        clear(expectedOutput);
        appendValues(expectedOutput,
            "< 2 , [T A A] >",
            "< 1 , [A A A] >");
        comparePipeStream(sampler, expectedOutput);


        string = "TAC";
        SEQAN_ASSERT_EQ(length(source), 3u);
        clear(expectedOutput);
        appendValues(expectedOutput,
            "< 2 , [A C A] >",
            "< 1 , [C A A] >");
        comparePipeStream(sampler, expectedOutput);


        string = "TACG";
        SEQAN_ASSERT_EQ(length(source), 4u);
        clear(expectedOutput);
        appendValues(expectedOutput,
            "< 4 , [T A C] >",
            "< 2 , [C G A] >",
            "< 1 , [G A A] >");
        comparePipeStream(sampler, expectedOutput);
    }
    {
        typedef StringSet<DnaString> TStringSet;
        typedef typename Concatenator<TStringSet>::Type TConcat;
        typedef Pipe<TConcat, Source<> > TSource;
        typedef Pipe<
            TSource,
            Multi<
                Sampler<7, BitPacked<> >,
                Pair<unsigned short, unsigned short, Pack>,
                StringSetLimits<TStringSet>::Type> >    TSamplerDC7;

        TStringSet set;
        appendValue(set, "CGGAAGGCC");
        appendValue(set, "A");
        appendValue(set, "");
        appendValue(set, "ACCGT");
        appendValue(set, "CGGAAGGCCA");
        appendValue(set, "TCGGAAGGCCA");
        appendValue(set, "TTCGGAAGGCCA");
        appendValue(set, "TTTCGGAAGGCCA");
        appendValue(set, "TTTTCGGAAGGCCA");
        appendValue(set, "TTTTTCGGAAGGCCA");
        appendValue(set, "TTTTTTCGGAAGGCCA");

        TSource source(concat(set));
        TSamplerDC7 sampler(source, stringSetLimits(set));

        SEQAN_ASSERT_EQ(&stringSetLimits(set), &sampler.limits);
        SEQAN_ASSERT_EQ(length(source), 106u);

        StringSet<CharString> expectedOutput;
        appendValues(expectedOutput,
            "< < 0 , 9 > , [C G G A A G G] >",
            "< < 0 , 8 > , [G G A A G G C] >",
            "< < 0 , 4 > , [G G C C A A C] >",
            "< < 0 , 2 > , [C C A A C C G] >",
            "< < 0 , 1 > , [C A A C C G T] >",
            "< < 1 , 1 > , [A A C C G T C] >",
            "< < 3 , 4 > , [C C G T C G G] >",
            "< < 3 , 2 > , [G T C G G A A] >",
            "< < 3 , 1 > , [T C G G A A G] >",
            "< < 4 , 9 > , [G G A A G G C] >",
            "< < 4 , 8 > , [G A A G G C C] >",
            "< < 4 , 4 > , [G C C A T C G] >",
            "< < 4 , 2 > , [C A T C G G A] >",
            "< < 4 , 1 > , [A T C G G A A] >",
            "< < 5 , 11 > , [T C G G A A G] >",
            "< < 5 , 9 > , [G G A A G G C] >",
            "< < 5 , 8 > , [G A A G G C C] >",
            "< < 5 , 4 > , [G C C A T T C] >",
            "< < 5 , 2 > , [C A T T C G G] >",
            "< < 5 , 1 > , [A T T C G G A] >",
            "< < 6 , 11 > , [T C G G A A G] >",
            "< < 6 , 9 > , [G G A A G G C] >",
            "< < 6 , 8 > , [G A A G G C C] >",
            "< < 6 , 4 > , [G C C A T T T] >",
            "< < 6 , 2 > , [C A T T T C G] >",
            "< < 6 , 1 > , [A T T T C G G] >",
            "< < 7 , 11 > , [T C G G A A G] >",
            "< < 7 , 9 > , [G G A A G G C] >",
            "< < 7 , 8 > , [G A A G G C C] >",
            "< < 7 , 4 > , [G C C A T T T] >",
            "< < 7 , 2 > , [C A T T T T C] >",
            "< < 7 , 1 > , [A T T T T C G] >",
            "< < 8 , 11 > , [T C G G A A G] >",
            "< < 8 , 9 > , [G G A A G G C] >",
            "< < 8 , 8 > , [G A A G G C C] >",
            "< < 8 , 4 > , [G C C A T T T] >",
            "< < 8 , 2 > , [C A T T T T T] >",
            "< < 8 , 1 > , [A T T T T T C] >",
            "< < 9 , 15 > , [T T T T T C G] >",
            "< < 9 , 11 > , [T C G G A A G] >",
            "< < 9 , 9 > , [G G A A G G C] >",
            "< < 9 , 8 > , [G A A G G C C] >",
            "< < 9 , 4 > , [G C C A T T T] >",
            "< < 9 , 2 > , [C A T T T T T] >",
            "< < 9 , 1 > , [A T T T T T T] >",
            "< < 10 , 16 > , [T T T T T T C] >",
            "< < 10 , 15 > , [T T T T T C G] >",
            "< < 10 , 11 > , [T C G G A A G] >",
            "< < 10 , 9 > , [G G A A G G C] >",
            "< < 10 , 8 > , [G A A G G C C] >",
            "< < 10 , 4 > , [G C C A A A A] >",
            "< < 10 , 2 > , [C A A A A A A] >",
            "< < 10 , 1 > , [A A A A A A A] >");
        comparePipeStream(sampler, expectedOutput);
    }
}

SEQAN_DEFINE_TEST(test_pipe_sampler)
{
    testPipeSampler<Pack>();
    testPipeSampler<BitPacked<> >();
}

template <bool omitLast>
void testPipeTupler()
{
    String<Peptide> testTexts;

    Peptide str;
    appendValue(testTexts, str);
    for (int i = 0; i <= 9; ++i)
    {
        appendValue(str, AminoAcid(i));
        appendValue(testTexts, str);
    }

    typedef Tupler<4, omitLast, Pack>           TTuplerSpec1;
    typedef Tupler<4, omitLast, BitPacked<> >   TTuplerSpec2;

    typedef Pipe<Peptide, Source<> >            TSource;
    typedef Pipe<TSource, TTuplerSpec1>         TTupler1;
    typedef Pipe<TSource, TTuplerSpec2>         TTupler2;

    TSource src(back(testTexts));
    TTupler1 tupler(src);

    StringSet<CharString> expectedOutput;
    if (omitLast)
    {
        appendValues(expectedOutput,
            "< 0 , [A B C D] >",
            "< 1 , [B C D E] >",
            "< 2 , [C D E F] >",
            "< 3 , [D E F G] >",
            "< 4 , [E F G H] >",
            "< 5 , [F G H I] >",
            "< 6 , [G H I J] >");
    }
    else
    {
        appendValues(expectedOutput,
            "< 0 , [A B C D] >",
            "< 1 , [B C D E] >",
            "< 2 , [C D E F] >",
            "< 3 , [D E F G] >",
            "< 4 , [E F G H] >",
            "< 5 , [F G H I] >",
            "< 6 , [G H I J] >",
            "< 7 , [H I J A] >",
            "< 8 , [I J A A] >",
            "< 9 , [J A A A] >");
    }
    comparePipeStream(tupler, expectedOutput);

    for (unsigned i = 0; i < length(testTexts); ++i)
    {
        TSource src1(testTexts[i]);
        TTupler1 tupler1(src1);

        TSource src2(testTexts[i]);
        TTupler2 tupler2(src2);

        SEQAN_ASSERT_EQ((int)length(tupler1), std::max(0, (int)length(testTexts[i]) - (omitLast? 3 : 0)));
        comparePipes(tupler1, tupler2);
    }
}

template <bool omitLast>
void testPipeMultiTupler()
{
    typedef StringSet<Peptide> TStringSet;
    String<TStringSet> testTexts;

    TStringSet set;
    appendValue(testTexts, set);
    appendValue(set, "");
    for (int i = 0; i <= 9; ++i)
    {
        appendValue(back(set), AminoAcid(i));
        appendValue(testTexts, set);
    }
    appendValue(set, "ABC");
    appendValue(testTexts, set);
    appendValue(set, "D");
    appendValue(testTexts, set);
    appendValue(set, "");
    appendValue(testTexts, set);
    appendValue(set, "IHGEQCDNRA");
    appendValue(testTexts, set);

    typedef typename Concatenator<TStringSet>::Type         TConcat;
    typedef typename StringSetLimits<TStringSet>::Type      TLimits;

    typedef Multi<
        Tupler<4, omitLast, Pack>,
        Pair<unsigned short, unsigned short, Pack>,
        TLimits>                                            TTuplerSpec1;

    typedef Multi<
        Tupler<4, omitLast, BitPacked<> >,
        Pair<unsigned, unsigned, Pack>,
        TLimits>                                            TTuplerSpec2;

    typedef Pipe<TConcat, Source<> >                        TSource;
    typedef Pipe<TSource, TTuplerSpec1>                     TTupler1;
    typedef Pipe<TSource, TTuplerSpec2>                     TTupler2;

    TSource src(concat(back(testTexts)));
    TTupler1 tupler(src, stringSetLimits(back(testTexts)));

    StringSet<CharString> expectedOutput;
    if (omitLast)
    {
        appendValues(expectedOutput,
            "< < 0 , 0 > , [A B C D] >",
            "< < 0 , 1 > , [B C D E] >",
            "< < 0 , 2 > , [C D E F] >",
            "< < 0 , 3 > , [D E F G] >",
            "< < 0 , 4 > , [E F G H] >",
            "< < 0 , 5 > , [F G H I] >",
            "< < 0 , 6 > , [G H I J] >",
            "< < 4 , 0 > , [I H G E] >",
            "< < 4 , 1 > , [H G E Q] >",
            "< < 4 , 2 > , [G E Q C] >",
            "< < 4 , 3 > , [E Q C D] >",
            "< < 4 , 4 > , [Q C D N] >",
            "< < 4 , 5 > , [C D N R] >",
            "< < 4 , 6 > , [D N R A] >");
    }
    else
    {
        appendValues(expectedOutput,
            "< < 0 , 0 > , [A B C D] >",
            "< < 0 , 1 > , [B C D E] >",
            "< < 0 , 2 > , [C D E F] >",
            "< < 0 , 3 > , [D E F G] >",
            "< < 0 , 4 > , [E F G H] >",
            "< < 0 , 5 > , [F G H I] >",
            "< < 0 , 6 > , [G H I J] >",
            "< < 0 , 7 > , [H I J A] >",
            "< < 0 , 8 > , [I J A A] >",
            "< < 0 , 9 > , [J A A A] >",
            "< < 1 , 0 > , [A B C A] >",
            "< < 1 , 1 > , [B C A A] >",
            "< < 1 , 2 > , [C A A A] >",
            "< < 2 , 0 > , [D A A A] >",
            "< < 4 , 0 > , [I H G E] >",
            "< < 4 , 1 > , [H G E Q] >",
            "< < 4 , 2 > , [G E Q C] >",
            "< < 4 , 3 > , [E Q C D] >",
            "< < 4 , 4 > , [Q C D N] >",
            "< < 4 , 5 > , [C D N R] >",
            "< < 4 , 6 > , [D N R A] >",
            "< < 4 , 7 > , [N R A A] >",
            "< < 4 , 8 > , [R A A A] >",
            "< < 4 , 9 > , [A A A A] >");
    }
    comparePipeStream(tupler, expectedOutput);

    for (unsigned i = 0; i < length(testTexts); ++i)
    {
        TSource src1(concat(testTexts[i]));
        TTupler1 tupler1(src1, stringSetLimits(testTexts[i]));

        TSource src2(concat(testTexts[i]));
        TTupler2 tupler2(src2, stringSetLimits(testTexts[i]));

        comparePipes(tupler1, tupler2);
    }
}

SEQAN_DEFINE_TEST(test_pipe_tupler)
{
    testPipeTupler<false>();
    testPipeTupler<true>();
}

SEQAN_DEFINE_TEST(test_pipe_tupler_multi)
{
    testPipeMultiTupler<false>();
    testPipeMultiTupler<true>();
}

SEQAN_BEGIN_TESTSUITE(test_pipe) {
	std::cerr << "";  // This line is an esoteric fix for an even more esoteric crash in MS VC++ 9/10.
    SEQAN_CALL_TEST(test_pipe_test_external_string);
    SEQAN_CALL_TEST(test_pipe_test_simple_pool);
    SEQAN_CALL_TEST(test_pipe_test_mapper);
    SEQAN_CALL_TEST(test_pipe_test_mapper_partially_filled);
    SEQAN_CALL_TEST(test_pipe_test_sorter);
    SEQAN_CALL_TEST(test_pipe_sampler);
    SEQAN_CALL_TEST(test_pipe_tupler);
    SEQAN_CALL_TEST(test_pipe_tupler_multi);
}
SEQAN_END_TESTSUITE

