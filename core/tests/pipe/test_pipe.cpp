// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#include <seqan/file.h>
#include <seqan/pipe.h>
#include <seqan/system.h>

#include "test_pipe.h"

using namespace std;
using namespace seqan;


// Maximum size of data to use for the external tests.
const size_t MAX_SIZE = 1u << 20;


// TODO(holtgrew): The following test* functions should actually be defined with SEQAN_DEFINE_TEST().


template <typename TStringSpec>
void testExternalString(unsigned maxSize = 16*1024*1024) 
{
	typedef String<unsigned, TStringSpec> TExtString;
    typedef typename Iterator<TExtString const, Standard>::Type TIter;

    Buffer<unsigned> buf;
    allocPage(buf, maxSize, buf);

    TExtString extString;
    for(unsigned i = 1; i <= maxSize; i = i << 1) 
	{
        // ::std::cout << i << " "; ::std::cout.flush();
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
        ::std::cout << i;
        ::std::cout.flush();
        */

        resize(buf, i);
        randomize(buf);

        Pipe<Buffer<unsigned>, Source<> > src(buf);
        pool << src;

        /*
        if (pool.memBuffer.begin)
            ::std::cout << "* ";
        else 
            ::std::cout << " ";
        ::std::cout.flush();
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
        ::std::cout << i;
        ::std::cout.flush();
        */

        resize(buf, i);
        permute(buf);

        Pipe<Buffer<unsigned>, Source<> > src(buf);
        mapper << src;

        /*(
        if (mapper.memBuffer.begin)
            ::std::cout << "* ";
        else 
            ::std::cout << " ";
        ::std::cout.flush();
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
        ::std::cout << i;
        ::std::cout.flush();
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
            ::std::cout << "* ";
        else 
            ::std::cout << " ";
        ::std::cout.flush();
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
        ::std::cout << i;
        ::std::cout.flush();
        */

        resize(buf, i);
        permute(buf);

        Pipe<Buffer<unsigned>, Source<> > src(buf);
        sorter << src;

        /*
        if (sorter.memBuffer.begin)
            ::std::cout << "* ";
        else 
            ::std::cout << " ";
        ::std::cout.flush();
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


SEQAN_BEGIN_TESTSUITE(test_pipe) {
	std::cerr << "";  // This line is an esoteric fix for an even more esoteric crash in MS VC++ 9/10.
    SEQAN_CALL_TEST(test_pipe_test_external_string);
    SEQAN_CALL_TEST(test_pipe_test_simple_pool);
    SEQAN_CALL_TEST(test_pipe_test_mapper);
    SEQAN_CALL_TEST(test_pipe_test_mapper_partially_filled);
    SEQAN_CALL_TEST(test_pipe_test_sorter);
}
SEQAN_END_TESTSUITE

