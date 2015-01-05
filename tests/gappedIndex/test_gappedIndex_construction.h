// ==========================================================================
//                      test_gappedIndex_construction.h
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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#ifndef CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_CONSTRUCTION_H_
#define CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_CONSTRUCTION_H_

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
// _createAndCompareSAs: only for fixed shapes, because of External Dislex
// ----------------------------------------------------------------------------

template <typename TStr, typename TShape>
bool _createAndCompareSAs(TStr const & text,
                     TShape const & shape)
{
    String<unsigned> sa1, sa2;
    resize(sa1, length(text));
    createGappedSuffixArray(sa1, text, shape, ModCyclicShape<TShape>(), SAQSort());
    
    clear(sa2);
    resize(sa2, length(text));
    createGappedSuffixArray(sa2, text, shape, ModCyclicShape<TShape>(), DislexExternal<TShape>());
    if (sa1 != sa2) return false;
    
//    clear(sa2);
//    resize(sa2, length(text));
//    createGappedSuffixArray(sa2, text, shape, ModCyclicShape<TShape>(), InplaceRadixSort());
//    if (sa1 != sa2) return false;
//    
//    clear(sa2);
//    resize(sa2, length(text));
//    createGappedSuffixArray(sa2, text, shape, ModCyclicShape<TShape>(), Dislex<Skew7>());
//    if (sa1 != sa2) return false;
    
    return true;

}

template <typename TStr, typename TSetSpec, typename TShape>
bool _createAndCompareSAs(StringSet<TStr, TSetSpec> const & text,
                          TShape const & shape)
{
    String<Pair<unsigned, unsigned> > sa1, sa2;
    resize(sa1, lengthSum(text));
    createGappedSuffixArray(sa1, text, shape, ModCyclicShape<TShape>(), SAQSort());
    
    clear(sa2);
    resize(sa2, lengthSum(text));
    createGappedSuffixArray(sa2, text, shape, ModCyclicShape<TShape>(), DislexExternal<TShape>());
    if (sa1 != sa2) return false;
    
//    clear(sa2);
//    resize(sa2, lengthSum(text));
//    createGappedSuffixArray(sa2, text, shape, ModCyclicShape<TShape>(), InplaceRadixSort());
//    if (sa1 != sa2) return false;
//    
//    clear(sa2);
//    resize(sa2, lengthSum(text));
//    createGappedSuffixArray(sa2, text, shape, ModCyclicShape<TShape>(), Dislex<>());
//    if (sa1 != sa2) return false;

    return true;
}



SEQAN_DEFINE_TEST(test_gappedIndex_construction_str)
{
    Rng<> rng(/*seed=*/100);
    
    typedef String<char> TText;
    typedef String<unsigned> TArray;
    
    TText   text;
    TArray  sa1;
    TArray  sa2;
    
    const int maxSize = 1024; //256 * 1024;	// max text size is 1 megabyte
    int runs = 5;
    
    for(int i = 0; i < runs; ++i) {
        
        std::cout << "*** RUN " << i << " ***";
        
        Pdf<Uniform<int> > pdf(0, maxSize);
        int size = pickRandomNumber(rng, pdf);
        
        //___randomize_text___________________________________________________________
        resize(text,size);
        textRandomize(text);
        std::cout << "   textSize: " << length(text) << std::endl;
        
        //___create_suffix_array______________________________________________________
        TestGappedIndexShapeDefs_ SD;
        
        if ( !_createAndCompareSAs(text, SD.S_10) )      std::cout << "SAs differ!!" << std::endl;
        if ( !_createAndCompareSAs(text, SD.S_11010) )   std::cout << "SAs differ!!" << std::endl;
        if ( !_createAndCompareSAs(text, SD.S_111100) )  std::cout << "SAs differ!!" << std::endl;
        if ( !_createAndCompareSAs(text, SD.S_10001) )   std::cout << "SAs differ!!" << std::endl;
        if ( !_createAndCompareSAs(text, SD.S_01) )      std::cout << "SAs differ!!" << std::endl;
        if ( !_createAndCompareSAs(text, SD.S_0011) )    std::cout << "SAs differ!!" << std::endl;
    }
        
}


SEQAN_DEFINE_TEST(test_gappedIndex_construction_strSet)
{
    Rng<> rng(/*seed=*/100);
    
    typedef StringSet<CharString> TText;
    typedef String<Pair<unsigned, unsigned> > TArray;
    
    TText   text;
    TArray  sa1;
    TArray  sa2;
    
    const int maxSize = 1024; // 16 * 1024;	// max text size is 30 kilobyte per String
    const int maxSet = 64; // 1024;
    int runs = 5;
    int totalSize =0;
    
    for(int i = 0; i < runs; ++i) {
        
        std::cout << "*** RUN " << i << " ***";
        
        Pdf<Uniform<int> > pdf(10, maxSet);
        int set = pickRandomNumber(rng, pdf);
        
        for (int j=0; j< set; ++j)
        {
            Pdf<Uniform<int> > pdf2(0, maxSize);
            int size = pickRandomNumber(rng, pdf);
            totalSize += size;
            CharString x;
            resize(x,size);
            textRandomize(x);
            appendValue(text, x);
        }
        
        std::cout << "   totalSize: " << totalSize << "\tsetSize: " << length(text) << std::endl;
        
        
        //___create_suffix_array______________________________________________________
        TestGappedIndexShapeDefs_ SD;
        
        if ( !_createAndCompareSAs(text, SD.S_10) )      std::cout << "SAs differ!!" << std::endl;
        if ( !_createAndCompareSAs(text, SD.S_11010) )   std::cout << "SAs differ!!" << std::endl;
        if ( !_createAndCompareSAs(text, SD.S_111100) )  std::cout << "SAs differ!!" << std::endl;
        if ( !_createAndCompareSAs(text, SD.S_10001) )   std::cout << "SAs differ!!" << std::endl;
        if ( !_createAndCompareSAs(text, SD.S_01) )      std::cout << "SAs differ!!" << std::endl;
        if ( !_createAndCompareSAs(text, SD.S_0011) )    std::cout << "SAs differ!!" << std::endl;
    }
    
}



#endif  // #ifndef CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_CONSTRUCTION_H_
