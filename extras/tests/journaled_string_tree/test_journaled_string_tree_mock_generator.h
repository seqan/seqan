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
// Implements a mock generator to generate the different test cases.
// ==========================================================================

#ifndef EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_DATA_PARALLEL_MOCK_GENERATOR_H_
#define EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_DATA_PARALLEL_MOCK_GENERATOR_H_

#include <iostream>

#include <seqan/basic.h>
#include <seqan/random.h>
#include <seqan/seq_io.h>
#include <seqan/journaled_string_tree.h>

namespace seqan
{

// ----------------------------------------------------------------------------
// Class MockVariantData
// ----------------------------------------------------------------------------

template <typename TValue>
struct MockVariantData
{
    unsigned _hostPos;
    unsigned _delSize;
    DeltaType::TValue _varType;
    String<TValue> _insBuff;

    MockVariantData() : _hostPos(0), _delSize(0), _varType(DeltaType::DELTA_TYPE_INDEL)
    {}
};

struct TestRandomGen_
{
    Rng<MersenneTwister> rng;

    TestRandomGen_()
    {
        rng = Rng<MersenneTwister>(42);
    }
};

static TestRandomGen_ globalRng;

// ----------------------------------------------------------------------------
// Class DataParallelTestConfig
// ----------------------------------------------------------------------------

template <typename TValue>
struct DataParallelTestConfig
{

    StringSet<String<unsigned> > _posConfig;
    StringSet<String<MockVariantData<TValue> > > _varConfig;
    StringSet<String<String<bool, Packed<> > > > _covConfig;

    DataParallelTestConfig()
    {
        _readConfigFile(*this);
    }

    void getTestConfiguration(String<MockVariantData<TValue> > & varData,
                              String<String<bool, Packed<> > > & covData,
                              unsigned posConfig,
                              unsigned varConfig,
                              unsigned covConfig)
    {

        varData = _varConfig[varConfig];
        for (unsigned i = 0; i < length(varData); ++i)
            varData[i]._hostPos = _posConfig[posConfig][i];
        covData = _covConfig[covConfig];
    }
};

template <typename TSize, typename TAlphabet>
struct MockGenerator_
{
    typedef DeltaMap<TSize, TAlphabet> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TStringTree;
    typedef typename GetStringSet<TStringTree>::Type TJournalSet;
    typedef typename Value<TJournalSet>::Type TJournalString;
    typedef typename Host<TJournalSet>::Type THost;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_DEL>::Type TDel;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INDEL>::Type TIndel;

    TJournalSet _seqData;
    TDeltaMap _varStore;

    MockGenerator_() : _seqData(), _varStore()
    {}

    template <typename TVar, typename TCov, typename TSize2>
    void generate(String<TVar> & varData, String<TCov> & covData, TSize2 const & refSize)
    {
        typedef typename Host<TJournalString>::Type THost;
        SEQAN_ASSERT(!empty(varData));
        SEQAN_ASSERT(!empty(covData));

        unsigned numSeq  = length(covData[0]);
        setCoverageSize(_varStore._deltaCoverageStore, numSeq);

        THost ref;
        generateRef(ref, refSize);

        clear(_seqData);
        // Generate the Journal Data.
        createHost(this->_seqData, ref);

        for (unsigned i = 0; i < numSeq; ++i)
        {
            TJournalString tmp;
            setHost(tmp, host(this->_seqData));
            int virtPos = 0;
            unsigned lastPhysPos = 0;
            // We need to parse the coverage. And insert the variant for all sequences at this coverage.
            for (unsigned j = 0; j < length(covData); ++j)
            {
                if (covData[j][i])  // Seq i covers the jth variant.
                {
                    unsigned offset = varData[j]._hostPos - lastPhysPos;
                    virtPos += offset;
                    lastPhysPos = varData[j]._hostPos;

                    switch(varData[j]._varType)
                    {
                    case DeltaType::DELTA_TYPE_SNP:
                        assignValue(tmp, virtPos, varData[j]._insBuff[0]);
                        break;
                    case DeltaType::DELTA_TYPE_DEL:
                        erase(tmp, virtPos, virtPos + varData[j]._delSize);
                        virtPos -= static_cast<int>(varData[j]._delSize);
                        break;
                    case DeltaType::DELTA_TYPE_INS:
                        insert(tmp, virtPos, varData[j]._insBuff);
                        virtPos += length(varData[j]._insBuff);
                        break;
                    case DeltaType::DELTA_TYPE_INDEL:
                        replace(tmp, virtPos, virtPos + varData[j]._delSize,  varData[j]._insBuff);
                        virtPos += static_cast<int>(length(varData[j]._insBuff)) - static_cast<int>(varData[j]._delSize);
                        break;
                    }
                }
            }
            appendValue(_seqData, tmp);
        }

        //Generate the DeltaMap.
        for (unsigned i = 0; i < length(varData); ++i)
        {
            switch(varData[i]._varType)
            {
            case DeltaType::DELTA_TYPE_SNP:
                insert(_varStore, varData[i]._hostPos, varData[i]._insBuff[0], covData[i]);
                break;
            case DeltaType::DELTA_TYPE_DEL:
                insert(_varStore, varData[i]._hostPos, static_cast<TDel>(varData[i]._delSize), covData[i]);
                break;
            case DeltaType::DELTA_TYPE_INS:
                insert(_varStore, varData[i]._hostPos, varData[i]._insBuff, covData[i]);
                break;
            case DeltaType::DELTA_TYPE_INDEL:
                insert(_varStore, varData[i]._hostPos, TIndel(varData[i]._delSize, varData[i]._insBuff), covData[i]);
                break;
            }
        }
    }
};

template <typename THost, typename TSize>
void generateRef(THost & seq,
                 TSize const & newLength)
{
    typedef typename Value<THost>::Type THostValue;
    typedef typename Iterator<THost>::Type TIter;

    Pdf<Uniform<unsigned> > pdf(65, 90);

    resize(seq, newLength, Exact());
    TIter it = begin(seq);
    TIter itEnd = end(seq);

    while(it != itEnd)
    {
        *it = static_cast<THostValue>(pickRandomNumber(globalRng.rng, pdf));
        ++it;
    }
}

// ----------------------------------------------------------------------------
// Function _readConfigFile
// ----------------------------------------------------------------------------

template <typename TValue>
int _readConfigFile(DataParallelTestConfig<TValue> & obj)
{
    CharString confFile = SEQAN_PATH_TO_ROOT();
    append(confFile, "/extras/tests/journaled_string_tree/test_journaled_string_tree_find.conf");

    std::ifstream inputFile;
    inputFile.open(toCString(confFile), std::ios_base::in);

    if (!inputFile.good())
    {
        std::cerr << "Cannot open config file <" << confFile << ">" << std::endl;
        return -1;
    }

    RecordReader<std::ifstream, SinglePass<> > reader(inputFile);

    // Read positions.
    if (value(reader) != '>')
        return -2;
    skipLine(reader);  // Skip the first id line.

    // Read the positions in the config file.
    while(value(reader) != '>')
    {
        unsigned currConf = length(obj._posConfig);
        resize(obj._posConfig, currConf + 1);
        while(value(reader) != '\n')
        {
            CharString buffer;
            readUntilChar(buffer, reader, ',');
            skipNChars(reader, 1);
            unsigned pos;
            lexicalCast2(pos, buffer);
            appendValue(obj._posConfig[currConf], pos);
        }
        skipLine(reader);
    }

    // Read the variant information
    SEQAN_ASSERT_EQ(value(reader), '>');
    skipLine(reader);
    while(value(reader) != '>')
    {
        unsigned currConf = length(obj._varConfig);
        resize(obj._varConfig, currConf + 1);
        while(value(reader) != '\n')
        {
            CharString buffer;
            readNChars(buffer, reader, 1);
            MockVariantData<TValue> varData;
            switch(buffer[0])
            {
            case 's':  // Read SNP data.
            {
                varData._varType = DeltaType::DELTA_TYPE_SNP;
                CharString snpBuff;
                readUntilChar(snpBuff, reader, ',');
                unsigned snpSize;
                lexicalCast2(snpSize, snpBuff);
                generateRef(varData._insBuff, snpSize);
                break;
            }
            case 'd':
            {
                varData._varType = DeltaType::DELTA_TYPE_DEL;
                CharString delBuff;
                readUntilChar(delBuff, reader, ',');
                lexicalCast2(varData._delSize, delBuff);
                break;
            }
            case 'i':
            {
                varData._varType = DeltaType::DELTA_TYPE_INS;
                CharString insBuff;
                readUntilChar(insBuff, reader, ',');
                unsigned insSize;
                lexicalCast2(insSize, insBuff);
                generateRef(varData._insBuff, insSize);
                break;
            }
            case 'r':
            {
                varData._varType = DeltaType::DELTA_TYPE_INDEL;
                CharString buff;
                readUntilChar(buff, reader, '-');
                lexicalCast2(varData._delSize, buff);
                skipNChars(reader, 1);
                clear(buff);
                readUntilChar(buff, reader, ',');
                unsigned insSize;
                lexicalCast2(insSize, buff);
                generateRef(varData._insBuff, insSize);
                break;
            }
            }
            skipNChars(reader, 1);
            appendValue(obj._varConfig[currConf], varData);
        }
        skipLine(reader);
    }

    // Read the coverage information
    SEQAN_ASSERT_EQ(value(reader), '>');
    skipLine(reader);

    while(!atEnd(reader) && value(reader) != '>')
    {
        unsigned currConf = length(obj._covConfig);
        resize(obj._covConfig, currConf + 1);
        while(value(reader) != '\n')
        {
            CharString buffer;
            readUntilChar(buffer, reader, ',');
            skipNChars(reader, 1);
            String<bool, Packed<> > tmpVal;
            resize(tmpVal, length(buffer), false, Exact());
            for (unsigned i = 0; i < length(tmpVal); ++i)
                if (buffer[i] == '1')
                    tmpVal[i] = true;

            appendValue(obj._covConfig[currConf], tmpVal);
        }
        skipLine(reader);
    }
    return 0;
    }

}  // namespace seqan

#endif  // EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_DATA_PARALLEL_MOCK_GENERATOR_H_
