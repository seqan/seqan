// ==========================================================================
//                                  FMIndex
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_BITSTRING_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_BITSTRING_H_

namespace seqan {

//forward declaration
template <typename TSpec>
struct RankSupportBitString;

struct SmallString_;
struct MediumString_;
struct LargeString_;

struct FibreRankSupportBitString_;
struct FibreRankSupportBucketString_;
struct FibreRankSupportSuperBucketString_;

typedef Tag<SmallString_> const SmallString;
typedef Tag<MediumString_> const MediumString;
typedef Tag<LargeString_> const LargeString;

typedef Tag<FibreRankSupportBitString_> const FibreRankSupportBitString;
typedef Tag<FibreRankSupportBucketString_> const FibreRankSupportBucketString;
typedef Tag<FibreRankSupportSuperBucketString_> const FibreRankSupportSuperBucketString;

template <typename TSpec>
struct Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>
{
    typedef String<unsigned long> Type;
};

template <typename TSpec>
struct Fibre<RankSupportBitString<TSpec>, FibreRankSupportBucketString>
{
    typedef String<unsigned short> Type;
};

template <typename TSpec>
struct Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>
{
    typedef String<unsigned int> Type;
};

template <typename TSpec>
struct Size<RankSupportBitString<TSpec> >
{
    typedef typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type>::Type Type;
};

template <typename TSpec>
struct RankSupportBitString
{
    typedef typename Fibre<RankSupportBitString, FibreRankSupportBitString>::Type             TBitString;
    typedef typename Fibre<RankSupportBitString, FibreRankSupportBucketString>::Type          TBucketString;
    typedef typename Fibre<RankSupportBitString, FibreRankSupportSuperBucketString>::Type     TSuperBucketString;

    TBitString                                      bitString;
    TBucketString                                   bucketString;
    TSuperBucketString                              superBucketString;
    typename Value<TSuperBucketString>::Type   	    length;

    RankSupportBitString() :
    	bitString(),
    	bucketString(),
    	superBucketString(),
    	length(0)
    {}

    template <typename TText>
    RankSupportBitString(TText & text) :
		bitString(),
		bucketString(),
		superBucketString(),
        length(length(text))
    {
        unsigned const bitsPerBucket = BitsPerValue<typename Value<TBitString>::Type>::VALUE;

        TSuperBucketString textLength = length(text);
        (length(text) % bitsPerBucket == 0) ? resize(bitString, textLength / bitsPerBucket) : resize(bitString, textLength / bitsPerBucket + 1);

        unsigned const bitsPerBucketStringEntrie = BitsPerValue<TBucketString>::VALUE;
        resize(bucketString, length(bitString) / bitsPerBucketStringEntrie);

        unsigned const bitsPerSuperBucketStringEntrie = bitsPerBucketStringEntrie * bitsPerBucketStringEntrie;
        resize(superBucketString, length(bucketString) / bitsPerSuperBucketStringEntrie);
    }

    inline RankSupportBitString & operator=(RankSupportBitString const & other)
    {
    	bitString = other.bitString;
    	bucketString = other.bucketString;
    	superBucketString = other.superBucketString;
    	length = other.length;
    	return *this;
    }

    inline bool operator==(const RankSupportBitString & b) const
    {
        return length == b.length &&
               bitString == b.bitString &&
               bucketString == b.bucketString &&
               superBucketString == b.superBucketString;
    }

    template <typename TPos>
    inline typename Value<TSuperBucketString>::Type const operator[](TPos pos)
    {
        unsigned const bitsPerBucket = BitsPerValue<typename Value<TBitString>::Type>::VALUE;
        return bitString[pos / bitsPerBucket];
    }

};

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type &
getFibre(RankSupportBitString<TSpec> & string, const FibreRankSupportBitString)
{
    return string.bitString;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type const &
getFibre(RankSupportBitString<TSpec> const & string, const FibreRankSupportBitString)
{
    return string.bitString;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBucketString>::Type &
getFibre(RankSupportBitString<TSpec> & string, const FibreRankSupportBucketString)
{
    return string.bucketString;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBucketString>::Type const &
getFibre(RankSupportBitString<TSpec> const & string, const FibreRankSupportBucketString)
{
    return string.bucketString;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type &
getFibre(RankSupportBitString<TSpec> & string, const FibreRankSupportSuperBucketString)
{
    return string.superBucketString;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type const &
getFibre(const RankSupportBitString<TSpec> & string, const FibreRankSupportSuperBucketString)
{
    return string.superBucketString;
}

template <typename TValue>
inline void printBits(TValue entrie)
{
    unsigned bitsPerValue = BitsPerValue<TValue>::VALUE;
    TValue one = 1;
    std::cerr << "entrie: " << entrie << std::endl;
    std::cerr << bitsPerValue << std::endl;
    for (TValue i = 0; i < bitsPerValue; ++i)
    {
        std::cerr << ((entrie >> i) & one);
    }
    std::cerr << std::endl;
}

template <typename TValue, typename TSize>
inline std::ostream & printBits(std::ostream & stream, TValue entrie, TSize blockSize)
{
    unsigned bitsPerValue = BitsPerValue<TValue>::VALUE;
    bool temp;
    for (int i = bitsPerValue - 1; i >= 0; --i)
    {
        temp = (entrie >> i) & 1;
        stream << temp;
        if ((bitsPerValue - i) % blockSize == 0)
            stream << " ";
    }
    return stream;
}

template <typename TSpec>
inline std::ostream & operator<<(std::ostream & stream, const RankSupportBitString<TSpec> & rankSupportBitString)
{
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type                TBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBucketString>::Type             TBucketString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type        TSuperBucketString;

    typedef typename Value<TBitString>::Type          TBitStringValue;
    typedef typename Value<TBucketString>::Type       TBucketStringValue;
    typedef typename Value<TSuperBucketString>::Type  TSuperBucketStringValue;

    unsigned bitsPerBucket = BitsPerValue<typename Value<TBitString>::Type>::VALUE;

    TBitString const & bitString = rankSupportBitString.bitString;
    TBucketString const & bucketString = rankSupportBitString.bucketString;
    TSuperBucketString const & superBucketString = rankSupportBitString.superBucketString;

    stream << "  ";
    for (TBitStringValue i = 0; i < length(bitString); i++)
    {
        printBits(stream, bitString[i], bitsPerBucket);
    }
    stream << std::endl;

    for (TBucketStringValue i = 0; i < length(bucketString); i++)
    {
        stream << bucketString[i] << " ";
    }
    stream << std::endl;

    for (TSuperBucketStringValue i = 0; i < length(superBucketString); i++)
    {
        stream << superBucketString[i] << " ";
    }
    return stream;
    //	return print(stream, rankSupportBitString);
}

template <typename TSpec, typename TSize>
inline void reserve(RankSupportBitString<TSpec> & rankSupportBitString, TSize size)
{
    //rankSupportBitString.length = 0;

    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type                TBitString;
    typedef typename Value<TBitString>::Type TBitStringValue;
    unsigned bitsPerBucket = BitsPerValue<TBitStringValue>::VALUE;
    unsigned long long numberOfBuckets;
    (size) ? numberOfBuckets = size / bitsPerBucket + 1 : numberOfBuckets = 0;

    resize(rankSupportBitString.bitString, numberOfBuckets, 0);
    resize(rankSupportBitString.bucketString, numberOfBuckets, 0);
    resize(rankSupportBitString.superBucketString, numberOfBuckets / bitsPerBucket + 1, 0);
}

template <typename TSpec, typename TSize>
inline void resize(RankSupportBitString<TSpec> & rankSupportBitString, TSize size)
{
    reserve(rankSupportBitString, size);
    rankSupportBitString.length = size;
}

template <typename TSpec, typename TSize, typename TValue>
inline void resize(RankSupportBitString<TSpec> & rankSupportBitString, TSize size, TValue value)
{
    rankSupportBitString.length = size;

    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type    TBitString;
    typedef typename Value<TBitString>::Type TBitStringValue;
    unsigned bitsPerBucket = BitsPerValue<TBitStringValue>::VALUE;
    unsigned long long numberOfBuckets;
    (size) ? numberOfBuckets = size / bitsPerBucket + 1 : numberOfBuckets = 0;

    resize(rankSupportBitString.bitString, numberOfBuckets, value);
    resize(rankSupportBitString.bucketString, numberOfBuckets, value);
    resize(rankSupportBitString.superBucketString, numberOfBuckets / bitsPerBucket + 1, value);
}

template <typename TSpec>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type>::Type
length(RankSupportBitString<TSpec> & rankSupportBitString)
{
    return rankSupportBitString.length;
}

template <typename TSpec>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type>::Type
length(RankSupportBitString<TSpec> const & rankSupportBitString)
{
    return rankSupportBitString.length;
}

template <typename TSpec>
inline void clear(RankSupportBitString<TSpec> & rankSupportBitString)
{
    clear(rankSupportBitString.bitString);
    clear(rankSupportBitString.bucketString);
    clear(rankSupportBitString.superBucketString);
}

template <typename TPos>
inline bool getBit(String<bool> & bitString, TPos & pos)
{
    return bitString[pos];
}

template <typename TSpec, typename TPos>
inline bool getBit(RankSupportBitString<TSpec> & bitString, TPos pos)
{
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type                TBitString;
    typedef typename Value<TBitString>::Type TValue;
    TValue bitsPerValue = BitsPerValue<TValue>::VALUE;
    TValue one = 1;
    TValue shiftValue = pos & (bitsPerValue - 1);
    return (bitString.bitString[(pos / bitsPerValue)] >> shiftValue) & one;
}

template <typename TSpec, typename TPos>
inline bool getBit(RankSupportBitString<TSpec> const & bitString, TPos pos)
{
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type                TBitString;
    typedef typename Value<TBitString>::Type TValue;
    TValue bitsPerValue = BitsPerValue<TValue>::VALUE;
    TValue one = 1;
    TValue shiftValue = pos & (bitsPerValue - 1);
    return (bitString.bitString[(pos / bitsPerValue)] >> shiftValue) & one;
}

template <typename TIndex, typename TPos, typename TBit>
inline void setBit(String<bool> & bitString, TPos pos, TBit setBit)
{
    bitString[pos] = setBit;
}

template <typename TSpec, typename TPos, typename TBit>
inline void setBit(RankSupportBitString<TSpec> & bitString, TPos pos, TBit setBit)
{
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type                TBitString;
    typedef typename Value<TBitString>::Type TValue;
    TValue const bitsPerValue = BitsPerValue<typename Value<TBitString>::Type>::VALUE;
    TValue const bucketPos = pos / bitsPerValue;
    TValue const posInBucket = pos % bitsPerValue;
    TValue const one = 1;
    TValue const shiftValue = one << posInBucket;
    if (!setBit)
    {
        bitString.bitString[bucketPos] &= ~(shiftValue);
        return;
    }
    bitString.bitString[bucketPos] |= shiftValue;
}
template <typename TSpec, typename TBit>
inline void appendBit(RankSupportBitString<TSpec> & rankSupportBitString, TBit bit)
{
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type                TBitString;
    typedef typename Value<TBitString>::Type TValue;
    unsigned const bitsPerValue = BitsPerValue<typename Value<TBitString>::Type>::VALUE;


	if((length(rankSupportBitString.bitString) * bitsPerValue) <= (length(rankSupportBitString)))
	{
		reserve(rankSupportBitString, 2 * (length(rankSupportBitString) + 1));
	}
    setBit(rankSupportBitString, length(rankSupportBitString), bit);
    ++rankSupportBitString.length;
}

template <typename TSpec, typename TBit>
inline void appendBitOnly(RankSupportBitString<TSpec> & rankSupportBitString, TBit bit)
{
    setBit(rankSupportBitString, length(rankSupportBitString), bit);
    ++rankSupportBitString.length;
}

template <typename TValue>
inline unsigned getRankInBucket(const TValue value)
{
    return getRankInBucket(value, typename Eval < (BitsPerValue<TValue>::VALUE > 32) > ::Type());
}

template <typename TValue>
inline unsigned getRankInBucket(const TValue value, False)
{
    return __builtin_popcountl(static_cast<int32_t>(value));
}

template <typename TValue>
inline unsigned getRankInBucket(const TValue value, True)
{
    return __builtin_popcountll(static_cast<int64_t>(value));
}

template <typename TValue, typename TPos>
inline TValue getRankInBucket(const String<TValue> & bitString, const TPos pos)
{

    unsigned short const bitsPerValue = BitsPerValue<TValue>::VALUE;
    TValue const one = -1;
   // TValue const mask = one >> (bitsPerValue - (pos % bitsPerValue) - 1);
    TValue const mask = one >> (bitsPerValue - (pos & (bitsPerValue - 1)) - 1);
    return getRankInBucket(bitString[pos / bitsPerValue] & mask);
}

template <typename TSpec>
inline void printBits(const TSpec entrie, const int blocks);

template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type>::Type
getRank(const RankSupportBitString<TSpec> & rankSupportBitString, const TPos pos)
{
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type 		TSuperBucketString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type                TBitString;

    typename Value<TSuperBucketString>::Type sum = 0;
    unsigned bitsPerBucket = BitsPerValue<typename Value<TBitString>::Type>::VALUE;
    unsigned long bucketPos = pos / bitsPerBucket;
    if (bucketPos)
    {
        sum += rankSupportBitString.bucketString[bucketPos];
        unsigned superBucketPos = (bucketPos - 1) / bitsPerBucket;
        if (superBucketPos)
        {
            --superBucketPos;
            sum += rankSupportBitString.superBucketString[superBucketPos];
        }
    }

    sum += getRankInBucket(rankSupportBitString.bitString, pos);
    return sum;
}

template <typename TBitString>
inline void completeRankSupportBitString(TBitString & bitString){}

template <typename TSpec>
//, typename TRankSupportBitString>
inline void completeRankSupportBitString(RankSupportBitString<TSpec> & bitString)
{
    if (length(bitString))
    {
        typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type                TBitString;
        typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBucketString>::Type         TBucketString;
        typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type        TSuperBucketString;

        typedef typename Value<TBucketString>::Type TBucketValue;
        typedef typename Value<TSuperBucketString>::Type TSuperBucketValue;

        resize(bitString, length(bitString));

        unsigned bitsPerBucket = BitsPerValue<typename Value<TBitString>::Type>::VALUE;
        TSuperBucketValue superBucketCounter = 0;

        TBucketValue tempSum = 0;
        TBucketValue bucketSum = 0;
        TSuperBucketValue superBucketSum = 0;

        TBitString & bitString_ = bitString.bitString;
        TBucketString & bucketString = bitString.bucketString;
        TSuperBucketString & superBucketString = bitString.superBucketString;

        typedef typename Value<TSuperBucketString>::Type TSize;
        for (TSize i = 0; i < length(bitString_) - 1; i++)
        {
            tempSum = getRankInBucket(bitString_[i]);
            bucketSum += tempSum;
            bucketString[i + 1] = bucketSum;
            if (!((i + 1) & (bitsPerBucket - 1)))
            {
                superBucketSum += bucketSum;
                superBucketString[superBucketCounter] = superBucketSum;
                bucketSum = 0;
                ++superBucketCounter;
            }
        }
    }
}

//Manuel Forwards
template <typename TText, typename TSpec>
struct WaveletTree;

template <typename TChar, typename TPointer, typename TSpec>
struct WaveletTreeStructure;

struct FibreSplitValues_;
typedef Tag<FibreSplitValues_> const FibreSplitValues;
//template <typename TText>
//struct WaveletTreeStructure;

struct FibreBitStrings_;
typedef Tag<FibreBitStrings_> const FibreBitStrings;

//fills ONE bit string with 0 for characters that will appear in the left branch and
//1 for characters that will appear in the right branch
//template < typename TBitString, typename TCharacterValue, typename TPosInSubTree, typename TText, typename TWaveletTreeSpec >//, typename TRankSupportBitString>
template <typename TBitString, typename TText>
//, typename TRankSupportBitString>
inline void fillBitString(
    const unsigned lowestValue,
    const unsigned splitValue,
    const unsigned highestValue,
    TBitString & bitString,        //TRankSupportBitString &counterBitString,
    const TText & text)
{
    unsigned short character;
    unsigned long long pos = 0;
    if (length(bitString) == 0)
    {
        return;
    }
    for (unsigned i = 0; i < length(text); ++i)
    {
        character = text[i];
        if (character >= lowestValue && character <= highestValue)
        {
            if (character >= splitValue)
            {
                setBit(bitString, pos, 1);
            }
            ++pos;
        }
    }
    completeRankSupportBitString(bitString);
}

template <typename TCharacterValue, typename TWaveletTreeSpec, typename TText>
//, typename TRankSupportBitString>
inline void fillBitString(
    WaveletTree<TText, TWaveletTreeSpec> & waveletTree,
    typename Iterator<typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Type>::Type & iter,
    const TText & text,
    const TCharacterValue lowerBound,
    const TCharacterValue upperBound)
{
    TCharacterValue lowestValue = lowerBound;
    TCharacterValue splitValue = iter.waveletTreeStructure->treeNodes[iter.position].i1;
    TCharacterValue highestValue = upperBound;
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreBitStrings>::Type TBitStrings;
    typedef typename Value<TBitStrings>::Type TBitString;
    TBitString & bitString = waveletTree.bitStrings[iter.position];

    TCharacterValue character;
    unsigned long long pos = 0;
    //	TBitString &bitString = tree.bitStrings[treePos];
    if (length(bitString) == 0)
    {
        return;
    }
    for (unsigned i = 0; i < length(text); ++i)
    {
        character = text[i];
        if ((character >= lowestValue) && (character <= highestValue))
        {
            if (character >= splitValue)
            {
                setBit(bitString, pos, 1);
            }
            ++pos;
        }
    }
    completeRankSupportBitString(bitString);
}

template <typename TText>
inline unsigned short nearestPowOfTwo(TText & text)
{
    unsigned l = length(text);
    --l;
    for (unsigned i = 1; i <= 64; i *= 2)
    {
        l |= (l >> i);
    }
    ++l;
    return l;
}

template <typename TSpec>
inline bool open(
    RankSupportBitString<TSpec> & string,
    const char * fileName,
    int openMode)
{
    typedef typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type>::Type TValue;
    String<TValue> lengthString;
    resize(lengthString, 1);

    String<char> name;
    name = fileName;    append(name, ".bit");
    if (!open(getFibre(string, FibreRankSupportBitString()), toCString(name), openMode))
        return false;

    name = fileName;    append(name, ".bucket");    open(getFibre(string, FibreRankSupportBucketString()), toCString(name), openMode);
    name = fileName;    append(name, ".sbucket");   open(getFibre(string, FibreRankSupportSuperBucketString()), toCString(name), openMode);
    name = fileName;    append(name, ".length");    open(lengthString, toCString(name), openMode);
    string.length = lengthString[0];
    return true;
}

template <typename TSpec>
inline bool open(
    RankSupportBitString<TSpec> & string,
    const char * fileName)
{
    return open(string, fileName, OPEN_RDONLY);
}

// ATTENTION:
// This implementation of open doesn't work with external memory StringSets (External<>, MMap<>)
// If you need a persistent external StringSet you have to use a Owner<ConcatDirect<> > StringSet.
template <typename TSpec>
inline bool open(StringSet<RankSupportBitString<TSpec> > & multi, const char * fileName, int openMode)
{
    SEQAN_CHECKPOINT

    typedef typename Size<RankSupportBitString<TSpec> >::Type TSize;
    typedef String<TSize> TSizeString;
    TSizeString sizeString;
    resize(sizeString, length(multi));
    CharString name = fileName;
    name = fileName;    append(name, ".ssize"); open(sizeString, toCString(name), openMode);

    char id[12];     // 2^32 has 10 decimal digits + 1 (0x00)
    unsigned i = 0;
    clear(multi);
    while (true)
    {
        sprintf(id, ".%u", i);
        name = fileName;
        append(name, id);
        {
            resize(multi, i + 1);
            if (!open(multi[i], toCString(name), (openMode & ~OPEN_CREATE) | OPEN_QUIET))
            {
                resize(multi, i);
                break;
            }
        }
        ++i;
    }

    for (TSize i = 0; i < length(multi); ++i)
    {
        multi[i].length = sizeString[i];
    }

    return i > 1;
}

template <typename TSpec>
inline bool open(
    StringSet<RankSupportBitString<TSpec> > & strings,
    const char * fileName)
{
    return open(strings, fileName, OPEN_RDONLY);
}

template <typename TSpec>
inline bool save(
    RankSupportBitString<TSpec> const & string,
    const char * fileName,
    int openMode)
{
    typedef typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type>::Type TValue;
    String<TValue> lengthString;
    resize(lengthString, 1);
    lengthString[0] = string.length;
    String<char> name;
    name = fileName;    append(name, ".length");    save(lengthString, toCString(name), openMode);
    name = fileName;    append(name, ".bit");       save(getFibre(string, FibreRankSupportBitString()), toCString(name), openMode);
    name = fileName;    append(name, ".bucket");    save(getFibre(string, FibreRankSupportBucketString()), toCString(name), openMode);
    name = fileName;    append(name, ".sbucket");   save(getFibre(string, FibreRankSupportSuperBucketString()), toCString(name), openMode);
    return true;
}

template <typename TSpec>
inline bool save(
    RankSupportBitString<TSpec> const & string,
    const char * fileName)
{
    return save(string, fileName, DefaultOpenMode<RankSupportBitString<TSpec> >::VALUE);
}

template <typename TSpec>
inline bool save(StringSet<RankSupportBitString<TSpec> > const & multi, const char * fileName, int openMode)
{
    SEQAN_CHECKPOINT

    typedef typename Size<RankSupportBitString<TSpec> >::Type TSize;
    typedef String<TSize> TSizeString;
    TSizeString sizeString;
    resize(sizeString, length(multi));
    for (TSize i = 0; i < length(multi); ++i)
    {
        sizeString[i] = length(multi[i]);
    }

    CharString name;
    name = fileName;    append(name, ".ssize"); save(sizeString, toCString(name), openMode);
    if (length(multi) == 0) return true;

    char id[12];     // 2^32 has 10 decimal digits + 2 ('.' and 0x00)
    for (unsigned i = 0; i < length(multi); ++i)
    {
        sprintf(id, ".%u", i);
        name = fileName;
        append(name, &(id[0]));
        if (!save(multi[i], toCString(name), openMode))
            return false;
    }
    return true;
}

template <typename TSpec>
inline bool save(
    StringSet<RankSupportBitString<TSpec> > const & strings,
    const char * fileName)
{
    return save(strings, fileName, OPEN_RDONLY);
}

}


#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_BITSTRING_H_
