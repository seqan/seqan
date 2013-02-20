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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// State: Tested but not used in real code.
// ==========================================================================
// Implementation of the LongWord class.  It provides an almost
// transparent interface to arbitrarily long words.  In contrast to
// std::bitset, it provides a specialization that allows to set the
// number of bits in the virtual word at runtime.
// ==========================================================================

// TODO(holtgrew): This should probably not be called LongWord, maybe BitSet/BitVector/BitString is more appropriate?
// TODO(holtgrew): Optimize implementation of operator>>(TWord, 1) and operator<<(TWord, 1).
// TODO(holtgrew): Optimized implementation of value(TWord, index), operator[] always seems to call the proxy version!

#ifndef SEQAN_MISC_MISC_LONG_WORD_H_
#define SEQAN_MISC_MISC_LONG_WORD_H_

namespace seqan {

// TODO(holtgrew): Document this.
template <typename TSpec>
struct LongWord;


template <typename TWord, typename TSpec>
struct LongWordBitProxy;


struct _NativeWidth;
typedef Tag<_NativeWidth> NativeWidth;


struct _NativeWideWidth;
typedef Tag<_NativeWideWidth> NativeWideWidth;


struct DynamicWidth_;
typedef Tag<DynamicWidth_> DynamicWidth;


template <unsigned LENGTH>
struct StaticWidth;


template <typename TWord>
struct LongWordBitProxy<TWord, NativeWidth> {
    TWord &_word;
    unsigned _bitIndex;

    LongWordBitProxy(TWord & word, unsigned bitIndex)
            : _word(word), _bitIndex(bitIndex) {
        SEQAN_CHECKPOINT;
    }

    LongWordBitProxy & operator=(unsigned x) {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LT(x, sizeof(unsigned) * 8);
        _word._value = (_word._value ^ (1 << _bitIndex)) | (x << _bitIndex);
        return *this;
    }

    operator unsigned() const {
        return (_word & (1u << _bitIndex)) >> _bitIndex;
    }
};


// Forward declaration.
inline LongWord<NativeWidth> & reset(LongWord<NativeWidth> & x);


template <>
struct LongWord<NativeWidth> {
    typedef LongWord<NativeWidth> TWord;
    unsigned _value;

    LongWord() {
        SEQAN_CHECKPOINT;
        reset(*this);
    }

    LongWord(unsigned const & value) : _value(value) {
        SEQAN_CHECKPOINT;
    }

    // Conversion operators must be member functions.
    operator unsigned() const {
        SEQAN_CHECKPOINT;
        return _value;
    }

    unsigned operator[](unsigned index) const {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LT(index, sizeof(unsigned) * 8);
        return (_value & (1u << index)) >> index;
    }

    LongWordBitProxy<TWord, NativeWidth> operator[](unsigned index) {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LT(index, sizeof(unsigned) * 8);
        return LongWordBitProxy<TWord, NativeWidth>(*this, index);
    }
};


inline LongWord<NativeWidth> & reset(LongWord<NativeWidth> & x) {
    SEQAN_CHECKPOINT;
    x._value = 0u;
    return x;
}


inline size_t length(LongWord<NativeWidth> const &) {
    SEQAN_CHECKPOINT;
    return sizeof(unsigned) * 8;
}


inline LongWord<NativeWidth> & operator>>=(LongWord<NativeWidth> & x, unsigned const shift) {
    SEQAN_CHECKPOINT;
    x._value >>= shift;
    return x;
}


inline LongWord<NativeWidth> & operator<<=(LongWord<NativeWidth> & x, unsigned const shift) {
    SEQAN_CHECKPOINT;
    x._value <<= shift;
    return x;
}


inline LongWord<NativeWidth> & operator&=(LongWord<NativeWidth> & x, LongWord<NativeWidth> const & mask) {
    SEQAN_CHECKPOINT;
    x._value &= mask;
    return x;
}


inline LongWord<NativeWidth> & operator|=(LongWord<NativeWidth> & x, LongWord<NativeWidth> const & mask) {
    SEQAN_CHECKPOINT;
    x._value |= mask;
    return x;
}


inline LongWord<NativeWidth> & operator^=(LongWord<NativeWidth> & x, LongWord<NativeWidth> const & mask) {
    SEQAN_CHECKPOINT;
    x._value ^= mask;
    return x;
}


template <>
struct LongWord<NativeWideWidth> {
};


template <typename TWord, unsigned LENGTH>
struct LongWordBitProxy<TWord, StaticWidth<LENGTH> > {
    TWord &_word;
    unsigned _bitIndex;

    LongWordBitProxy(TWord & word, unsigned bitIndex)
            : _word(word), _bitIndex(bitIndex) {
        SEQAN_CHECKPOINT;
    }

    LongWordBitProxy & operator=(unsigned x) {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LT(x, LENGTH);
        unsigned & block = _word._data[TWord::UNSIGNED_COUNT - 1 - _bitIndex / BitsPerValue<unsigned>::VALUE];
        unsigned bitIndex = _bitIndex % BitsPerValue<unsigned>::VALUE;
        block = (block ^ (1 << bitIndex)) | (x << bitIndex);
        return *this;
    }

    operator unsigned() const {
        SEQAN_CHECKPOINT;
        unsigned const & block = _word._data[TWord::UNSIGNED_COUNT - 1 - _bitIndex / BitsPerValue<unsigned>::VALUE];
        unsigned bitIndex = _bitIndex % BitsPerValue<unsigned>::VALUE;
//         std::cerr << "LongWordBitProxy<... StaticWidth<" << LENGTH << "> >" << std::endl;
//         std::cerr << "  _word      = " << _word << std::endl;
//         std::cerr << "  _bitIndex  = " << _bitIndex << std::endl;
//         std::cerr << "  block      = " << block << std::endl;
//         std::cerr << "  bitIndex   = " << bitIndex << std::endl;
//         std::cerr << "  (block & (1 << bitIndex)) >> bitIndex = " << ((block & (1 << bitIndex)) >> bitIndex) << std::endl;
//         std::cerr << std::endl;
        return (block & (1 << bitIndex)) >> bitIndex;
    }
};


template <unsigned LENGTH>
struct LongWord<StaticWidth<LENGTH> > {
    typedef LongWord<StaticWidth<LENGTH> > TWord;
    enum TDummy { UNSIGNED_COUNT = (LENGTH / BitsPerValue<unsigned>::VALUE + ((LENGTH % BitsPerValue<unsigned>::VALUE > 0) ? 1 : 0)) };
    unsigned _data[UNSIGNED_COUNT];

    LongWord() {
        SEQAN_CHECKPOINT;
        reset(*this);
    }

    LongWord(LongWord const & other) {
        SEQAN_CHECKPOINT;
        for (size_t i = 0; i < UNSIGNED_COUNT; ++i)
            _data[i] = other._data[i];
    }

    LongWord & operator=(LongWord const & other) {
        SEQAN_CHECKPOINT;
        for (size_t i = 0; i < UNSIGNED_COUNT; ++i)
            _data[i] = other._data[i];
        return *this;
    }

    unsigned operator[](unsigned index) const {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LT(index, LENGTH);
        unsigned const & block = _data[UNSIGNED_COUNT - 1 - index / BitsPerValue<unsigned>::VALUE];
        unsigned bitIndex = index % BitsPerValue<unsigned>::VALUE;
        return (block & (1 << bitIndex)) >> bitIndex;
    }

    LongWordBitProxy<TWord, StaticWidth<LENGTH> > operator[](unsigned index) {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LT(index, LENGTH);
        return LongWordBitProxy<TWord, StaticWidth<LENGTH> >(*this, index);
    }
};


template <unsigned LENGTH>
inline LongWord<StaticWidth<LENGTH> > & reset(LongWord<StaticWidth<LENGTH> > & a) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    for (size_t i = 0; i < TLongWord::UNSIGNED_COUNT; ++i)
        a._data[i] = 0u;
    return a;
}


template <unsigned LENGTH>
size_t length(LongWord<StaticWidth<LENGTH> > const &) {
    SEQAN_CHECKPOINT;
    return LENGTH;
}


template <unsigned LENGTH>
bool operator==(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i) {
        if (a._data[i] != b._data[i]) {
            return false;
        }
    }
    return true;    
}


template <unsigned LENGTH>
bool operator!=(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Roll out loop?
//     for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i) {
//         if (a._data[i] == b._data[i]) {
//             return false;
//         }
//     }
//     return true;
    return !(a == b);
}


template <unsigned LENGTH>
bool operator<=(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i) {
        if (a._data[TLongWord::UNSIGNED_COUNT - 1 - i] > b._data[TLongWord::UNSIGNED_COUNT - 1 - i])
            return false;
    }
    return true;    
}


template <unsigned LENGTH>
bool operator>=(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i) {
        if (a._data[TLongWord::UNSIGNED_COUNT - 1 - i] <= b._data[TLongWord::UNSIGNED_COUNT - 1 - i])
            return false;
    }
    return true;    
}


template <unsigned LENGTH>
bool operator<(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    return !(a >= b);
//     typedef LongWord<StaticWidth<LENGTH> > TLongWord;
//     // TODO(holtgrew): Roll out loop?
//     for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i) {
//         if (a._data[i] >= b._data[i])
//             return false;
//     }
//     return true;    
}


template <unsigned LENGTH>
bool operator>(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    return !(a <= b);
//     typedef LongWord<StaticWidth<LENGTH> > TLongWord;
//     // TODO(holtgrew): Roll out loop?
//     for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i) {
//         if (a._data[i] <= b._data[i])
//             return false;
//     }
//     return true;    
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > operator|(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    TLongWord result;
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i)
        result._data[i] = a._data[i] | b._data[i];
    return result;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > & operator|=(LongWord<StaticWidth<LENGTH> > & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    TLongWord result;
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i)
        a._data[i] |= b._data[i];
    return a;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > operator&(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    TLongWord result;
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i)
        result._data[i] = a._data[i] & b._data[i];
    return result;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > & operator&=(LongWord<StaticWidth<LENGTH> > & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    TLongWord result;
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i)
        a._data[i] &= b._data[i];
    return a;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > operator^(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    TLongWord result;
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i)
        result._data[i] = a._data[i] ^ b._data[i];
    return result;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > & operator^=(LongWord<StaticWidth<LENGTH> > & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    TLongWord result;
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i)
        a._data[i] ^= b._data[i];
    return a;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > operator>>(LongWord<StaticWidth<LENGTH> > const & a, unsigned shift) {
    SEQAN_CHECKPOINT;
    LongWord<StaticWidth<LENGTH> > x(a);
    x >>= shift;
    return x;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > & operator>>=(LongWord<StaticWidth<LENGTH> > & a, unsigned shift) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // Clear if shift is larger than the number of bits.
    if (shift > LENGTH)
        return reset(a);
    // Do nothing if shift is zero.
    if (shift == 0u)
        return a;

    // Otherwise, perform shift operation.
    size_t const last = TLongWord::UNSIGNED_COUNT - 1;
    size_t const blockShift = shift / BitsPerValue<unsigned>::VALUE;    // Number of blocks to shift values from.
    size_t const indexInBlock = shift % BitsPerValue<unsigned>::VALUE;
    if (indexInBlock != 0) {
        // More common case: shift is not multiple of bit count of unsigned.
        unsigned const x = BitsPerValue<unsigned>::VALUE - indexInBlock;
        for (size_t i = last - blockShift; i < last; ++i)
            a._data[i + blockShift] = (a._data[i] >> indexInBlock) | (a._data[i + 1] << x);
        a._data[last - blockShift] = a._data[0] >> indexInBlock;
    } else {
        // Less common case: We can simply copy these unsigneds.
        for (size_t i = blockShift; i <= last; ++i)
            a._data[i + blockShift] = a._data[0];
    }
    std::fill_n(&(a._data[TLongWord::UNSIGNED_COUNT - 1 - blockShift]), blockShift, 0u);
    return a;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > operator<<(LongWord<StaticWidth<LENGTH> > const & a, unsigned shift) {
    SEQAN_CHECKPOINT;
    LongWord<StaticWidth<LENGTH> > x(a);
    x <<= shift;
    return x;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > & operator<<=(LongWord<StaticWidth<LENGTH> > & a, unsigned shift) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // Clear if shift is larger than the number of bits.
    if (shift > LENGTH)
        return reset(a);
    // Do nothing if shift is zero.
    if (shift == 0u)
        return a;

    // Otherwise, perform shift operation.
    size_t const last = TLongWord::UNSIGNED_COUNT - 1;
    size_t const blockShift = shift / BitsPerValue<unsigned>::VALUE;    // Number of blocks to shift values from.
    size_t const indexInBlock = shift % BitsPerValue<unsigned>::VALUE;
    if (indexInBlock != 0) {
        // More common case: shift is not multiple of bit count of unsigned.
        unsigned const x = BitsPerValue<unsigned>::VALUE - indexInBlock;
        for (size_t i = last - blockShift; i > 0; --i)
            a._data[i + blockShift] = (a._data[i] << indexInBlock) | (a._data[i - 1] >> x);
        a._data[blockShift] = a._data[0] << indexInBlock;
    } else {
        // Less common case: We can simply copy these unsigneds.
        for (size_t i = last - blockShift; i > 0; --i)
            a._data[i + blockShift] = a._data[0];
        a._data[blockShift] = a._data[0];
    }
    std::fill_n(&(a._data[0]), blockShift, 0u);
    return a;
}


template <typename TWord>
struct LongWordBitProxy<TWord, DynamicWidth> {
    TWord &_word;
    unsigned _bitIndex;

    LongWordBitProxy(TWord & word, unsigned bitIndex)
            : _word(word), _bitIndex(bitIndex) {
        SEQAN_CHECKPOINT;
    }

    LongWordBitProxy & operator=(unsigned x) {
        SEQAN_CHECKPOINT;
        unsigned & block = _word._data[_word._unsignedCount() - 1 - _bitIndex / BitsPerValue<unsigned>::VALUE];
        unsigned bitIndex = _bitIndex % BitsPerValue<unsigned>::VALUE;
        block = (block ^ (1 << bitIndex)) | (x << bitIndex);
        return *this;
    }

    operator unsigned() const {
        SEQAN_CHECKPOINT;
        unsigned const & block = _word._data[_word._unsignedCount() - 1 - _bitIndex / BitsPerValue<unsigned>::VALUE];
        unsigned bitIndex = _bitIndex % BitsPerValue<unsigned>::VALUE;
//         std::cerr << "LongWordBitProxy<... StaticWidth<" << LENGTH << "> >" << std::endl;
//         std::cerr << "  _word      = " << _word << std::endl;
//         std::cerr << "  _bitIndex  = " << _bitIndex << std::endl;
//         std::cerr << "  block      = " << block << std::endl;
//         std::cerr << "  bitIndex   = " << bitIndex << std::endl;
//         std::cerr << "  (block & (1 << bitIndex)) >> bitIndex = " << ((block & (1 << bitIndex)) >> bitIndex) << std::endl;
//         std::cerr << std::endl;
        return (block & (1 << bitIndex)) >> bitIndex;
    }
};


// Forward Declarations.
inline LongWord<DynamicWidth> & reset(LongWord<DynamicWidth> & a);


template <>
struct LongWord<DynamicWidth> {
    typedef LongWord<DynamicWidth> TWord;
    unsigned _length;  // Number of bits in the word.
    unsigned * _data;

    explicit
    LongWord(unsigned const & length)
            : _length(length) {
        SEQAN_CHECKPOINT;
        _data = new unsigned[_unsignedCount()];
        reset(*this);
    }

    LongWord(LongWord const & other) {
        SEQAN_CHECKPOINT;
        _length = other._length;
        _data = new unsigned[_unsignedCount()];
        for (size_t i = 0; i < _unsignedCount(); ++i)
            _data[i] = other._data[i];
    }

    ~LongWord() {
        delete [] _data;
    }

    LongWord & operator=(LongWord const & other) {
        SEQAN_CHECKPOINT;
        if (_length != other._length) {
            _length = other._length;
            _data = new unsigned[_unsignedCount()];
        }
        for (size_t i = 0; i < _unsignedCount(); ++i)
            _data[i] = other._data[i];
        return *this;
    }

    unsigned operator[](unsigned index) const {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LT(index, _length);
        unsigned const & block = _data[_unsignedCount() - 1 - index / BitsPerValue<unsigned>::VALUE];
        unsigned bitIndex = index % BitsPerValue<unsigned>::VALUE;
        return (block & (1 << bitIndex)) >> bitIndex;
    }

    LongWordBitProxy<TWord, DynamicWidth> operator[](unsigned index) {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LT(index, _length);
        return LongWordBitProxy<TWord, DynamicWidth>(*this, index);
    }

    unsigned _unsignedCount() const {
        SEQAN_CHECKPOINT;
        return _length / BitsPerValue<unsigned>::VALUE + ((_length % BitsPerValue<unsigned>::VALUE > 0) ? 1 : 0);
    }
};


inline LongWord<DynamicWidth> & reset(LongWord<DynamicWidth> & a) {
    SEQAN_CHECKPOINT;
    for (size_t i = 0; i < a._unsignedCount(); ++i)
        a._data[i] = 0u;
    return a;
}


size_t length(LongWord<DynamicWidth> const & a) {
    SEQAN_CHECKPOINT;
    return a._length;
}


bool operator==(LongWord<DynamicWidth> const & a, LongWord<DynamicWidth> const & b) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_EQ(a._length, b._length);
    // TODO(holtgrew): Roll out loop?
    for (unsigned i = 0; i < a._unsignedCount(); ++i) {
        if (a._data[i] != b._data[i]) {
            return false;
        }
    }
    return true;    
}


bool operator!=(LongWord<DynamicWidth> const & a, LongWord<DynamicWidth> const & b) {
    SEQAN_CHECKPOINT;
    return !(a == b);
}


bool operator<=(LongWord<DynamicWidth> const & a, LongWord<DynamicWidth> const & b) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_EQ(a._length, b._length);
    // TODO(holtgrew): Roll out loop?
    for (unsigned i = 0; i < a._unsignedCount(); ++i) {
        if (a._data[a._unsignedCount() - 1 - i] > b._data[a._unsignedCount() - 1 - i])
            return false;
    }
    return true;    
}


bool operator>=(LongWord<DynamicWidth> const & a, LongWord<DynamicWidth> const & b) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_EQ(a._length, b._length);
    // TODO(holtgrew): Roll out loop?
    for (unsigned i = 0; i < a._unsignedCount(); ++i) {
        if (a._data[a._unsignedCount() - 1 - i] <= b._data[a._unsignedCount() - 1 - i])
            return false;
    }
    return true;    
}


bool operator<(LongWord<DynamicWidth> const & a, LongWord<DynamicWidth> const & b) {
    SEQAN_CHECKPOINT;
    return !(a >= b);
}


bool operator>(LongWord<DynamicWidth> const & a, LongWord<DynamicWidth> const & b) {
    SEQAN_CHECKPOINT;
    return !(a <= b);
}


LongWord<DynamicWidth> & operator|=(LongWord<DynamicWidth> & a, LongWord<DynamicWidth> const & b) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_EQ(a._length, b._length);
    // TODO(holtgrew): Roll out loop?
    for (unsigned i = 0; i < a._unsignedCount(); ++i)
        a._data[i] |= b._data[i];
    return a;
}


LongWord<DynamicWidth> operator|(LongWord<DynamicWidth> const & a, LongWord<DynamicWidth> const & b) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_EQ(a._length, b._length);
    LongWord<DynamicWidth> result(a);
    return result |= b;
}


LongWord<DynamicWidth> & operator&=(LongWord<DynamicWidth> & a, LongWord<DynamicWidth> const & b) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_EQ(a._length, b._length);
    // TODO(holtgrew): Roll out loop?
    for (unsigned i = 0; i < a._unsignedCount(); ++i)
        a._data[i] &= b._data[i];
    return a;
}


LongWord<DynamicWidth> operator&(LongWord<DynamicWidth> const & a, LongWord<DynamicWidth> const & b) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_EQ(a._length, b._length);
    LongWord<DynamicWidth> result(a);
    return result &= b;
}


LongWord<DynamicWidth> & operator^=(LongWord<DynamicWidth> & a, LongWord<DynamicWidth> const & b) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_EQ(a._length, b._length);
    // TODO(holtgrew): Roll out loop?
    for (unsigned i = 0; i < a._unsignedCount(); ++i)
        a._data[i] ^= b._data[i];
    return a;
}


LongWord<DynamicWidth> operator^(LongWord<DynamicWidth> const & a, LongWord<DynamicWidth> const & b) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_EQ(a._length, b._length);
    LongWord<DynamicWidth> result(a);
    return result ^= b;
}


LongWord<DynamicWidth> & operator>>=(LongWord<DynamicWidth> & a, unsigned shift) {
    SEQAN_CHECKPOINT;
    // Clear if shift is larger than the number of bits.
    if (shift > a._length)
        return reset(a);
    // Do nothing if shift is zero.
    if (shift == 0u)
        return a;

    // Otherwise, perform shift operation.
    size_t const last = a._unsignedCount() - 1;
    size_t const blockShift = shift / BitsPerValue<unsigned>::VALUE;    // Number of blocks to shift values from.
    size_t const indexInBlock = shift % BitsPerValue<unsigned>::VALUE;
    if (indexInBlock != 0) {
        // More common case: shift is not multiple of bit count of unsigned.
        unsigned const x = BitsPerValue<unsigned>::VALUE - indexInBlock;
        for (size_t i = last - blockShift; i < last; ++i)
            a._data[i + blockShift] = (a._data[i] >> indexInBlock) | (a._data[i + 1] << x);
        a._data[last - blockShift] = a._data[0] >> indexInBlock;
    } else {
        // Less common case: We can simply copy these unsigneds.
        for (size_t i = blockShift; i <= last; ++i)
            a._data[i + blockShift] = a._data[0];
    }
    std::fill_n(&(a._data[a._unsignedCount() - 1 - blockShift]), blockShift, 0u);
    return a;
}


LongWord<DynamicWidth> operator>>(LongWord<DynamicWidth> const & a, unsigned shift) {
    SEQAN_CHECKPOINT;
    LongWord<DynamicWidth> result(a);
    return result >>= shift;
}


LongWord<DynamicWidth> & operator<<=(LongWord<DynamicWidth> & a, unsigned shift) {
    SEQAN_CHECKPOINT;
    // Clear if shift is larger than the number of bits.
    if (shift > a._length)
        return reset(a);
    // Do nothing if shift is zero.
    if (shift == 0u)
        return a;

    // Otherwise, perform shift operation.
    size_t const last = a._unsignedCount() - 1;
    size_t const blockShift = shift / BitsPerValue<unsigned>::VALUE;    // Number of blocks to shift values from.
    size_t const indexInBlock = shift % BitsPerValue<unsigned>::VALUE;
    if (indexInBlock != 0) {
        // More common case: shift is not multiple of bit count of unsigned.
        unsigned const x = a._unsignedCount() - indexInBlock;
        for (size_t i = last - blockShift; i > 0; --i)
            a._data[i + blockShift] = (a._data[i] << indexInBlock) | (a._data[i - 1] >> x);
        a._data[blockShift] = a._data[0] << indexInBlock;
    } else {
        // Less common case: We can simply copy these unsigneds.
        for (size_t i = last - blockShift; i > 0; --i)
            a._data[i + blockShift] = a._data[0];
        a._data[blockShift] = a._data[0];
    }
    std::fill_n(&(a._data[0]), blockShift, 0u);
    return a;
}


LongWord<DynamicWidth> operator<<(LongWord<DynamicWidth> const & a, unsigned shift) {
    SEQAN_CHECKPOINT;
    LongWord<DynamicWidth> x(a);
    x <<= shift;
    return x;
}

}  // namespace seqan


namespace std {

template <typename TStream, unsigned LENGTH>
TStream & operator<<(TStream & stream, seqan::LongWord<seqan::StaticWidth<LENGTH> > const & a) {
    SEQAN_CHECKPOINT;
    for (size_t i = 0; i < LENGTH; ++i) {
        unsigned x = a[LENGTH - 1 - i];
        stream << x;
    }
    typedef seqan::LongWord<seqan::StaticWidth<LENGTH> > TLongWord;
    stream << "  ( ";
    for (size_t i = 0; i < TLongWord::UNSIGNED_COUNT; ++i) {
        stream << a._data[i] << " ";
    }
    stream << ")";
    return stream;
}


template <typename TStream>
TStream & operator<<(TStream & stream, seqan::LongWord<seqan::DynamicWidth> const & a) {
    SEQAN_CHECKPOINT;
    for (size_t i = 0; i < a._length; ++i) {
        unsigned x = a[a._length - 1 - i];
        stream << x;
    }
    //typedef seqan::LongWord<seqan::DynamicWidth> TLongWord;
    for (size_t i = 0; i < a._unsignedCount(); ++i) {
        stream << a._data[i] << " ";
    }
    stream << ")";
    return stream;
}

}  // namespace std

#endif  // SEQAN_MISC_MISC_LONG_WORD_H_
