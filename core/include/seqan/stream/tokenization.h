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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Adaptions for std::ios streams.
// ==========================================================================

#ifndef SEQAN_STREAM_TOKENIZATION_H_
#define SEQAN_STREAM_TOKENIZATION_H_

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TValue>
struct IsInAlphabet
{
    template <typename TInValue>
    bool operator() (TInValue inVal) const
    {
        TValue val = inVal;
        return convert<TInValue>(val) == toUpperValue(inVal);
    }

    bool operator() (TValue) const
    {
        return true;
    }
};

template <char FIRST_CHAR, char LAST_CHAR>
struct IsInRange
{
    template <typename TValue>
    bool operator() (TValue val) const
    {
        return FIRST_CHAR <= val && val <= LAST_CHAR;
    }
};

template <char VALUE>
struct EqualsChar
{
    template <typename TValue>
    bool operator() (TValue val) const
    {
        return val == VALUE;
    }
};

template <typename TFunctor1, typename TFunctor2>
struct OrFunctor
{
    TFunctor1 func1;
    TFunctor2 func2;

    OrFunctor()
    {}

    OrFunctor(TFunctor1 const &func1, TFunctor2 const &func2):
        func1(func1), func2(func2)
    {}

    template <typename TValue>
    bool operator() (TValue val) const
    {
        return func1(val) || func2(val);
    }
};

template <typename TFunctor1, typename TFunctor2>
struct AndFunctor
{
    TFunctor1 func1;
    TFunctor2 func2;

    AndFunctor()
    {}

    AndFunctor(TFunctor1 const &func1, TFunctor2 const &func2):
        func1(func1), func2(func2)
    {}

    template <typename TValue>
    bool operator() (TValue val) const
    {
        return func1(val) && func2(val);
    }
};

template <typename TFunctor>
struct NotFunctor
{
    TFunctor func;

    NotFunctor()
    {}

    NotFunctor(TFunctor const &func):
        func(func)
    {}

    template <typename TValue>
    bool operator() (TValue val) const
    {
        return !func(val);
    }
};

template <typename TIgnoreFunctor, typename TAssertFunctor, typename TException>
struct IgnoreOrAssertFunctor
{
    TIgnoreFunctor ignoreFunc;
    TAssertFunctor assertFunc;
    const char *message;

    IgnoreOrAssertFunctor(const char *message):
        message(message)
    {}

    template <typename TValue>
    bool operator() (TValue val) const
    {
        if (ignoreFunc(val))
            return true;
        if (!assertFunc(val))
            throw TException(message);
    }
};

typedef OrFunctor<EqualsChar<' '>, EqualsChar<'\t'> >           IsBlank;
typedef OrFunctor<EqualsChar<'\n'>, EqualsChar<'\r'> >          IsNewline;
typedef OrFunctor<IsBlank, IsNewline>                           IsWhitespace;
typedef IsInRange<'!', '~'>                                     IsGraph;
typedef OrFunctor<IsInRange<'a', 'z'>, IsInRange<'A', 'Z'> >    IsAlpha;
typedef IsInRange<'0', '9'>                                     IsDigit;
typedef OrFunctor<IsAlpha, IsDigit>                             IsAlphaNum;

// By default, iterators don't support chunking.
template <typename TIterator>
struct Chunk
{
    typedef Nothing Type;
};

// Chunk interface for rooted iterators
template <typename TContainer, typename TValue, typename TSpec>
struct Chunk<Iter<TContainer, AdaptorIterator<TValue*, TSpec> > >
{
    typedef Buffer<TValue, Simple> Type;
};

template <typename TValue, typename TSpec>
struct Chunk<String<TValue, TSpec> >:
    Chunk<typename Iterator<String<TValue, TSpec>, Rooted>::Type> {};

// ============================================================================
// Functions
// ============================================================================

template <class charT, class traits>
struct Value<std::istreambuf_iterator<charT, traits> >
{
    typedef typename std::istreambuf_iterator<charT, traits>::char_type Type;
};

template <class charT, class traits>
struct Position<std::istreambuf_iterator<charT, traits> >
{
    typedef typename std::istreambuf_iterator<charT, traits>::difference_type Type;
};

template <class charT, class traits>
struct Reference<std::istreambuf_iterator<charT, traits> >
{
    typedef typename std::istreambuf_iterator<charT, traits>::char_type Type;
//    typedef typename std::istreambuf_iterator<charT, traits>::reference Type;
};

template <class charT, class traits>
inline bool
atEnd(std::istreambuf_iterator<charT, traits> const &it)
{
    return *it == traits::eof();
}




template <typename TIterator, typename TSize>
inline void
reserveChunk(TIterator &, TSize)
{
}

template <typename TValue, typename TSpec, typename TSize>
inline void
reserveChunk(String<TValue, TSpec> &str, TSize size)
{
    reserve(str, length(str) + size);
}

template <typename TContainer, typename TSpec, typename TSize>
inline void
reserveChunk(Iter<TContainer, TSpec> &iter, TSize size)
{
    typedef Iter<TContainer, TSpec> TIter;

    TContainer &cont = container(iter);
    typename Size<TIter>::Type newCap = length(cont) + size;

    if (newCap <= capacity(cont))
        return;

    typename Position<TIter>::Type pos = position(iter);
    reserve(cont, newCap);
    setPosition(iter, pos);
}



template <typename TIterator, typename TSize>
inline void
advanceChunk(TIterator &iter, TSize size)
{
    iter += size;
}

template <typename TContainer, typename TSpec, typename TSize>
inline void
advanceChunk(Iter<TContainer, TSpec> &iter, TSize size)
{
    typedef Iter<TContainer, TSpec> TIter;

    TContainer &cont = container(iter);
    typename Position<TIter>::Type pos = position(iter);

    iter += size;
    if (pos > length(cont))
        _setLength(cont, pos);
}

// extend target string size
template <typename TValue, typename TSpec, typename TSize>
inline void
advanceChunk(String<TValue, TSpec> &str, TSize size)
{
    _setLength(str, length(str) + size);
}


template <typename TContainer, typename TValue>
inline void
writeValue(TContainer &cont, TValue val)
{
    appendValue(cont, val)
}

template <typename TValue, typename TTraits, typename TValue2>
inline void
writeValue(std::ostreambuf_iterator<TValue, TTraits> &iter, TValue2 val)
{
    *iter = val;
    ++iter;
}

template <typename TContainer, typename TSpec, typename TValue>
inline void
writeValue(Iter<TContainer, TSpec> &iter, TValue val)
{
    typedef Iter<TContainer, TSpec> TIter;

    TContainer &cont = container(iter);
    typename Position<TIter>::Type pos = position(iter);
    typename Size<TIter>::Type len = length(cont);

    if (pos < len)
    {
        assignValue(iter, val);
        ++iter;
    }
    else
    {
        if (pos > len)
            resize(cont, pos - 1)
        appendValue(cont, val);
        setPosition(iter, pos + 1);
    }
}

template <typename TContainer, typename TValue, typename TSpec>
typename Chunk<Iter<TContainer, AdaptorIterator<TValue*, TSpec> > >::Type
getChunk(Iter<TContainer, AdaptorIterator<TValue*, TSpec> > const &rootedIter)
{
    TContainer &cont = container(rootedIter);
    return typename Chunk<Iter<TContainer, AdaptorIterator<TValue*, TSpec> > >::Type(
        hostIterator(rootedIter),
        end(cont, Standard()),
        begin(cont, Standard()) + capacity(cont));
}

template <typename TValue, typename TSpec>
typename Chunk<String<TValue, TSpec> >::Type
getChunk(String<TValue, TSpec> const &cont)
{
    return typename Chunk<String<TValue, TSpec> >::Type(
        end(cont, Standard()),
        end(cont, Standard()),
        begin(cont, Standard()) + capacity(cont));
}


// _readUntil() - non-chunked version
template <typename TTarget, typename TFwdIterator, typename TStopFunctor, typename TIgnoreFunctor, typename TIChunk, typename TOChunk>
inline void
_readUntil(TTarget &target, TFwdIterator &iter, TStopFunctor &stopFunctor, TIgnoreFunctor &ignoreFunctor, TIChunk, TOChunk)
{
    typename Value<TFwdIterator>::Type val;
    for (; !atEnd(iter); ++iter)
    {
        if (stopFunctor(val = *iter))
            return;
        if (!ignoreFunctor(val))
            writeValue(target, val);
    }
}

// _readUntil() - chunked version
template <typename TTarget, typename TFwdIterator, typename TStopFunctor, typename TIgnoreFunctor, typename TIValue, typename TOValue>
inline void
_readUntil(
    TTarget &target,
    TFwdIterator &iter,
    TStopFunctor &stopFunctor,
    TIgnoreFunctor &ignoreFunctor,
    Buffer<TIValue, Simple> *,
    Buffer<TOValue, Simple> *)
{
    Buffer<TIValue, Simple> ichunk;
    Buffer<TOValue, Simple> ochunk;

    for (; !atEnd(iter); )
    {
        ichunk = getChunk(iter);

        // reserve memory for the worst-case
        // TODO(weese):Document worst-case behavior
        reserveChunk(target, ichunk.end - ichunk.begin);
        ochunk = getChunk(end(target, Rooted()));

        TIValue* istart = ichunk.begin;
        TOValue* ostart = ochunk.begin;

        for (; ichunk.begin != ichunk.end; ++ichunk.begin)
        {
            if (stopFunctor(*ichunk.begin))
            {
                iter += ichunk.begin - istart;              // advance input iterator
                advanceChunk(target, ochunk.begin - ostart);     // extend target string size
                return;
            }
            if (!ignoreFunctor(*ichunk.begin))
            {
                // construct values in reserved memory
                valueConstruct(ochunk.begin, getValue(ichunk.begin));
                if (++ochunk.begin == ochunk.resEnd)
                    break;
            }
        }

        iter += ichunk.begin - istart;                      // advance input iterator
        advanceChunk(target, ochunk.begin - ostart);
    }
}

// readUntil()
template <typename TTarget, typename TFwdIterator, typename TStopFunctor, typename TIgnoreFunctor>
inline void
readUntil(TTarget &target, TFwdIterator &iter, TStopFunctor &stopFunctor, TIgnoreFunctor &ignoreFunctor)
{
    typedef typename Chunk<TFwdIterator>::Type*         TIChunk;
    typedef typename Iterator<TTarget, Rooted>::Type    TTargetIter;
    typedef typename Chunk<TTargetIter>::Type*          TOChunk;

    _readUntil(target, iter, stopFunctor, ignoreFunctor, TIChunk(), TOChunk());
}

template <typename TTarget, typename TFwdIterator, typename TStopFunctor, typename TIgnoreFunctor>
inline void
readUntil(TTarget &target, TFwdIterator &iter, TStopFunctor const &stopFunctor, TIgnoreFunctor const &ignoreFunctor)
{
    TStopFunctor stopFunctor_ = stopFunctor;
    TIgnoreFunctor ignoreFunctor_ = ignoreFunctor;
    readUntil(target, iter, stopFunctor_, ignoreFunctor_);
}


// _skipUntil() - non-chunked version
template <typename TFwdIterator, typename TStopFunctor, typename TChunk>
inline void
_skipUntil(TFwdIterator &iter, TStopFunctor &stopFunctor, TChunk)
{
    for (; !atEnd(iter) && !stopFunctor(*iter); ++iter) ;
}

// _skipUntil() - chunked version
template <typename TFwdIterator, typename TStopFunctor, typename TValue>
inline void
_skipUntil(TFwdIterator &iter, TStopFunctor &stopFunctor, Buffer<TValue, Simple>)
{
    typedef typename Value<TFwdIterator>::Type TIValue;
    Buffer<TIValue, Simple> ichunk;

    for (; !atEnd(iter); )
    {
        ichunk = getChunk(iter);
        TIValue* istart = ichunk.begin;

        for (; ichunk.begin != ichunk.end; ++ichunk.begin)
        {
            if (stopFunctor(*ichunk.begin))
            {
                iter += ichunk.begin - istart;      // advance input iterator
                return;
            }
        }

        iter += ichunk.begin - istart;    // advance input iterator
    }
}

// skipUntil()
template <typename TFwdIterator, typename TStopFunctor>
inline void
skipUntil(TFwdIterator &iter, TStopFunctor &stopFunctor)
{
    _skipUntil(iter, stopFunctor, typename Chunk<TFwdIterator>::Type());
}

template <typename TFwdIterator, typename TStopFunctor>
inline void
skipUntil(TFwdIterator &iter, TStopFunctor const &stopFunctor)
{
    TStopFunctor stopFunctor_ = stopFunctor;
    _skipUntil(iter, stopFunctor_, typename Chunk<TFwdIterator>::Type());
}

// Shortcuts

template <typename TTarget, typename TFwdIterator, typename TStopFunctor>
inline void
readUntil(TTarget &target, TFwdIterator &iter, TStopFunctor const &stopFunctor)
{
    readUntil(target, iter, stopFunctor, False());
}

template <typename TTarget, typename TFwdIterator, typename TStopFunctor>
inline void
readUntil(TTarget &target, TFwdIterator &iter, TStopFunctor &stopFunctor)
{
    readUntil(target, iter, stopFunctor, False());
}

template <typename TTarget, typename TFwdIterator>
inline void
readLine(TTarget &target, TFwdIterator &iter)
{
    readUntil(target, iter, OrFunctor<EqualsChar<'\n'>, EqualsChar<'\r'> >());

    // consume "\r\n.", "\r[!\n]" or "\n."

    if (atEnd(iter))
        return;

    // If the current character is Line Feed ('\r') then this can be an ANSI or a Mac line ending.
    if (*iter == '\r')
    {
        ++iter;     // consume the found newline
        if (atEnd(iter))
            return;
    }

    // Unix Carriage Return ('\n') is the simplest case.
    if (*iter == '\n')
        ++iter;     // consume the found newline
}

template <typename TFwdIterator>
inline void
skipLine(TFwdIterator &iter)
{
    skipUntil(iter, EqualsChar<'\n'>());

    // consume "\r\n.", "\r[!\n]" or "\n."

    if (atEnd(iter))
        return;

    // If the current character is Line Feed ('\r') then this can be an ANSI or a Mac line ending.
    if (*iter == '\r')
    {
        ++iter;     // consume the found newline
        if (atEnd(iter))
            return;
    }

    // Unix Carriage Return ('\n') is the simplest case.
    if (*iter == '\n')
        ++iter;     // consume the found newline
}

template <typename TTarget, typename TFwdIterator, typename TSize>
inline TSize
readN(TTarget &target, TFwdIterator &iter, TSize n)
{
    TSize i;
    for (i = 0; !atEnd(iter) && i < n; ++i, ++iter)
        writeValue(target, value(iter));
    return i;
}





template <typename TTarget, typename TFwdIterator, typename TSize>
inline TSize
_writeN(TTarget &target, TFwdIterator &iter, TSize n)
{
    for (; n > (TSize)0; --n)
        writeValue(target, value(iter));
}

template <typename TTarget, typename TFwdIterator, typename TSize, typename TIValue, typename TOValue>
inline void
_writeN(
    TTarget &target,
    TFwdIterator &iter,
    TSize n,
    Buffer<TIValue, Simple> *,
    Buffer<TOValue, Simple> *)
{
    Buffer<TIValue, Simple> ichunk;
    Buffer<TOValue, Simple> ochunk;

    TSize minChunkSize;
    for (; n > (TSize)0; n -= minChunkSize)
    {
        ichunk = getChunk(iter);
        minChunkSize = ichunk.end - ichunk.begin;

        reserveChunk(target, minChunkSize);
        ochunk = getChunk(end(target, Rooted()));

        typename Difference<TTarget>::Type olen = ochunk.endCap - ochunk.begin;

        if (minChunkSize > olen)
        {
            minChunkSize = olen;
            ichunk.end = ichunk.begin + olen;
        }

        for (; ichunk.begin != ichunk.end; ++ichunk.begin, ++ochunk.begin)
            assignValue(ochunk.begin, getValue(ichunk.begin));

        iter += minChunkSize;                      // advance input iterator
        advanceChunk(target, minChunkSize);
    }
}







template <typename TTarget, typename TFwdIterator, typename TSize>
inline TSize
writeN(TTarget &target, TFwdIterator &iter, TSize n)
{
    typedef typename Chunk<TFwdIterator>::Type* TIChunk;
    typedef typename Chunk<TTarget>::Type*      TOChunk;

    _writeN(target, iter, n, TIChunk(), TOChunk());
}

template <typename TTarget, typename TValue, typename TSize>
inline TSize
writeN(TTarget &target, TValue *ptr, TSize n)
{
    typedef Range<TValue*>                          TRange;
    typedef typename Iterator<TRange, Rooted>::Type TIterator;
    typedef typename Chunk<TIterator>::Type*        TIChunk;
    typedef typename Chunk<TTarget>::Type*          TOChunk;

    TRange range(ptr, ptr + n);
    TIterator iter = begin(range, Rooted());
    _writeN(target, iter, n, TIChunk(), TOChunk());
}

template <typename TTarget, typename TContainer>
inline TSize
write3(TTarget &target, TContainer &cont)
{
    writeN(target, cont, length(cont));
}


// This reads Meta and Sequence
template <typename TIdString,
          typename TSeqString,
          typename TFwdIterator,
          typename TTag>
inline void
readRecord(TIdString &meta,
           TSeqString &seq,
           TFwdIterator &iter,
           Fasta)
{
    EqualsChar<'>'> fastaBegin;
    IgnoreOrAssertFunctor<
        IsWhitespace,
        IsInAlphabet<typename Value<TSeqString>::Type>(),
        std::exception > ignoreWhiteSpaceAndAssertAlphabet("Invalid character in Fasta sequence!");

    clear(seq);

    skipUntil(iter, fastaBegin);    // forward to the next '>'
    ++iter;                         // skip '>'
    readLine(meta, iter);           // read Fasta id
    readUntil(seq, iter, fastaBegin, ignoreWhiteSpaceAndAssertAlphabet);    // read Fasta sequence
}

template <typename TIdString,
          typename TSeqString,
          typename TQualString,
          typename TFwdIterator,
          typename TTag>
inline void
readRecord(TIdString &meta,
           TSeqString &seq,
           TQualString &qual,
           TFwdIterator &iter,
           Fastq)
{
    EqualsChar<'@'> fastqBegin;
    EqualsChar<'+'> qualsBegin;
    IgnoreOrAssertFunctor<
        IsWhitespace,
        IsInAlphabet<typename Value<TSeqString>::Type>(),
        std::exception > ignoreWhiteSpaceAndAssertAlphabet("Invalid character in Fastq sequence!");

    clear(seq);
    clear(qual);

    skipUntil(iter, fastqBegin);    // forward to the next '@'
    ++iter;                         // skip '@'
    readLine(meta, iter);           // read Fastq id

    // TODO(weese): Actually, we have to search for the sequence "\n+" and then for "\n@".

    readUntil(seq, iter, qualsBegin, ignoreWhiteSpaceAndAssertAlphabet);    // read Fastq sequence
    ++iter;                         // skip '+'
    skipLine(iter);                 // skip 2nd Fastq id
    readUntil(qual, iter, fastqBegin, ignoreWhiteSpaceAndAssertAlphabet);    // read Fastq qualities
}



}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_TOKENIZATION_H_
