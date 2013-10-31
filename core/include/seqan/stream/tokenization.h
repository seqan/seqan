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
// Tokenization.
// ==========================================================================

#ifndef SEQAN_STREAM_TOKENIZATION_H_
#define SEQAN_STREAM_TOKENIZATION_H_

namespace seqan {

// ============================================================================
// Exceptions
// ============================================================================

struct ParseError : RuntimeError
{
    template <typename TString>
    ParseError(TString const &message):
        RuntimeError(message)
    {}
};

struct UnexpectedEnd : ParseError
{
    UnexpectedEnd():
        ParseError("Unexpected end of file.")
    {}
};

// ============================================================================
// Functors
// ============================================================================
// TODO(esiragusa): move these functors into basic

template <typename TFunctor, typename TContext = void>
struct FunctorErrorMessage
{
    static const std::string VALUE;
};

// ----------------------------------------------------------------------------
// Functor OrFunctor
// ----------------------------------------------------------------------------

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
    bool operator() (TValue const & val) const
    {
        return func1(val) || func2(val);
    }
};

// ----------------------------------------------------------------------------
// Functor AndFunctor
// ----------------------------------------------------------------------------

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
    bool operator() (TValue const & val) const
    {
        return func1(val) && func2(val);
    }
};

// ----------------------------------------------------------------------------
// Functor NotFunctor
// ----------------------------------------------------------------------------

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
    bool operator() (TValue const & val) const
    {
        return !func(val);
    }
};

// ----------------------------------------------------------------------------
// Functor IsInAlphabet
// ----------------------------------------------------------------------------

template <typename TValue>
struct IsInAlphabet
{
    template <typename TInValue>
    bool operator() (TInValue const & inVal) const
    {
        TValue val = inVal;
        return convert<TInValue>(val) == toUpperValue(inVal);
    }

    bool operator() (TValue const &) const
    {
        return true;
    }
};

// ----------------------------------------------------------------------------
// Functor IsInRange
// ----------------------------------------------------------------------------

template <char FIRST_CHAR, char LAST_CHAR>
struct IsInRange
{
    template <typename TValue>
    bool operator() (TValue const & val) const
    {
        return FIRST_CHAR <= val && val <= LAST_CHAR;
    }
};

template <char FIRST_CHAR, char LAST_CHAR, typename TContext>
struct FunctorErrorMessage<IsInRange<FIRST_CHAR, LAST_CHAR>, TContext>
{
    static const std::string VALUE;
};

template <char FIRST_CHAR, char LAST_CHAR, typename TContext>
const std::string FunctorErrorMessage<IsInRange<FIRST_CHAR, LAST_CHAR>, TContext>::VALUE = std::string("Character in range'") + FIRST_CHAR + "' to '" + LAST_CHAR + "' expected.";

// ----------------------------------------------------------------------------
// Functor EqualsChar
// ----------------------------------------------------------------------------

template <char VALUE>
struct EqualsChar
{
    template <typename TValue>
    bool operator() (TValue const & val) const
    {
        return val == VALUE;
    }
};

template <char CHAR, typename TContext>
struct FunctorErrorMessage<EqualsChar<CHAR>, TContext>
{
    static const std::string VALUE;
};

template <char CHAR, typename TContext>
const std::string FunctorErrorMessage<EqualsChar<CHAR>, TContext>::VALUE = std::string("Character '") + CHAR + "' expected.";

// ----------------------------------------------------------------------------
// Functor AssertFunctor
// ----------------------------------------------------------------------------

template <typename TFunctor, typename TException, bool RETURN_VALUE = false>
struct AssertFunctor
{
    TFunctor func;

    AssertFunctor(TFunctor &func):
        func(func)
    {}

    template <typename TValue>
    bool operator() (TValue const & val) const
    {
        if (SEQAN_UNLIKELY(!func(val)))
            throw TException(std::string("Character '") + val + "' causes an error. " + FunctorErrorMessage<TFunctor>::VALUE);
        return RETURN_VALUE;
    }
};

// ============================================================================
// Functors
// ============================================================================
// Don't use isblank() or isspace() as it they seem to be slower than our functors (due to inlining)

typedef OrFunctor<EqualsChar<' '>, EqualsChar<'\t'> >           IsBlank;
typedef OrFunctor<EqualsChar<'\n'>, EqualsChar<'\r'> >          IsNewline;
typedef OrFunctor<IsBlank, IsNewline>                           IsWhitespace;
typedef IsInRange<'!', '~'>                                     IsGraph;
typedef OrFunctor<IsInRange<'a', 'z'>, IsInRange<'A', 'Z'> >    IsAlpha;
typedef IsInRange<'0', '9'>                                     IsDigit;
typedef OrFunctor<IsAlpha, IsDigit>                             IsAlphaNum;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _skipUntil(); Element-wise
// ----------------------------------------------------------------------------

template <typename TFwdIterator, typename TStopFunctor, typename TChunk>
inline void _skipUntil(TFwdIterator &iter, TStopFunctor &stopFunctor, TChunk)
{
    for (; !atEnd(iter) && !stopFunctor(*iter); ++iter) ;
}

// ----------------------------------------------------------------------------
// Function _skipUntil(); Chunked
// ----------------------------------------------------------------------------

template <typename TFwdIterator, typename TStopFunctor, typename TValue>
inline void _skipUntil(TFwdIterator &iter, TStopFunctor &stopFunctor, Range<TValue*> *)
{
    typedef typename Value<TFwdIterator>::Type TIValue;

    for (; !atEnd(iter); )
    {
        Range<TIValue*> const ichunk = getChunk(iter, Input());
        SEQAN_ASSERT(begin(ichunk, Standard()) < end(ichunk, Standard()));

        register const TIValue* __restrict__ ptr = begin(ichunk, Standard());

        for (; ptr != end(ichunk, Standard()); ++ptr)
        {
            if (SEQAN_UNLIKELY(stopFunctor(*ptr)))
            {
                iter += ptr - begin(ichunk, Standard());    // advance input iterator
                return;
            }
        }

        iter += ptr - begin(ichunk, Standard());            // advance input iterator
    }
}

// ----------------------------------------------------------------------------
// Function skipUntil()
// ----------------------------------------------------------------------------

template <typename TFwdIterator, typename TStopFunctor>
inline void skipUntil(TFwdIterator &iter, TStopFunctor &stopFunctor)
{
    typedef typename Chunk<TFwdIterator>::Type* TIChunk;

    _skipUntil(iter, stopFunctor, TIChunk());
}

template <typename TFwdIterator, typename TStopFunctor>
inline void skipUntil(TFwdIterator &iter, TStopFunctor const &stopFunctor)
{
    typedef typename Chunk<TFwdIterator>::Type* TIChunk;

    TStopFunctor stopFunctor_ = stopFunctor;
    _skipUntil(iter, stopFunctor_, TIChunk());
}

// ----------------------------------------------------------------------------
// Function skipOne()
// ----------------------------------------------------------------------------

template <typename TException, typename TFwdIterator, typename TFunctor>
inline void skipOne(TFwdIterator &iter, TFunctor &functor)
{
    AssertFunctor<TFunctor, TException> assertFunctor(functor);

    if (atEnd(iter))
        throw UnexpectedEnd();

    assertFunctor(*iter);
    ++iter;
}

template <typename TException, typename TFwdIterator>
inline void skipOne(TFwdIterator &iter)
{
    True func;
    skipOne(iter, func);
}

// ----------------------------------------------------------------------------
// Function _readUntil(); Element-wise
// ----------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------
// Function _readUntil(); Chunked
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TStopFunctor, typename TIgnoreFunctor, typename TIValue, typename TOValue>
inline void _readUntil(TTarget &target,
                       TFwdIterator &iter,
                       TStopFunctor &stopFunctor,
                       TIgnoreFunctor &ignoreFunctor,
                       Range<TIValue*> *,
                       Range<TOValue*> *)
{
    for (; !atEnd(iter); )
    {
        Range<TIValue*> const ichunk = getChunk(iter, Input());
        SEQAN_ASSERT(begin(ichunk, Standard()) < end(ichunk, Standard()));

        // reserve memory for the worst-case
        // TODO(weese):Document worst-case behavior
        reserveChunk(target, length(ichunk));

        Range<TOValue*> const ochunk = getChunk(end(target, Rooted()), Output());
        SEQAN_ASSERT(begin(ochunk, Standard()) < end(ochunk, Standard()));

        register const TIValue* __restrict__ iptr = begin(ichunk, Standard());
        register       TOValue* __restrict__ optr = begin(ochunk, Standard());

        for (; iptr != end(ichunk, Standard()); ++iptr)
        {
            if (SEQAN_UNLIKELY(stopFunctor(*iptr)))
            {
                iter += iptr - begin(ichunk, Standard());               // advance input iterator
                advanceChunk(target, optr - begin(ochunk, Standard())); // extend target string size
                return;
            }
            if (SEQAN_LIKELY(!ignoreFunctor(*iptr)))
            {
                // construct values in reserved memory
                *optr = *iptr;
                if (++optr == end(ochunk, Standard()))
                    break;
            }
        }

        iter += iptr - begin(ichunk, Standard());                       // advance input iterator
        advanceChunk(target, optr - begin(ochunk, Standard()));
    }
}

// ----------------------------------------------------------------------------
// Function readUntil()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TStopFunctor, typename TIgnoreFunctor>
inline void readUntil(TTarget &target, TFwdIterator &iter, TStopFunctor &stopFunctor, TIgnoreFunctor &ignoreFunctor)
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

// ----------------------------------------------------------------------------
// Function readUntil(); Not ignoring
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TStopFunctor>
inline void readUntil(TTarget &target, TFwdIterator &iter, TStopFunctor const &stopFunctor)
{
    readUntil(target, iter, stopFunctor, False());
}

template <typename TTarget, typename TFwdIterator, typename TStopFunctor>
inline void readUntil(TTarget &target, TFwdIterator &iter, TStopFunctor &stopFunctor)
{
    readUntil(target, iter, stopFunctor, False());
}

// ----------------------------------------------------------------------------
// Function readLine()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator>
inline void readLine(TTarget &target, TFwdIterator &iter)
{
    readUntil(target, iter, IsNewline());

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

// ----------------------------------------------------------------------------
// Function skipLine()
// ----------------------------------------------------------------------------

template <typename TFwdIterator>
inline void skipLine(TFwdIterator &iter)
{
    skipUntil(iter, IsNewline());

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

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_TOKENIZATION_H_
