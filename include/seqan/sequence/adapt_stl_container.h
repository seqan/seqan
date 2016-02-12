// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// Adaptions for STL containers to SeqAn sequences.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_ADAPT_STL_CONTAINER_H_
#define SEQAN_SEQUENCE_ADAPT_STL_CONTAINER_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ===========================================================================
// Concepts
// ===========================================================================

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::vector<TChar, TAlloc>), (StlContainerConcept));

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::vector<TChar, TAlloc> const), (StlContainerConcept));

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::deque<TChar, TAlloc>), (StlContainerConcept));

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::deque<TChar, TAlloc> const), (StlContainerConcept));

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::list<TChar, TAlloc>), (StlContainerConcept));

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::list<TChar, TAlloc> const), (StlContainerConcept));

template<typename TChar, typename TTraits, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::basic_string<TChar, TTraits, TAlloc>), (StlContainerConcept));

template<typename TChar, typename TTraits, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::basic_string<TChar, TTraits, TAlloc> const), (StlContainerConcept));

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::forward_list<TChar, TAlloc>), (StlContainerConcept));

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::forward_list<TChar, TAlloc> const), (StlContainerConcept));

template <typename TChar, std::size_t N>
SEQAN_CONCEPT_IMPL((std::array<TChar, N>), (StlContainerConcept));

template <typename TChar, std::size_t N>
SEQAN_CONCEPT_IMPL((std::array<TChar, N> const), (StlContainerConcept));

/* NOTE(h-2) on ConceptChecking and universal references
 *
 * When type deduction takes place before && then && can act as both && or &
 * [https://isocpp.org/blog/2012/11/universal-references-in-c11-scott-meyers]
 *
 * This has the effect that in
 *      foobar(T &&)
 * T actually can be "int &". Consequently the above becomes "int & &&" which
 * collapses to "int &". The important thing is that in
 *      template<typename EnableIf<Is<SuperConcept<T> >
 * this will result in the ConceptCheck always failing (because the Concept is
 * not defined for reference types).
 *
 * So in the following you will see a RemoveReference< > around the type
 * whenever the function takes a && argument AND it shall act as univeral
 * reference (e.g. length()). But in some cases it shall not (because there is a
 * different overload for the l-value-reference arg) and it is  intentionally
 * omitted, e.g. value().

 */

// ===========================================================================
// Metafunctions
// ===========================================================================

// ----------------------------------------------------------------------------
// Metafunctions that are the same for all STL containers
// ----------------------------------------------------------------------------

#define COMMA ,

#define SUPERMACRO__(CONT, TMPL) \
template <TMPL> \
struct Value<CONT> \
{ \
    typedef typename CONT::value_type Type; \
};\
template <TMPL> \
struct Value<CONT const> \
{ \
    typedef typename CONT::value_type Type; \
};\
template <TMPL> \
struct Reference<CONT> \
{ \
    typedef typename CONT::reference Type; \
};\
template <TMPL> \
struct Reference<CONT const> \
{ \
    typedef typename CONT::const_reference Type; \
};\
template <TMPL> \
struct GetValue<CONT> \
{ \
    typedef typename CONT::const_reference Type; \
};\
template <TMPL> \
struct GetValue<CONT const> \
{ \
    typedef typename CONT::const_reference Type; \
};\
template <TMPL> \
struct Position<CONT> \
{ \
    typedef typename CONT::size_type Type; \
};\
template <TMPL> \
struct Position<CONT const> \
{ \
    typedef typename CONT::size_type Type; \
};\
template <TMPL> \
struct Size<CONT> \
{ \
    typedef typename CONT::size_type Type; \
};\
template <TMPL> \
struct Size<CONT const> \
{ \
    typedef typename CONT::size_type Type; \
};\
template <TMPL> \
struct Iterator<CONT, Standard> \
{ \
    typedef Iter<CONT, StdIteratorAdaptor> Type;\
};\
template <TMPL> \
struct Iterator<CONT const, Standard> \
{ \
    typedef Iter<CONT const, StdIteratorAdaptor> Type;\
};\
template <TMPL> \
struct Iterator<CONT, Rooted> \
{ \
    typedef Iter<CONT, AdaptorIterator<Iter<CONT, StdIteratorAdaptor> > > Type;\
}; \
template <TMPL> \
struct Iterator<CONT const, Rooted> \
{ \
    typedef Iter<CONT const, AdaptorIterator<Iter<CONT const, StdIteratorAdaptor> > > Type;\
};

SUPERMACRO__(std::vector<TChar COMMA TAlloc>,      typename TChar COMMA typename TAlloc)
SUPERMACRO__(std::deque<TChar COMMA TAlloc>,       typename TChar COMMA typename TAlloc)
SUPERMACRO__(std::list<TChar COMMA TAlloc>,        typename TChar COMMA typename TAlloc)
SUPERMACRO__(std::forward_list<TChar COMMA TAlloc>,typename TChar COMMA typename TAlloc)
SUPERMACRO__(std::array<TChar COMMA N>,            typename TChar COMMA std::size_t N)
SUPERMACRO__(std::basic_string<TChar COMMA TTraits COMMA TAlloc>, typename TChar COMMA typename TTraits COMMA typename TAlloc)

// ----------------------------------------------------------------------------
// Mfn IsContiguous (default is False)
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
struct IsContiguous<std::vector<TChar, TAlloc> > :
    public True
{};

template <typename TChar, typename TAlloc>
struct IsContiguous<std::vector<TChar, TAlloc> const> :
    public True
{};

template <typename TChar, typename TTraits, typename TAlloc>
struct IsContiguous<std::basic_string<TChar, TTraits, TAlloc> > :
    public True
{};

template <typename TChar, typename TTraits, typename TAlloc>
struct IsContiguous<std::basic_string<TChar, TTraits, TAlloc> const> :
    public True
{};

template <typename TChar, size_t N>
struct IsContiguous<std::array<TChar, N> > :
    public True
{};

template <typename TChar, size_t N>
struct IsContiguous<std::array<TChar, N> const> :
    public True
{};

// ----------------------------------------------------------------------------
// Mfn HasSubscriptOperator (different for e.g. std::deque)
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
struct HasSubscriptOperator<std::deque<TChar, TAlloc> > :
    public True
{};

template <typename TChar, typename TAlloc>
struct HasSubscriptOperator<std::deque<TChar, TAlloc> const> :
    public True
{};

// ----------------------------------------------------------------------------
// Mfn IsSequence
// ----------------------------------------------------------------------------

//NOTE(h-2): shouldnt we be using Is<StringConcept> or IsContiguous<> or ...?

template <typename TChar, typename TAlloc>
struct IsSequence<std::vector<TChar, TAlloc> > : True {};

template <typename TChar, typename TAlloc>
struct IsSequence<std::deque<TChar, TAlloc> > : True {};

template <typename TChar, typename TAlloc>
struct IsSequence<std::list<TChar, TAlloc> > : True {};

template <typename TChar, typename TAlloc>
struct IsSequence<std::forward_list<TChar, TAlloc> > : True {};

//NOTE(h-2): array isn't by SeqAn definition, right?

template <typename TChar, typename TCharTraits, typename TAlloc>
struct IsSequence<std::basic_string<TChar, TCharTraits, TAlloc> > : True {};

// ----------------------------------------------------------------------------
// Mfn FixedSize_ (private)
// ----------------------------------------------------------------------------

template <typename TContainer>
struct FixedSize_ : public False {};

template <typename TChar, size_t N>
struct FixedSize_<std::array<TChar, N> > : public True {};

// ===========================================================================
// Functions
// ===========================================================================

// ----------------------------------------------------------------------------
// Function getObjectId
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, void const *)
getObjectId(TContainer && me)
{
    if (me.empty())
        return NULL;
    else
        return (& *(me.end() - 1)) + 1;
}

// ----------------------------------------------------------------------------
// Function begin (standard)
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer, Standard>::Type)
begin(TContainer & me, Standard const &)
{
    return me.begin();
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer const, Standard>::Type)
begin(TContainer const & me, Standard const &)
{
    return me.begin();
}

// ----------------------------------------------------------------------------
// Function begin (rooted)
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer, Rooted>::Type)
begin(TContainer & me, Rooted const &)
{
    return _beginDefault(me, Rooted());
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer const, Rooted>::Type)
begin(TContainer const & me, Rooted const &)
{
    return _beginDefault(me, Rooted());
}

// ----------------------------------------------------------------------------
// Function end (standard)
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer, Standard>::Type)
end(TContainer & me, Standard const &)
{
    return me.end();
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer const, Standard>::Type)
end(TContainer const & me, Standard const &)
{
    return me.end();
}

// ----------------------------------------------------------------------------
// Function end (rooted)
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer, Rooted>::Type)
end(TContainer & me, Rooted const &)
{
    return _endDefault(me, Rooted());
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer const, Rooted>::Type)
end(TContainer const & me, Rooted const &)
{
    return _endDefault(me, Rooted());
}

// ----------------------------------------------------------------------------
// Function _iterStl (like seqan's iter() but for stl-iterators)
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TPos>
inline typename TContainer::iterator
_iterStl(TContainer & me,
         TPos const pos)
{
    SEQAN_ASSERT_LEQ_MSG(pos,
                         static_cast<TPos>(length(me)),
                         "Trying to get an iterator behind a container through _iterStl().");
    return std::next(me.begin(), pos);
}

template <typename TContainer,
          typename TPos>
inline typename TContainer::const_iterator
_iterStl(TContainer const & me,
         TPos const pos)
{
    SEQAN_ASSERT_LEQ_MSG(pos,
                         static_cast<TPos>(length(me)),
                         "Trying to get an iterator behind a container through _iterStl().");
    return std::next(me.begin(), pos);
}

// ----------------------------------------------------------------------------
// Function iter (Standard) (only overloaded for non-contiguous)
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::list<TChar, TAlloc>, Standard>::Type
iter(std::list<TChar, TAlloc> & me,
     TPos const pos,
     Standard const &)
{
    return _iterStl(me, pos);
}

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::list<TChar, TAlloc> const, Standard>::Type
iter(std::list<TChar, TAlloc> const & me,
     TPos const pos,
     Standard const &)
{
    return _iterStl(me, pos);
}

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::forward_list<TChar, TAlloc>, Standard const>::Type
iter(std::forward_list<TChar, TAlloc> & me,
     TPos const pos,
     Standard const &)
{
    return _iterStl(me, pos);
}

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::forward_list<TChar, TAlloc> const, Standard>::Type
iter(std::forward_list<TChar, TAlloc> const & me,
     TPos const pos,
     Standard const &)
{
    return _iterStl(me, pos);
}

// ----------------------------------------------------------------------------
// Function iter (Rooted) (only overloaded for non-contiguous)
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::list<TChar, TAlloc>, Rooted>::Type
iter(std::list<TChar, TAlloc> & me,
     TPos const pos,
     Rooted const &)
{
    typedef typename Iterator<std::list<TChar, TAlloc>, Rooted>::Type TIterator;
    return TIterator(me, iter(me, pos, Standard()));
}

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::list<TChar, TAlloc> const, Rooted>::Type
iter(std::list<TChar, TAlloc> const & me,
     TPos const pos,
     Rooted const &)
{
    typedef typename Iterator<std::list<TChar, TAlloc> const, Rooted>::Type TIterator;
    return TIterator(me, iter(me, pos, Standard()));
}

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::forward_list<TChar, TAlloc>, Rooted>::Type
iter(std::forward_list<TChar, TAlloc> & me,
     TPos const pos,
     Rooted const &)
{
    typedef typename Iterator<std::forward_list<TChar, TAlloc>, Rooted>::Type TIterator;
    return TIterator(me, iter(me, pos, Standard()));
}

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::forward_list<TChar, TAlloc> const, Rooted>::Type
iter(std::forward_list<TChar, TAlloc> const & me,
     TPos const pos,
     Rooted const &)
{
    typedef typename Iterator<std::forward_list<TChar, TAlloc> const, Rooted>::Type TIterator;
    return TIterator(me, iter(me, pos, Standard()));
}

// ----------------------------------------------------------------------------
// Function length
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Size<TContainer>::Type)
length(TContainer const & me)
{
    return me.size();
}

template <typename TChar, typename TAlloc>
inline typename Size<std::forward_list<TChar, TAlloc> >::Type
length(std::forward_list<TChar, TAlloc> const & me)
{
    typename Size<std::forward_list<TChar, TAlloc> >::Type l = 0;

    for (auto it = me.begin(), itEnd = me.end(); it != itEnd; ++it, ++l);

    return l;
}

template <typename TChar, std::size_t N>
constexpr std::size_t
length(std::array<TChar, N> const & me)
{
    return me.size();
}

// ----------------------------------------------------------------------------
// Function capacity
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
inline typename Size<std::vector<TChar, TAlloc> >::Type
capacity(std::vector<TChar, TAlloc> const & me)
{
    return me.capacity();
}

template <typename TChar, typename TTraits, typename TAlloc>
inline typename Size<std::basic_string<TChar, TTraits, TAlloc> >::Type
capacity(std::basic_string<TChar, TTraits, TAlloc> const & me)
{
    return me.capacity();
}

template <typename TChar, std::size_t N>
constexpr std::size_t
capacity(std::array<TChar, N> const & me)
{
    return me.size();
}

// for other types the default overload to length holds

// ----------------------------------------------------------------------------
// Function empty
// ----------------------------------------------------------------------------

// VC2015 implements some C++17 functions which would collide for
// applications that do using namespace std
// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n4280.pdf
#ifndef _MSC_VER
template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, bool)
empty(TContainer const & me)
{
    return me.empty();
}

template <typename TChar, std::size_t N>
constexpr bool
empty(std::array<TChar, N> const & me)
{
    return me.empty();
}
#endif

// ----------------------------------------------------------------------------
// Function reserve
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc, typename TSize, typename TExpand>
inline typename Size<std::vector<TChar, TAlloc> >::Type
reserve(std::vector<TChar, TAlloc> & me,
        TSize const s,
        Tag<TExpand> const &)
{
    me.reserve(s);
    return capacity(me);
}

template <typename TChar, typename TTraits, typename TAlloc, typename TSize, typename TExpand>
inline typename Size<std::basic_string<TChar, TTraits, TAlloc> >::Type
reserve(std::basic_string<TChar, TTraits, TAlloc> & me,
        TSize const s,
        Tag<TExpand> const &)
{
    me.reserve(s);
    return capacity(me);
}

template <typename TChar, typename TAlloc, typename TSize, typename TExpand>
inline typename Size<std::vector<TChar, TAlloc> >::Type
reserve(std::vector<TChar, TAlloc> && me,
        TSize const s,
        Tag<TExpand> const &)
{
    me.reserve(s);
    return capacity(me);
}

template <typename TChar, typename TTraits, typename TAlloc, typename TSize, typename TExpand>
inline typename Size<std::basic_string<TChar, TTraits, TAlloc> >::Type
reserve(std::basic_string<TChar, TTraits, TAlloc> && me,
        TSize const s,
        Tag<TExpand> const &)
{
    me.reserve(s);
    return capacity(me);
}

// for other types the default overload holds

// ----------------------------------------------------------------------------
// Function resize
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSize,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >,
                            typename Size<typename RemoveReference<TContainer>::Type>::Type)
resize(TContainer && me,
       TSize const s,
       Tag<TExpand> const &)
{
    me.resize(s);
    return length(me);
}

template <typename TContainer,
          typename TSize,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >,
                            typename Size<typename RemoveReference<TContainer>::Type>::Type)
resize(TContainer && me,
       TSize const s,
       typename Value<typename RemoveReference<TContainer>::Type>::Type const & val,
       Tag<TExpand> const &)
{
    me.resize(s, val);
    return length(me);
}

// ----------------------------------------------------------------------------
// Function clear
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
clear(TContainer && me)
{
    me.clear();
}

// ----------------------------------------------------------------------------
// Function value
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                HasSubscriptOperator<TContainer> >, typename Reference<TContainer>::Type)
value(TContainer & me, TPos const pos)
{
    return me[pos];
}

template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                HasSubscriptOperator<TContainer> >, typename Reference<TContainer const>::Type)
value(TContainer const & me, TPos const pos)
{
    return me[pos];
}

template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                HasSubscriptOperator<TContainer> >, typename Value<TContainer>::Type)
value(TContainer && me, TPos const pos)
{
    return me[pos];
}

// linear complexity for list and fwd list
template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                Not<HasSubscriptOperator<TContainer> > >, typename Reference<TContainer>::Type)
value(TContainer & me, TPos const pos)
{
    return *(_iterStl(me, pos));
}

template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                Not<HasSubscriptOperator<TContainer> > >, typename Reference<TContainer const>::Type)
value(TContainer const & me, TPos const pos)
{
    return *(_iterStl(me, pos));
}

template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                Not<HasSubscriptOperator<TContainer> > >, typename Value<TContainer>::Type)
value(TContainer && me, TPos const pos)
{
    return *(_iterStl(me, pos));
}

// --------------------------------------------------------------------------
// Function getValue
// --------------------------------------------------------------------------

template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                HasSubscriptOperator<TContainer> >, typename Reference<TContainer const>::Type)
getValue(TContainer const & me, TPos const pos)
{
    return me[pos];
}

template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                HasSubscriptOperator<TContainer> >, typename Value<TContainer>::Type)
getValue(TContainer && me, TPos const pos)
{
    return me[pos];
}

// linear complexity for list and fwd list
template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                Not<HasSubscriptOperator<TContainer> > >, typename Reference<TContainer const>::Type)
getValue(TContainer const & me, TPos const pos)
{
    return *(_iterStl(me, pos));
}

template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                Not<HasSubscriptOperator<TContainer> > >, typename Value<TContainer>::Type)
getValue(TContainer && me, TPos const pos)
{
    return *(_iterStl(me, pos));
}

// ----------------------------------------------------------------------------
// Function front
// ----------------------------------------------------------------------------

// default implementation used

// ----------------------------------------------------------------------------
// Function back
// ----------------------------------------------------------------------------

// default implementation used

// forward_list doesnt have back, we achieve it in linear time
template <typename TChar, typename TAlloc>
inline TChar &
back(std::forward_list<TChar, TAlloc> & me)
{
    return *(std::next(me.before_begin(), length(me)));
}

template <typename TChar, typename TAlloc>
inline TChar const &
back(std::forward_list<TChar, TAlloc> const & me)
{
    return *(std::next(me.before_begin(), length(me)));
}

template <typename TChar, typename TAlloc>
inline TChar
back(std::forward_list<TChar, TAlloc> && me)
{
    return *(std::next(me.before_begin(), length(me)));
}

// ----------------------------------------------------------------------------
// Function assign
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, void)
assign(TContainer && me,
       typename RemoveReference<TContainer>::Type source)
{
    std::swap(me, source);
}

template <typename TContainer,
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Not<IsSameType<typename RemoveReference<TContainer>::Type, TSource> >,
                                Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> > >, void)
assign(TContainer && me,
       TSource const & source)
{
    me.assign(begin(source, Standard()), end(source, Standard()));
}

template <typename TContainer,
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, void)
assign(TContainer && me,
       TSource const & source,
       typename Size<TSource>::Type limit)
{
    me.assign(begin(source, Standard()),
              iter(source, std::min(length(source), limit), Standard()));
}

// ----------------------------------------------------------------------------
// Function insert
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
insert(TContainer && me,
       typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
       TSource const & source,
       Tag<TExpand> const &)
{
    me.insert(_iterStl(std::forward<TContainer>(me), pos),
              begin(source, Standard()),
              end(source, Standard()));
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
insert(TContainer && me,
       typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
       TSource const & source,
       typename Size<TSource>::Type const limit,
       Tag<TExpand> const &)
{
    me.insert(_iterStl(std::forward<TContainer>(me), pos),
              begin(source, Standard()),
              iter(source, std::min(length(source), limit), Standard()));
}

// forward_list doesnt have insert, we achieve it slower
template <typename TChar, typename TAlloc, typename TSource, typename TExpand>
inline void
insert(std::forward_list<TChar, TAlloc> & me,
       typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
       TSource const & source,
       Tag<TExpand> const &)
{
    me.insert_after(std::next(me.before_begin(), pos), begin(source, Standard()), end(source, Standard()));
}

template <typename TChar, typename TAlloc, typename TSource, typename TExpand>
inline void
insert(std::forward_list<TChar, TAlloc> && me,
       typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
       TSource const & source,
       Tag<TExpand> const &)
{
    insert(me, pos, source);
}

template <typename TChar, typename TAlloc, typename TSource, typename TExpand>
inline void
insert(std::forward_list<TChar, TAlloc> & me,
       typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
       TSource const & source,
       typename Size<std::forward_list<TChar, TAlloc> >::Type const limit,
       Tag<TExpand> const &)
{
    me.insert_after(std::next(me.before_begin(), pos),
                    begin(source, Standard()),
                    iter(source, std::min(length(source), limit), Standard()));
}

template <typename TChar, typename TAlloc, typename TSource, typename TExpand>
inline void
insert(std::forward_list<TChar, TAlloc> && me,
       typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
       TSource const & source,
       typename Size<std::forward_list<TChar, TAlloc> >::Type const limit,
       Tag<TExpand> const &)
{
    insert(me, pos, source, limit);
}

// ----------------------------------------------------------------------------
// Function append
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
append(TContainer && me,
       TSource const & source,
       Tag<TExpand> const &)
{
    insert(std::forward<TContainer>(me), length(me), source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
append(TContainer && me,
       TSource const & source,
       typename Size<typename RemoveReference<TContainer>::Type>::Type const limit,
       Tag<TExpand> const &)
{
    insert(std::forward<TContainer>(me), length(me), source, limit);
}

// forward list is handled downstream in insert

// ----------------------------------------------------------------------------
// Function prepend
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
prepend(TContainer && me,
        TSource const & source)
{
    insert(std::forward<TContainer>(me), 0, source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
prepend(TContainer && me,
        TSource const & source,
        Tag<TExpand> const &)
{
    insert(std::forward<TContainer>(me), 0, source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
prepend(TContainer && me,
        TSource const & source,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const limit)
{
    insert(std::forward<TContainer>(me), 0, source, limit);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
prepend(TContainer && me,
        TSource const & source,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const limit,
        Tag<TExpand> const &)
{
    insert(std::forward<TContainer>(me), 0, source, limit);
}

// forward list is handled downstream in insert

// ----------------------------------------------------------------------------
// Function insertValue
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
insertValue(TContainer && me,
            typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
            TSource && source,
            Tag<TExpand> const &)
{
    me.insert(_iterStl(std::forward<TContainer>(me), pos), std::forward<TSource>(source));
}

template <typename TChar, typename TAlloc, typename TSource, typename TExpand>
inline void
insertValue(std::forward_list<TChar, TAlloc> & me,
            typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
            TSource && source,
            Tag<TExpand> const &)
{
    me.insert_after(std::next(me.before_begin(), pos), std::forward<TSource>(source));
}

template <typename TChar, typename TAlloc, typename TSource, typename TExpand>
inline void
insertValue(std::forward_list<TChar, TAlloc> && me,
            typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
            TSource && source,
            Tag<TExpand> const &)
{
    insertValue(me, pos, source);
}

// ----------------------------------------------------------------------------
// Function appendValue
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
appendValue(TContainer && me,
            TSource && source,
            Tag<TExpand> const &)
{
    me.push_back(std::forward<TSource>(source));
}

template <typename TChar, typename TAlloc, typename TSource, typename TExpand>
inline void
appendValue(std::forward_list<TChar, TAlloc> & me,
            TSource && source,
            Tag<TExpand> const &)
{
    insertValue(me, length(me), std::forward<TSource>(source));
}

template <typename TChar, typename TAlloc, typename TSource, typename TExpand>
inline void
appendValue(std::forward_list<TChar, TAlloc> && me,
            TSource && source,
            Tag<TExpand> const &)
{
    appendValue(me, source);
}

// ----------------------------------------------------------------------------
// Function prependValue
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
prependValue(TContainer && me,
             TSource && source)
{
    me.push_front(std::forward<TSource>(source));
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
prependValue(TContainer && me,
             TSource && source,
             Tag<TExpand> const &)
{
    prependValue(std::forward<TContainer>(me), std::forward<TSource>(source));
}

// vector and string dont have push_front, we achieve it in linear time
template <typename TChar, typename TAlloc, typename TSource>
inline void
prependValue(std::vector<TChar, TAlloc> & me,
             TSource && source)
{
    insertValue(me, 0, std::forward<TSource>(source));
}

template <typename TChar, typename TTraits, typename TAlloc, typename TSource>
inline void
prependValue(std::basic_string<TChar, TTraits, TAlloc> & me,
             TSource && source)
{
    insertValue(me, 0, std::forward<TSource>(source));
}

// && to &
template <typename TChar, typename TAlloc, typename TSource>
inline void
prependValue(std::vector<TChar, TAlloc> && me,
             TSource && source)
{
    prependValue(me, std::forward<TSource>(source));
}

template <typename TChar, typename TTraits, typename TAlloc, typename TSource>
inline void
prependValue(std::basic_string<TChar, TTraits, TAlloc> && me,
             TSource && source)
{
    prependValue(me, std::forward<TSource>(source));
}

// ----------------------------------------------------------------------------
// Function erase
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
erase(TContainer && me,
      typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
      typename Size<typename RemoveReference<TContainer>::Type>::Type const posEnd)
{
    me.erase(_iterStl(me, pos), _iterStl(me, posEnd));
}

template <typename TChar, typename TAlloc>
inline void
erase(std::forward_list<TChar, TAlloc> & me,
      typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
      typename Size<std::forward_list<TChar, TAlloc> >::Type const posEnd)
{
    me.erase_after(std::next(me.before_begin(), pos), std::next(me.begin(), posEnd));
}

template <typename TChar, typename TAlloc>
inline void
erase(std::forward_list<TChar, TAlloc> && me,
      typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
      typename Size<std::forward_list<TChar, TAlloc> >::Type const posEnd)
{
    erase(me, pos, posEnd);
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
erase(TContainer && me,
      typename Size<typename RemoveReference<TContainer>::Type>::Type const pos)
{
    erase(std::forward<TContainer>(me), pos, pos + 1);
}

// ----------------------------------------------------------------------------
// Function eraseFront
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
eraseFront(TContainer && me)
{
    me.pop_front();
}

template <typename TChar, typename TAlloc>
inline void
eraseFront(std::vector<TChar, TAlloc> & me)
{
    me.erase(me.begin());
}

template <typename TChar, typename TTraits, typename TAlloc>
inline void
eraseFront(std::basic_string<TChar, TTraits, TAlloc> & me)
{
    me.erase(me.begin());
}

template <typename TChar, typename TAlloc>
inline void
eraseFront(std::vector<TChar, TAlloc> && me)
{
    me.erase(me.begin());
}

template <typename TChar, typename TTraits, typename TAlloc>
inline void
eraseFront(std::basic_string<TChar, TTraits, TAlloc> && me)
{
    me.erase(me.begin());
}

// ----------------------------------------------------------------------------
// Function eraseBack
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
eraseBack(TContainer && me)
{
    me.pop_back();
}

template <typename TChar, typename TAlloc>
inline void
eraseBack(std::forward_list<TChar, TAlloc> & me)
{
    if (!empty(me))
    {
        auto it = me.before_begin(), itN = std::next(me.begin()), itEnd = me.end();
        while (itN != itEnd)
        {
            ++it;
            ++itN;
        }
        me.erase_after(it);
    }
}

template <typename TChar, typename TAlloc>
inline void
eraseBack(std::forward_list<TChar, TAlloc> && me)
{
    eraseBack(me);
}

// ----------------------------------------------------------------------------
// Function replace
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
replace(TContainer && me,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const posEnd,
        TSource const & source,
        typename Size<TSource>::Type const limit,
        Tag<TExpand> const &)
{
    erase(me, pos, posEnd);
    insert(me, pos, source, limit);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
replace(TContainer && me,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const posEnd,
        TSource const & source,
        Tag<TExpand> const &)
{
    erase(me, pos, posEnd);
    insert(me, pos, source, length(source));
}

// ----------------------------------------------------------------------------
// Function reverse
// ----------------------------------------------------------------------------

// NOTE(h-2): seqan's reverse doesn't work for some stl containers
template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >, HasSubscriptOperator<TContainer> >)
reverse(TContainer && me)
{
    std::reverse(me.begin(), me.end());
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >, Not<HasSubscriptOperator<TContainer> > >)
reverse(TContainer && me)
{
    me.reverse();
}

// ----------------------------------------------------------------------------
// Function toCString
// ----------------------------------------------------------------------------

template <typename TChar, typename TTraits, typename TAlloc>
inline const TChar*
toCString(std::basic_string<TChar, TTraits, TAlloc> const & me)
{
    return me.c_str();
}

// ----------------------------------------------------------------------------
// Function assign (from CString)
// ----------------------------------------------------------------------------

// NOTE(h-2): is this optimal if container is not Contiguous?
template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, void)
assign(TContainer & me,
       char const * source)
{
    resize(me, length(source));
    std::copy(source, source + std::strlen(source), me.begin());
}

template <typename TTraits, typename TAlloc>
inline void
assign(std::basic_string<char, TTraits, TAlloc> & me,
       char const * source)
{
    me = source;
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, void)
assign(TContainer && me,
       char const * source)
{
    assign(me, source);
}

template <typename TTraits, typename TAlloc>
inline void
assign(std::basic_string<char, TTraits, TAlloc> && me,
       char const * source)
{
    me = source;
}

}

#undef SUPERMACRO__
#undef COMMA

#endif
