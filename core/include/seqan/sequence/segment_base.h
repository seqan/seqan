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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Declarations related to and implementation of the Segment class.
// ==========================================================================

#ifndef SEQAN_HEADER_SEGMENT_BASE_H
#define SEQAN_HEADER_SEGMENT_BASE_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Segment
//////////////////////////////////////////////////////////////////////////////

/**
.Class.Segment:
..cat:Sequences
..summary:A contiguous part of a sequence.
..signature:Segment<THost, TSpec>
..param.THost:Type of the whole sequence.
...metafunction:Metafunction.Host
...text:Instances of $Segment<THost, TSpec>$ are subsequences of $THost$ objects.
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:@Spec.InfixSegment@.
*/

struct InfixSegment {};

template <typename THost, typename TSpec = InfixSegment>
class Segment
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Host.param.T.type:Class.Segment
///.Metafunction.Host.class:Class.Segment

template <typename THost, typename TSpec>
struct Host<Segment<THost, TSpec> >
{
    typedef THost Type;
};

template <typename THost, typename TSpec>
struct Host<Segment<THost, TSpec> const >
{
    typedef THost Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.Segment
///.Metafunction.Spec.class:Class.Segment

template <typename THost, typename TSpec>
struct Spec<Segment<THost, TSpec> >
{
    typedef TSpec Type;
};
template <typename THost, typename TSpec>
struct Spec<Segment<THost, TSpec> const>
{
    typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Segment
///.Metafunction.Value.class:Class.Segment

template <typename THost, typename TSpec>
struct Value<Segment<THost, TSpec> >
{
    typedef typename Value<THost>::Type Type;
};

template <typename THost, typename TSpec>
struct Value<Segment<THost, TSpec> const >
{
    typedef typename Value<THost const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.GetValue.param.T.type:Class.Segment
///.Metafunction.GetValue.class:Class.Segment

template <typename THost, typename TSpec>
struct GetValue<Segment<THost, TSpec> >
{
    typedef typename GetValue<THost>::Type Type;
};

template <typename THost, typename TSpec>
struct GetValue<Segment<THost, TSpec> const >
{
    typedef typename GetValue<THost const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Reference.param.T.type:Class.Segment
///.Metafunction.Reference.class:Class.Segment

template <typename THost, typename TSpec>
struct Reference<Segment<THost, TSpec> >
{
	typedef typename Reference<THost>::Type Type;
};

template <typename THost, typename TSpec>
struct Reference<Segment<THost, TSpec> const >
{
	typedef typename Reference<THost>::Type Type;
};
	
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Class.Segment
///.Metafunction.Iterator.class:Class.Segment

template <typename THost, typename TSpec>
struct Iterator<Segment<THost, TSpec>, Rooted>
{
    typedef Segment<THost, TSpec> TSequence_;
    typedef typename Iterator<THost, Standard>::Type TIterator_;
    typedef Iter<TSequence_, AdaptorIterator<TIterator_> > Type;
};
template <typename THost, typename TSpec>
struct Iterator<Segment<THost, TSpec> const, Rooted>
{
    typedef Segment<THost, TSpec> const TSequence_;
    typedef typename Iterator<THost const, Standard>::Type TIterator_;
    typedef Iter<TSequence_, AdaptorIterator<TIterator_> > Type;
};

template <typename THost, typename TSpec>
struct Iterator<Segment<THost, TSpec>, Standard>:
    Iterator<THost, Standard>
{
};
template <typename THost, typename TSpec>
struct Iterator<Segment<THost, TSpec> const, Standard>:
    Iterator<THost, Standard>
{
};



//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Size.param.T.type:Class.Segment
///.Metafunction.Size.class:Class.Segment

template <typename THost, typename TSpec>
struct Size<Segment<THost, TSpec> >
{
    typedef typename Size<THost>::Type Type;
};

template <typename THost, typename TSpec>
struct Size<Segment<THost, TSpec> const >
{
    typedef typename Size<THost>::Type Type;
};

template <typename THost, typename TSpec>
struct Position<Segment<THost, TSpec> >
{
    typedef typename Position<THost>::Type Type;
};

template <typename THost, typename TSpec>
struct Position<Segment<THost, TSpec> const >
{
    typedef typename Position<THost>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.DefaultOverflowImplicit.param.T.type:Class.Segment
///.Metafunction.DefaultOverflowImplicit.class:Class.Segment

template <typename THost, typename TSpec>
struct DefaultOverflowImplicit<Segment<THost, TSpec > >:
    DefaultOverflowImplicit<THost>
{
};

template <typename THost, typename TSpec>
struct DefaultOverflowImplicit<Segment<THost, TSpec > const >:
    DefaultOverflowImplicit<THost>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.DefaultOverflowExplicit.param.T.type:Class.Segment
///.Metafunction.DefaultOverflowExplicit.class:Class.Segment

template <typename THost, typename TSpec>
struct DefaultOverflowExplicit<Segment<THost, TSpec > >:
    DefaultOverflowExplicit<THost>
{
};

template <typename THost, typename TSpec>
struct DefaultOverflowExplicit<Segment<THost, TSpec > const >:
    DefaultOverflowExplicit<THost>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IsContiguous.param.T.type:Class.Segment
///.Metafunction.IsContiguous.class:Class.Segment

template <typename THost, typename TSpec>
struct IsContiguous< Segment<THost, TSpec> >:
    public IsContiguous<THost> {};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IsSequence.param.T.type:Class.Segment
///.Metafunction.IsSequence.class:Class.Segment

template <typename THost, typename TSpec>
struct IsSequence< Segment<THost, TSpec> > {
    typedef True Type;
    enum { VALUE = true };
};

//////////////////////////////////////////////////////////////////////////////

///.Function.atBegin.param.iterator.type:Class.Segment
///.Function.atBegin.class:Class.Segment
///.Function.atEnd.param.iterator.type:Class.Segment
///.Function.atEnd.class:Class.Segment
///.Function.goBegin.param.iterator.type:Class.Segment
///.Function.goBegin.class:Class.Segment
///.Function.goEnd.param.iterator.type:Class.Segment
///.Function.goEnd.class:Class.Segment
///.Function.goNext.param.iterator.type:Class.Segment
///.Function.goNext.class:Class.Segment
///.Function.goPrevious.param.iterator.type:Class.Segment
///.Function.goPrevious.class:Class.Segment
///.Function.value.param.container.type:Class.Segment
///.Function.value.class:Class.Segment

///.Function.shareResources.param.sequence1, sequence2.type:Class.Segment
///.Function.shareResources.class:Class.Segment

//////////////////////////////////////////////////////////////////////////////
// functions for all Segment classes
//////////////////////////////////////////////////////////////////////////////

///.Function.getObjectId.param.object.type:Class.Segment
///.Function.getObjectId.class:Class.Segment

template <typename THost, typename TSpec>
inline void const *
getObjectId(Segment<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return getObjectId(host(me));
}

//////////////////////////////////////////////////////////////////////////////

///.Function.empty.param.object.type:Class.Segment

template <typename THost, typename TSpec>
inline bool
empty(Segment<THost, TSpec> const & me)
{
    SEQAN_CHECKPOINT;
    return (beginPosition(me) == endPosition(me));
}

//////////////////////////////////////////////////////////////////////////////

///.Function.length.param.object.type:Class.Segment
///.Function.length.class:Class.Segment

template <typename THost, typename TSpec>
inline typename Size<Segment<THost, TSpec> const>::Type
length(Segment<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return endPosition(me) - beginPosition(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.capacity.param.object.type:Class.Segment
///.Function.capacity.class:Class.Segment

template <typename THost, typename TSpec>
inline typename Size< Segment<THost, TSpec> const>::Type
capacity(Segment<THost, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return capacity(host(me)) + length(me) - length(host(me));
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
inline bool
hasNoHost(Segment<THost, TSpec> const & target)
{
    return !_toPointer(host(target));
}


//////////////////////////////////////////////////////////////////////////////

// operationToSet testet, ob statt einer assign/append-Funktion
// eine set Funktion aufgerufen werden soll. Das ist nur dann der Fall,
// wenn target keinen Host hat und source einen kompatiblen Host
// anbietet.
// returns true:  set Funktion verwendet,
//                wurde bereits von test_operation_2_set gemacht
// returns false: keine set Funktion verwenden
//                muss noch assign/append gemacht werden.

template <typename TSameSpec, typename TTargetInfix>
struct SegmentSetImpl_
{
    template <typename TTarget, typename TSource>
    static inline bool operationToSet(TTarget &, TSource &)
    {
        return false;
    }
};

template <typename TTargetInfix>
struct SegmentSetImpl_<True, TTargetInfix>
{
    template <typename TTarget, typename TSource>
    static inline bool operationToSet(TTarget & target, TSource & source)
    {
        set(target, source);
        return true;
    }
};

template <>
struct SegmentSetImpl_<False, True>
{
    template <typename TTarget, typename TSource>
    static inline bool operationToSet(TTarget & target, TSource & source)
    {
        set(target, host(source), beginPosition(source), endPosition(source));
        return true;
    }
};

template <typename THost, typename TSpec, typename TSource>
inline bool
operationToSet(Segment<THost, TSpec> &,
                TSource &)
{
    return false;
}
template <typename THost, typename TSpec, typename TSpec2>
inline bool
operationToSet(Segment<THost, TSpec> & target,
                Segment<THost, TSpec2> & source)
{
    if (IsConst_<THost>::VALUE || hasNoHost(target))
        return SegmentSetImpl_<
            typename IsSameType<TSpec, TSpec2>::Type,
            typename IsSameType<TSpec, InfixSegment>::Type
        >::operationToSet(target, source);
    return false;
}
template <typename THost, typename TSpec, typename TSpec2>
inline bool
operationToSet(Segment<THost const, TSpec> & target,
                Segment<THost, TSpec2> & source)
{
    // TODO(weese): make the compiler break if this function returns false
    return SegmentSetImpl_<
        typename IsSameType<TSpec, TSpec2>::Type,
        typename IsSameType<TSpec, InfixSegment>::Type
    >::operationToSet(target, source);
}
template <typename THost, typename TSpec, typename TSpec2>
inline bool
operationToSet(Segment<THost const, TSpec> & target,
                Segment<THost, TSpec2> const & source)
{
    // TODO(weese): make the compiler break if this function returns false
    return SegmentSetImpl_<
        typename IsSameType<TSpec, TSpec2>::Type,
        typename IsSameType<TSpec, InfixSegment>::Type
    >::operationToSet(target, source);
}
template <typename THost, typename TSpec, typename TSpec2>
inline bool
operationToSet(Segment<THost, TSpec> & target,
                Segment<THost, TSpec2> const & source)
{
    if (IsConst_<THost>::VALUE || hasNoHost(target))
        return SegmentSetImpl_<
            typename IsSameType<TSpec, TSpec2>::Type,
            typename IsSameType<TSpec, InfixSegment>::Type
        >::operationToSet(target, source);
    return false;
}
template <typename THost, typename TSpec>
inline bool
operationToSet(Segment<THost, TSpec> & target,
                THost & source)
{
    if (hasNoHost(target))
    {
        set(target, source);
        return true;
    }
    return false;
}


template <typename THost, typename TSpec, typename TSource, typename TSize>
inline bool
operationToSet(Segment<THost, TSpec> & /*target*/,
                TSource & /*source*/,
                TSize)
{
    return false;
}
template <typename THost, typename TSpec, typename TSpec2, typename TSize>
inline bool
operationToSet(Segment<THost, TSpec> & target,
                Segment<THost, TSpec2> & source,
                TSize limit)
{
    if (hasNoHost(target))
    {
        TSize beginpos = beginPosition(source);
        TSize endpos = endPosition(source);
        if (endpos - beginpos > limit)
        {
            endpos = beginpos + limit;
        }
        set(target, host(source), beginpos, endpos);
        return true;
    }
    return false;
}
template <typename THost, typename TSpec, typename TSize>
inline bool
operationToSet(Segment<THost, TSpec> & target,
                THost & source,
                TSize limit)
{
    if (hasNoHost(target))
    {
        TSize size = length(source);
        if (size > limit)
        {
            size = limit;
        }
        set(target, source, 0, size);
        return true;
    }
    return false;
}


//////////////////////////////////////////////////////////////////////////////
// assign
//////////////////////////////////////////////////////////////////////////////

/**
.Function.assign:
..class:Class.Segment
..remarks:If $target$ is a @Class.Segment@ object, then
$limit$ denotes the maximal length of @Function.host.$host(target)$@ after the operation.
..param.target.type:Class.Segment
..param.source.type:Class.Segment
*/

// TODO(holtgrew): We'd rather only have one version.

template<typename THost, typename TSpec>
inline void
assign(Segment<THost, TSpec> & target,
       Segment<THost, TSpec> const & source)
{
    typedef Segment<THost, TSpec> TSegment;
    assign(target, source, typename DefaultOverflowImplicit<TSegment>::Type());
}

template<typename THost, typename TSpec>
inline void
assign(Segment<THost, TSpec> & target,
       Segment<THost, TSpec> & source)
{
    typedef Segment<THost, TSpec> TSegment;
    assign(target, source, typename DefaultOverflowImplicit<TSegment>::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct AssignSegment_
{
    template <typename THost, typename TSpec, typename TSource>
    static inline void
    assign_(
        Segment<THost const, TSpec> & target,
        TSource & source)
    {
        set(target, source);
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    assign_(
        Segment<THost const, TSpec> & target,
        TSource & source,
        typename Size< Segment<THost const, TSpec> >::Type /*limit*/)
    {
        set(target, source);
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    assign_(
        Segment<THost, TSpec> & target,
        TSource & source)
    {
        set(target, source);
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    assign_(
        Segment<THost, TSpec> & target,
        TSource & source,
        typename Size< Segment<THost, TSpec> >::Type limit)
    {
        (void)limit;
        set(target, source);
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    assign_(
        Segment<THost, TSpec> const & target,
        TSource & source)
    {
        set(target, source);
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    assign_(
        Segment<THost, TSpec> const & target,
        TSource & source,
        typename Size< Segment<THost, TSpec> >::Type limit)
    {
        (void)limit;
        set(target, source);
    }
};

//____________________________________________________________________________

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> & target,
       TSource & source,
       Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> const>::assign_(target, source);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> & target,
       TSource const & source,
       Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> const>::assign_(target, source);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> & target,
       TSource & source,
       typename Size< Segment<THost, TSpec> >::Type limit,
       Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> const>::assign_(target, source, limit);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> & target,
       TSource const & source,
       typename Size< Segment<THost, TSpec> >::Type limit,
       Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> const>::assign_(target, source, limit);
}

//(for temporary targets)

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> const & target,
       TSource & source,
       Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> const>::assign_(target, source);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> const & target,
       TSource const & source,
       Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> const>::assign_(target, source);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> const & target,
       TSource & source,
       typename Size< Segment<THost, TSpec> >::Type limit,
       Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> const>::assign_(target, source, limit);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
assign(Segment<THost, TSpec> const & target,
       TSource const & source,
       typename Size< Segment<THost, TSpec> >::Type limit,
       Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AssignSegment_<Tag<TExpand> const>::assign_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
// move
//////////////////////////////////////////////////////////////////////////////

//overload of binary version:

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
move(Segment<THost, TSpec> & target,
     TSource & source)
{
SEQAN_CHECKPOINT
    typedef Segment<THost, TSpec> TTarget;
    move(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
move(Segment<THost, TSpec> & target,
     TSource const & source)
{
SEQAN_CHECKPOINT
    typedef Segment<THost, TSpec> TTarget;
    move(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

//(for temporary targets)

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
move(Segment<THost, TSpec> const & target,
     TSource & source)
{
SEQAN_CHECKPOINT
    typedef Segment<THost, TSpec> const TTarget;
    move(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
move(Segment<THost, TSpec> const & target,
     TSource const & source)
{
SEQAN_CHECKPOINT
    typedef Segment<THost, TSpec> const TTarget;
    move(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}



//////////////////////////////////////////////////////////////////////////////
// append
//////////////////////////////////////////////////////////////////////////////

/**
.Function.append:
..class:Class.Segment
..remarks:If $target$ is a @Class.Segment@ object, then
$limit$ denotes the maximal length of @Function.host.$host(target)$@ after the operation.
..param.target.type:Class.Segment
..param.source.type:Class.Segment
*/


template <typename TExpand>
struct AppendSequenceToSegment_
{
    template <typename THost, typename TSpec, typename TSource>
    static inline void
    append_(
        Segment<THost, TSpec> & target,
        TSource & source)
    {
SEQAN_CHECKPOINT
        typedef Segment<THost, TSpec> Target;

        if (!operationToSet(target, source))
        {
            replace(host(target), endPosition(target), endPosition(target), source, TExpand());

            typename Iterator<Target, Standard>::Type new_end = end(target, Standard()) + length(source);
            typename Iterator<THost, Standard>::Type host_end = end(host(target), Standard());
            if (new_end > host_end) new_end = host_end;
            setEnd(target, new_end);
        }
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    append_(
        Segment<THost, TSpec> & target,
        TSource const & source,
        typename Size< Segment<THost, TSpec> >::Type limit)
    {
SEQAN_CHECKPOINT
        typedef Segment<THost, TSpec> Target;

        if (!operationToSet(target, source))
        {
            replace(host(target), endPosition(target), endPosition(target), source, limit, TExpand());
            typename Iterator<Target, Standard>::Type new_end = end(target, Standard()) + length(source);
            typename Iterator<THost, Standard>::Type host_end = end(host(target), Standard());
            if (begin(target) > host_end) setBegin(target, host_end);
            if (new_end > host_end) new_end = host_end;
            setEnd(target, new_end);
        }
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    append_(
        Segment<THost, TSpec> const & target,
        TSource & source)
    {
SEQAN_CHECKPOINT
        SEQAN_ASSERT_NOT(hasNoHost(target));
        replace(host(target), endPosition(target), endPosition(target), source, TExpand());
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    append_(
        Segment<THost, TSpec> const & target,
        TSource const & source,
        typename Size< Segment<THost, TSpec> >::Type limit)
    {
SEQAN_CHECKPOINT
        SEQAN_ASSERT_NOT(hasNoHost(target));
        replace(host(target), endPosition(target), endPosition(target), source, limit, TExpand()); //??? INSERT
    }
};
//____________________________________________________________________________


template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
append(
    Segment<THost, TSpec> & target,
    TSource & source,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AppendSequenceToSegment_<Tag<TExpand> const>::append_(target, source);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
append(
    Segment<THost, TSpec> & target,
    TSource const & source,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AppendSequenceToSegment_<Tag<TExpand> const>::append_(target, source);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
append(
    Segment<THost, TSpec> & target,
    TSource & source,
    typename Size< Segment<THost, TSpec> >::Type limit,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AppendSequenceToSegment_<Tag<TExpand> const>::append_(target, source, limit);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
append(
    Segment<THost, TSpec> & target,
    TSource const & source,
    typename Size< Segment<THost, TSpec> >::Type limit,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AppendSequenceToSegment_<Tag<TExpand> const>::append_(target, source, limit);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
append(
    Segment<THost, TSpec> const & target,
    TSource & source,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AppendSequenceToSegment_<Tag<TExpand> const>::append_(target, source);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
append(
    Segment<THost, TSpec> const & target,
    TSource const & source,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AppendSequenceToSegment_<Tag<TExpand> const>::append_(target, source);
}

template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
append(
    Segment<THost, TSpec> const & target,
    TSource & source,
    typename Size< Segment<THost, TSpec> >::Type limit,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AppendSequenceToSegment_<Tag<TExpand> const>::append_(target, source, limit);
}
template <typename THost, typename TSpec, typename TSource, typename TExpand>
inline void
append(
    Segment<THost, TSpec> const & target,
    TSource const & source,
    typename Size< Segment<THost, TSpec> >::Type limit,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AppendSequenceToSegment_<Tag<TExpand> const>::append_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
// appendValue
//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct AppendValueToSegment_
{
    template <typename T, typename TValue>
    static inline void
    appendValue_(T & me,
                TValue & /*_value*/)
    {
SEQAN_CHECKPOINT
        insertValue(host(me), endPosition(me), TExpand());
        if (endPosition(me) < length(host(me)) ) //this could be false for some TExpand
        {
            setEndPosition(me, endPosition(me) + 1);
        }
    }
};

//____________________________________________________________________________

template <typename THost, typename TSpec, typename TValue, typename TExpand>
inline void
appendValue(Segment<THost, TSpec> & me,
            TValue const & _value,
            Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AppendValueToSegment_<Tag<TExpand> const>::appendValue_(me, _value);
}
template <typename THost, typename TSpec, typename TValue, typename TExpand>
inline void
appendValue(Segment<THost, TSpec> const & me,
            TValue const & _value,
            Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AppendValueToSegment_<Tag<TExpand> const>::appendValue_(me, _value);
}

//////////////////////////////////////////////////////////////////////////////
// insertValue
//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct InsertValueToSegment_
{
    template <typename T, typename TPosition, typename TValue>
    static inline void
    insertValue_(T & me,
                TPosition pos,
                TValue & /*_value*/)
    {
SEQAN_CHECKPOINT
        insertValue(host(me), beginPosition(me) + pos, TExpand());
        if (endPosition(me) < length(host(me)) ) //this could be false for some TExpand
        {
            setEndPosition(me, endPosition(me) + 1);
        }
    }
};

//____________________________________________________________________________

template <typename THost, typename TSpec, typename TPosition, typename TValue, typename TExpand>
inline void
insertValue(Segment<THost, TSpec> & me,
            TPosition pos,
            TValue const & _value,
            Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    InsertValueToSegment_<Tag<TExpand> const>::insertValue_(me, pos, _value);
}
template <typename THost, typename TSpec, typename TPosition, typename TValue, typename TExpand>
inline void
insertValue(Segment<THost, TSpec> const & me,
            TPosition pos,
            TValue const & _value,
            Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    InsertValueToSegment_<Tag<TExpand> const>::insertValue_(me, pos, _value);
}

//////////////////////////////////////////////////////////////////////////////
// replace
//////////////////////////////////////////////////////////////////////////////

/**
.Function.replace:
..class:Class.Segment
..remarks:If $target$ is a @Class.Segment@ object, then
$limit$ denotes the maximal length of @Function.host.$host(target)$@ after the operation.
..param.target.type:Class.Segment
..param.source.type:Class.Segment
*/

template <typename TExpand>
struct ReplaceSequenceToSegment_
{
    template <typename THost, typename TSpec, typename TSource>
    static inline void
    replace_(
        Segment<THost, TSpec> & target,
        typename Position< Segment<THost, TSpec> >::Type pos_begin,
        typename Position< Segment<THost, TSpec> >::Type pos_end,
        TSource & source)
    {
SEQAN_CHECKPOINT
        SEQAN_ASSERT_NOT(hasNoHost(target));

        typedef Segment<THost, TSpec> Target;

        replace(host(target), beginPosition(target) + pos_begin, beginPosition(target) + pos_end, source, TExpand());

        typename Iterator<Target, Standard>::Type new_end = begin(target, Standard()) + length(target) - pos_end + pos_begin + length(source);
        typename Iterator<THost, Standard>::Type host_end = end(host(target), Standard());
        if (new_end > host_end) new_end = host_end;
        setEnd(target, new_end);
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    replace_(
        Segment<THost, TSpec> & target,
        typename Position< Segment<THost, TSpec> >::Type pos_begin,
        typename Position< Segment<THost, TSpec> >::Type pos_end,
        TSource & source,
        typename Size< Segment<THost, TSpec> >::Type limit)
    {
SEQAN_CHECKPOINT
        SEQAN_ASSERT_NOT(hasNoHost(target));

        typedef Segment<THost, TSpec> Target;

        replace(host(target), beginPosition(target) + pos_begin, beginPosition(target) + pos_end, source, limit, TExpand());

        typename Iterator<Target, Standard>::Type new_end = begin(target, Standard()) + length(target) - pos_end + pos_begin + length(source);
        typename Iterator<THost, Standard>::Type host_end = end(host(target), Standard());
        if (begin(target, Standard()) > host_end) setBegin(target, host_end);
        if (new_end > host_end) new_end = host_end;
        setEnd(target, new_end);
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    replace_(
        Segment<THost, TSpec> const & target,
        typename Position< Segment<THost, TSpec> const>::Type pos_begin,
        typename Position< Segment<THost, TSpec> const>::Type pos_end,
        TSource & source)
    {
SEQAN_CHECKPOINT
        SEQAN_ASSERT_NOT(hasNoHost(target));

        replace(host(target), beginPosition(target) + pos_begin, beginPosition(target) + pos_end, source, TExpand());
    }

    template <typename THost, typename TSpec, typename TSource>
    static inline void
    replace_(
        Segment<THost, TSpec> const & target,
        typename Position< Segment<THost, TSpec> const>::Type pos_begin,
        typename Position< Segment<THost, TSpec> const>::Type pos_end,
        TSource & source,
        typename Size< Segment<THost, TSpec> >::Type limit)
    {
SEQAN_CHECKPOINT
        SEQAN_ASSERT_NOT(hasNoHost(target));

        replace(host(target), beginPosition(target) + pos_begin, beginPosition(target) + pos_end, source, limit, TExpand()); //??? INSERT
    }
};
//____________________________________________________________________________


template <typename THost, typename TSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(
    Segment<THost, TSpec> & target,
    TPositionBegin pos_begin,
    TPositionEnd pos_end,
    TSource & source,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    ReplaceSequenceToSegment_<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}
template <typename THost, typename TSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(
    Segment<THost, TSpec> & target,
    TPositionBegin pos_begin,
    TPositionEnd pos_end,
    TSource const & source,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    ReplaceSequenceToSegment_<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}

template <typename THost, typename TSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(
    Segment<THost, TSpec> & target,
    TPositionBegin pos_begin,
    TPositionEnd pos_end,
    TSource & source,
    typename Size< Segment<THost, TSpec> >::Type limit,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    ReplaceSequenceToSegment_<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}
template <typename THost, typename TSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(
    Segment<THost, TSpec> & target,
    TPositionBegin pos_begin,
    TPositionEnd pos_end,
    TSource const & source,
    typename Size< Segment<THost, TSpec> >::Type limit,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    ReplaceSequenceToSegment_<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}

template <typename THost, typename TSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(
    Segment<THost, TSpec> const & target,
    TPositionBegin pos_begin,
    TPositionEnd pos_end,
    TSource & source,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    ReplaceSequenceToSegment_<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}
template <typename THost, typename TSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(
    Segment<THost, TSpec> const & target,
    TPositionBegin pos_begin,
    TPositionEnd pos_end,
    TSource const & source,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    ReplaceSequenceToSegment_<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}

template <typename THost, typename TSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(
    Segment<THost, TSpec> const & target,
    TPositionBegin pos_begin,
    TPositionEnd pos_end,
    TSource & source,
    typename Size< Segment<THost, TSpec> >::Type limit,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    ReplaceSequenceToSegment_<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}
template <typename THost, typename TSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(
    Segment<THost, TSpec> const & target,
    TPositionBegin pos_begin,
    TPositionEnd pos_end,
    TSource const & source,
    typename Size< Segment<THost, TSpec> >::Type limit,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    ReplaceSequenceToSegment_<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}


//////////////////////////////////////////////////////////////////////////////
///.Function.resize.param.object.type:Class.Segment
///.Function.resize.class:Class.Segment

template <typename THost, typename TSpec, typename TExpand>
inline typename Size< Segment<THost, TSpec> >::Type
resize(
    Segment<THost, TSpec> & me,
    typename Size< Segment<THost, TSpec> >::Type new_length,
    Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT

    typename Size<Segment<THost, TSpec> >::Type me_length = length(me);
    typename Position<THost>::Type me_end_pos = endPosition(me);
    if (new_length > me_length)
    {
        new_length = me_length + resizeSpace(host(me), new_length - me_length, me_end_pos, me_end_pos, tag);
    }
    else if (new_length < me_length)
    {
        new_length = resizeSpace(host(me), 0, me_end_pos - (me_length - new_length), me_end_pos, tag);
    }
    _setLength(me, new_length);
    return new_length;
}

//////////////////////////////////////////////////////////////////////////////
//??? TODO: fill (kopie von resize anpassen)

//////////////////////////////////////////////////////////////////////////////
///.Function.clear.param.object.type:Class.Segment
///.Function.clear.class:Class.Segment

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TLeftSpec, typename TRight>
Segment<TLeftValue, TLeftSpec> const &
operator += (Segment<TLeftValue, TLeftSpec> & left,
             TRight const & right)
{
SEQAN_CHECKPOINT
    append(left, right);
    return left;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight >
inline bool
operator == (Segment<TLeftHost, TLeftSpec> const & left,
            TRight const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<Segment<TLeftHost, TLeftSpec> >::Type _lex(left, right);
    return isEqual(_lex);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight >
inline bool
operator != (Segment<TLeftHost, TLeftSpec> const & left,
            TRight const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<Segment<TLeftHost, TLeftSpec> >::Type _lex(left, right);
    return isNotEqual(_lex);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight>
inline bool
operator < (Segment<TLeftHost, TLeftSpec> const & left,
            TRight const & right)
{
SEQAN_CHECKPOINT
    return isLess(left, right, typename DefaultPrefixOrder<Segment<TLeftHost, TLeftSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight>
inline bool
operator <= (Segment<TLeftHost, TLeftSpec> const & left,
             TRight const & right)
{
SEQAN_CHECKPOINT
    return isLessOrEqual(left, right, typename DefaultPrefixOrder<Segment<TLeftHost, TLeftSpec> >::Type());
}
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight>
inline bool
operator > (Segment<TLeftHost, TLeftSpec> const & left,
            TRight const & right)
{
SEQAN_CHECKPOINT
    return isGreater(left, right, typename DefaultPrefixOrder<Segment<TLeftHost, TLeftSpec> >::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftHost, typename TLeftSpec, typename TRight>
inline bool
operator >= (Segment<TLeftHost, TLeftSpec> const & left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<Segment<TLeftHost, TLeftSpec> >::Type());
}


//////////////////////////////////////////////////////////////////////////////
// stream operators
//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename THost, typename TSpec>
inline TStream &
operator << (TStream & target,
             Segment<THost, TSpec> const & source)
{
SEQAN_CHECKPOINT
    write(target, source);
    return target;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename THost, typename TSpec>
inline TStream &
operator >> (TStream & source,
             Segment<THost, TSpec> & target)
{
SEQAN_CHECKPOINT
    read(source, target);
    return source;
}
template <typename TStream, typename THost, typename TSpec>
inline TStream &
operator >> (TStream & source,
             Segment<THost, TSpec> const & target)
{
SEQAN_CHECKPOINT
    read(source, target);
    return source;
}

//////////////////////////////////////////////////////////////////////////////


} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
