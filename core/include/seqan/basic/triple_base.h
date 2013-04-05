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
// Triple base class.
// ==========================================================================

// TODO(holtgrew): What about move construction? Useful for pairs of strings and such. Tricky to implement since ints have no move constructor, for example.

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_TRIPLE_BASE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_TRIPLE_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.Triple:
..cat:Aggregates
..concept:Concept.AggregateConcept
..summary:Stores three arbitrary objects.
..signature:Triple<T1[, T2[, T3[, TSpec]]]>
..param.T1:The type of the first object.
..param.T2:The type of the second object.
...default:$T1$
..param.T3:The type of the third object.
...default:$T1$
..param.TSpec:The specializing type.
...default:$void$, no packing (faster access).

.Memfunc.Triple#Triple:
..class:Class.Triple
..summary:Constructor
..signature:Triple<T1, T2, T3[, TSpec]> ()
..signature:Triple<T1, T2, T3[, TSpec]> (triple)
..signature:Triple<T1, T2, T3[, TSpec]> (i1, i2, i3)
..param.triple:Other Triple object. (copy constructor)
..param.i1:T1 object.
..param.i2:T2 object.
..param.i3:T3 object.

.Memvar.Triple#i1:
..class:Class.Triple
..summary:T1 object

.Memvar.Triple#i2:
..class:Class.Triple
..summary:T2 object

.Memvar.Triple#i3:
..class:Class.Triple
..summary:T3 object
..include:seqan/basic.h
*/

template <typename T1, typename T2 = T1, typename T3 = T1, typename TSpec = void>
struct Triple
{
    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

    T1 i1;
    T2 i2;
    T3 i3;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    inline Triple() : i1(T1()), i2(T2()), i3(T3()) {}
    
    inline Triple(Triple const & _p)
            : i1(_p.i1), i2(_p.i2), i3(_p.i3) {}
    
    inline Triple(T1 const & _i1, T2 const & _i2, T3 const & _i3)
            : i1(_i1), i2(_i2), i3(_i3) {}
    
    template <typename T1_, typename T2_, typename T3_, typename TSpec__>
    inline Triple(Triple<T1_, T2_, T3_, TSpec__> const & _p)
            : i1(getValueI1(_p)), i2(getValueI2(_p)), i3(getValueI3(_p)) {}

    // TODO(holtgrew): Move comparison operators to global functions?
    inline bool
    operator==(Triple const & other) const
    {
        return i1 == other.i1 && i2 == other.i2 && i3 == other.i3;
    }
    
    inline bool
    operator<(Triple const & other) const
    {
        if (i1 < other.i1)
            return true;
        if (i1 == other.i1 && i2 < other.i2)
            return true;
        if (i1 == other.i1 && i2 == other.i2 && i3 < other.i3)
                return true;
        return false;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// -----------------------------------------------------------------------
// Metafunction LENGTH
// -----------------------------------------------------------------------

///.Metafunction.LENGTH.param.T.type:Class.Triple
///.Metafunction.LENGTH.class:Class.Triple

template <typename T1, typename T2, typename T3, typename TSpec>
struct LENGTH<Triple<T1, T2, T3, TSpec> >
{
    enum { VALUE = 3 };
};

// Const variant is mapped to non-const.

// -----------------------------------------------------------------------
// Metafunction Value
// -----------------------------------------------------------------------

///.Metafunction.Value.param.T.type:Class.Triple
///.Metafunction.Value.class:Class.Triple

template <typename T1, typename T2, typename T3, typename TSpec>
struct Value<Triple<T1, T2, T3, TSpec>, 1>
{
    typedef T1 Type;
};

template <typename T1, typename T2, typename T3, typename TSpec>
struct Value<Triple<T1, T2, T3, TSpec>, 2>
{
    typedef T2 Type;
};

template <typename T1, typename T2, typename T3, typename TSpec>
struct Value<Triple<T1, T2, T3, TSpec>, 3 >
{
    typedef T3 Type;
};

// -----------------------------------------------------------------------
// Metafunction Spec
// -----------------------------------------------------------------------

///.Metafunction.Spec.param.T.type:Class.Triple
///.Metafunction.Spec.class:Class.Triple

template <typename T1, typename T2, typename T3, typename TSpec>
struct Spec<Triple<T1, T2, T3, TSpec> >
{
    typedef TSpec Type;
};

// ============================================================================
// Functions
// ============================================================================

// -----------------------------------------------------------------------
// Function operator<<();  Stream Output.
// -----------------------------------------------------------------------

template <typename T1, typename T2, typename T3, typename TSpec>
std::ostream & operator<<(std::ostream & out, Triple<T1,T2,T3,TSpec> const & t)
{
    out << "< " << getValueI1(t) << " , " << getValueI2(t) << " , " << getValueI3(t) << " >";
    return out;
}

// -----------------------------------------------------------------------
// Function getValueIX()
// -----------------------------------------------------------------------

template <typename T1, typename T2, typename T3, typename TSpec>
inline T1
getValueI1(Triple<T1, T2, T3, TSpec> const & triple)
{
    return triple.i1;
}

template <typename T1, typename T2, typename T3, typename TSpec>
inline T2
getValueI2(Triple<T1, T2, T3, TSpec> const & triple)
{
    return triple.i2;
}

template <typename T1, typename T2, typename T3, typename TSpec>
inline T3
getValueI3(Triple<T1, T2, T3, TSpec> const & triple)
{
    return triple.i3;
}

// -----------------------------------------------------------------------
// Function assignValueIX()
// -----------------------------------------------------------------------

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void assignValueI1(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    triple.i1 = _i;
}

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void assignValueI2(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    triple.i2 = _i;
}

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void assignValueI3(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    triple.i3 = _i;
}

// -----------------------------------------------------------------------
// Function setValueIX()
// -----------------------------------------------------------------------

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void setValueI1(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    set(triple.i1, _i);
}

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void setValueI2(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    set(triple.i2, _i);
}

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void setValueI3(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    set(triple.i3, _i);
}

// -----------------------------------------------------------------------
// Function moveValueIX()
// -----------------------------------------------------------------------

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void moveValueI1(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    move(triple.i1, _i);
}

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void moveValueI2(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    move(triple.i2, _i);
}

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void moveValueI3(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    move(triple.i3, _i);
}

// -----------------------------------------------------------------------
// Function operator<()
// -----------------------------------------------------------------------

template <
    typename L1, typename L2, typename L3, typename LPack, 
    typename R1, typename R2, typename R3, typename RPack>
inline bool
operator<(Triple<L1, L2, L3, LPack> const & _left,
          Triple<R1, R2, R3, RPack> const & _right)
{
    return _left.i1 < _right.i1 || (_left.i1 == _right.i1 && _left.i2 < _right.i2) || (_left.i1 == _right.i1 && _left.i2 == _right.i2 && _left.i3 < _right.i3);
}

// -----------------------------------------------------------------------
// Function operator>()
// -----------------------------------------------------------------------

template <
    typename L1, typename L2, typename L3, typename LPack, 
    typename R1, typename R2, typename R3, typename RPack>
inline bool
operator>(Triple<L1, L2, L3, LPack> const & _left,
          Triple<R1, R2, R3, RPack> const & _right)
{
    return _left.i1 > _right.i1 || (_left.i1 == _right.i1 && _left.i2 > _right.i2) || (_left.i1 == _right.i1 && _left.i2 == _right.i2 && _left.i3 > _right.i3);
}

// -----------------------------------------------------------------------
// Function operator<=()
// -----------------------------------------------------------------------

template <
    typename L1, typename L2, typename L3, typename LPack, 
    typename R1, typename R2, typename R3, typename RPack>
inline bool
operator<=(Triple<L1, L2, L3, LPack> const & _left,
           Triple<R1, R2, R3, RPack> const & _right)
{
    return !operator>(_left, _right);
}

// -----------------------------------------------------------------------
// Function operator==()
// -----------------------------------------------------------------------

template <
    typename L1, typename L2, typename L3, typename LPack, 
    typename R1, typename R2, typename R3, typename RPack>
inline bool
operator==(Triple<L1, L2, L3, LPack> const & _left,
           Triple<R1, R2, R3, RPack> const & _right)
{
    return _left.i1 == _right.i1 && _left.i2 == _right.i2 && _left.i3 == _right.i3;
}

// -----------------------------------------------------------------------
// Function operator>()
// -----------------------------------------------------------------------

template <
    typename L1, typename L2, typename L3, typename LPack, 
    typename R1, typename R2, typename R3, typename RPack>
inline bool
operator>=(Triple<L1, L2, L3, LPack> const & _left,
           Triple<R1, R2, R3, RPack> const & _right)
{
    return !operator<(_left, _right);
}

// -----------------------------------------------------------------------
// Function operator!=()
// -----------------------------------------------------------------------

template <
    typename L1, typename L2, typename L3, typename LPack, 
    typename R1, typename R2, typename R3, typename RPack>
inline bool
operator!=(Triple<L1, L2, L3, LPack> const & _left,
           Triple<R1, R2, R3, RPack> const & _right)
{
    return !operator==(_left, _right);
}
}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_TRIPLE_BASE_H_
