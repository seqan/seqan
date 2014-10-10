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
// Author: Tobias Rausch <rausch@embl.de>
// ==========================================================================

#ifndef SEQAN_HEADER_GRAPH_PROPERTY_H
#define SEQAN_HEADER_GRAPH_PROPERTY_H

namespace seqan
{

// --------------------------------------------------------------------------
// Class ExternalPropertyMap
// --------------------------------------------------------------------------

/*!
 * @class ExternalPropertyMap External Property Map
 * @brief An external property map.
 * @signature  template<TValue, TSpec>
 *             class String<TValue, TSpec>
 * @tparam TValue The type of information stored in the property map.
 * @tparam TSpec The specializing type. Default: @link AllocString @endlink.
 * @section Remarks
 *
 * The external property map is assumed to be an instance of @link String @endlink. It is indexed via @link
 * VertexDescriptor @endlink  or @link Graph#EdgeDescriptor @endlink.
 */

// --------------------------------------------------------------------------
// Class Graph#resizeVertexMap
// --------------------------------------------------------------------------

/*!
 * @fn Graph#resizeVertexMap
 * @brief Initializes a vertex map.
 * @signature void resizeVertexMap(pm g, [, prototype])
 *
 * @param[in,out] pm        An External Property Map. Type: @link ExternalPropertyMap @endlink
 * @param[in]     g         A Graph.
 * @param[in]     prototype An optional prototype that is used for initializing the property map.
 */

template <typename TSpec, typename TPropertyMap>
void resizeVertexMap(TPropertyMap & pm, Graph<TSpec> const & g)
{
    typedef typename Value<TPropertyMap>::Type TValue;
    resize(pm, getIdUpperBound(_getVertexIdManager(g)), TValue(), Generous());
}

template <typename TSpec, typename TPropertyMap, typename TPrototype>
void resizeVertexMap(TPropertyMap & pm, Graph<TSpec> const & g, TPrototype const & prototype)
{
    resize(pm, getIdUpperBound(_getVertexIdManager(g)), prototype, Generous());
}

// --------------------------------------------------------------------------
// Class Graph#resizeEdgeMap
// --------------------------------------------------------------------------

/*!
 * @fn Graph#resizeEdgeMap
 * @brief Initializes an edge map
 * @signature void resizeEdgeMap(pm, g[, prototype]);
 * @param[in,out] pm        An External or Internal Property Map. Types: @link ExternalPropertyMap @endlink,
 *                          @link InternalMap @endlink, @link InternalPointerMap @endlink.
 * @param[in]     g         A Graph.
 * @param[in]     prototype An optional prototype that is used for initializing the property map.
 */

template <typename TSpec, typename TPropertyMap>
void resizeEdgeMap(TPropertyMap & pm, Graph<TSpec> const & g)
{
    typedef typename Value<TPropertyMap>::Type TValue;
    resize(pm, getIdUpperBound(_getEdgeIdManager(g)), TValue(), Generous());
}

template <typename TSpec, typename TPropertyMap, typename TPrototype>
void resizeEdgeMap(TPropertyMap & pm, Graph<TSpec> const & g, TPrototype const & prototype)
{
    resize(pm, getIdUpperBound(_getEdgeIdManager(g)), prototype, Generous());
}

// --------------------------------------------------------------------------
// Class ExternalPropertyMap#assignProperty
// --------------------------------------------------------------------------

/*!
 * @fn ExternalPropertyMap#assignProperty
 * @brief Assigns a property to an item in the property map.
 * @signature void assignProperty(pm, d, val)
 * @param pm  The property map
 * @param d   A vertex or edge descriptor that identifies the item in the property map.
 *            Types: @link VertexDescriptor @endlink, @link Graph#EdgeDescriptor @endlink
 * @param val The new value, where the type of the new value must match the value type of the property map.
 */

template <typename TPropertyMap, typename TDescriptor, typename TValue>
void assignProperty(TPropertyMap & pm, TDescriptor const d, TValue const val)
{
    assignValue(pm, _getId(d), val);
}

// --------------------------------------------------------------------------
// Class ExternalPropertyMap#property
// --------------------------------------------------------------------------

/*!
 * @fn ExternalPropertyMap#property
 * @brief Accesses the property of an item in the property map.
 * @signature TRef property(pm, d)
 * @param pm  The property map.
 * @param d   A vertex or edge descriptor that identifies the item in the property map.
 *            Types: @link VertexDescriptor @endlink, @link Graph#EdgeDescriptor @endlink
 * @return TRef Reference to the item in the property map of type @link Reference @endlink.
 */

template <typename TPropertyMap, typename TDescriptor>
typename Reference<TPropertyMap>::Type
property(TPropertyMap & pm, TDescriptor const d)
{
    return value(pm, _getId(d));
}

template <typename TPropertyMap, typename TDescriptor>
typename Reference<TPropertyMap const>::Type
property(TPropertyMap const & pm, TDescriptor const d)
{
    return value(pm, _getId(d));
}

// --------------------------------------------------------------------------
// Class ExternalPropertyMap#getProperty
// --------------------------------------------------------------------------

/*!
 * @fn ExternalPropertyMap#getProperty
 * @brief Get method for an item's property.
 * @signature TValue getProperty(pm, d)
 * @param pm  The property map.
 * @param d   A vertex or edge descriptor that identifies the item in the property map.
 *            Types: @link VertexDescriptor @endlink, @link Graph#EdgeDescriptor @endlink
 * @return TValue Reference to the item in the property map of type @link GetValue @endlink.
 */

template <typename TPropertyMap, typename TDescriptor>
typename GetValue<TPropertyMap const>::Type
getProperty(TPropertyMap const & pm, TDescriptor const d)
{
    return getValue(pm, _getId(d));
}

// --------------------------------------------------------------------------
// Class InternalMap
// --------------------------------------------------------------------------

/*!
 * @class InternalMap
 * @brief An internal property map using member ids.
 * @signature template<typename TContainer, unsigned const MEMBER_ID>
 *            class InternalMap<TContainer, MEMBER_ID>
 * @tparam TContainer The cargo type which can be determined with @link Cargo @endlink.
 * @tparam MEMBER_ID   An unsigned to specify the position of the member in the cargo. Note: If zero it is assumed
 *                    that the cargo is a simple type (e.g., int). Default <tt>0</tt>.
 * @section Remarks
 *
 * Internal property maps are used to access internal edge cargos.
 */
template <typename TContainer, unsigned const MEMBER_ID = 0>
struct InternalMap
{};

// --------------------------------------------------------------------------
// Metafunction InternalMap#Value
// --------------------------------------------------------------------------

template <typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 1> const>
{
    typedef T1 const Type;
};

template <typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 1> >
{
    typedef T1 Type;
};

template <typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 2> const>
{
    typedef T2 const Type;
};

template <typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 2> >
{
    typedef T2 Type;
};

template <typename T>
struct Value<InternalMap<T, 0> const>
{
    typedef T const Type;
};

template <typename T>
struct Value<InternalMap<T, 0> >
{
    typedef T Type;
};

// --------------------------------------------------------------------------
// Function InternalMap#resizeEdgeMap
// --------------------------------------------------------------------------

template <typename TSpec, typename TContainer, unsigned const MEMBER_ID>
void resizeEdgeMap(InternalMap<TContainer, MEMBER_ID> &,
                   Graph<TSpec> const &)
{}

template <typename TSpec, typename TContainer, unsigned const MEMBER_ID, typename TSource>
void assignEdgeMap(InternalMap<TContainer, MEMBER_ID> &,
                   Graph<TSpec> const &,
                   TSource const &)
{}

// --------------------------------------------------------------------------
// Function InternalMap#assignProperty
// --------------------------------------------------------------------------

/*!
 * @fn InternalMap#assignProperty:
 * @brief Assigns a property to an item in the property map.
 * @signature void assignProperty(pm, d, val)
 * @param pm  The property map
 * @param d   A vertex or edge descriptor that identifies the item in the property map.
 *            Types: @link VertexDescriptor @endlink, @link Graph#EdgeDescriptor @endlink
 * @param val The new value, where thg type of the new value must match the value type of the property map.
 */

template <typename T1, typename T2, typename TEdgeDescriptor, typename TValue>
void assignProperty(InternalMap<Pair<T1, T2>, 1> &,
                    TEdgeDescriptor const e,
                    TValue const val)
{
    (cargo(e)).i1 = val;
}

template <typename T1, typename T2, typename TEdgeDescriptor, typename TValue>
void assignProperty(InternalMap<Pair<T1, T2>, 2> &,
                    TEdgeDescriptor const e,
                    TValue const val)
{
    (cargo(e)).i2 = val;
}

template <typename T, typename TEdgeDescriptor, typename TValue>
void assignProperty(InternalMap<T, 0> &,
                    TEdgeDescriptor const e,
                    TValue const val)
{
    assignCargo(e, val);
}

// --------------------------------------------------------------------------
// Function InternalMap#property
// --------------------------------------------------------------------------

/*!
 * @fn InternalMap#property:
 * @brief Accesses the property of an item in the property map.
 * @signature TRef property(pm, d)
 * @param pm  The property map.
 * @param d   A vertex or edge descriptor that identifies the item in the property map.
 *            Types: @link VertexDescriptor @endlink, @link Graph#EdgeDescriptor @endlink
 * @return TRef Reference to the item in the property map of type @link Reference @endlink.
 */

template <typename T1, typename T2, typename TEdgeDescriptor>
typename Value<InternalMap<Pair<T1, T2>, 2> >::Type &
property(InternalMap<Pair<T1, T2>, 2> &,
         TEdgeDescriptor e)
{
    return (cargo(e)).i2;
}

template <typename T1, typename T2, typename TEdgeDescriptor>
typename Value<InternalMap<Pair<T1, T2>, 2> const>::Type &
property(InternalMap<Pair<T1, T2>, 2> const &,
         TEdgeDescriptor e)
{
    return (cargo(e)).i2;
}

template <typename T1, typename T2, typename TEdgeDescriptor>
typename Value<InternalMap<Pair<T1, T2>, 1> >::Type &
property(InternalMap<Pair<T1, T2>, 1> &,
         TEdgeDescriptor e)
{
    return (cargo(e)).i1;
}

template <typename T1, typename T2, typename TEdgeDescriptor>
typename Value<InternalMap<Pair<T1, T2>, 1> const>::Type &
property(InternalMap<Pair<T1, T2>, 1> const &,
         TEdgeDescriptor e)
{
    return (cargo(e)).i1;
}

template <typename T, typename TEdgeDescriptor>
typename Value<InternalMap<T, 0> >::Type &
property(InternalMap<T, 0>&,
         TEdgeDescriptor e)
{
    return cargo(e);
}

template <typename T, typename TEdgeDescriptor>
typename Value<InternalMap<T, 0> const>::Type &
property(InternalMap<T, 0> const &,
         TEdgeDescriptor e)
{
    return cargo(e);
}

// --------------------------------------------------------------------------
// Function InternalMap#getProperty
// --------------------------------------------------------------------------

/*!
 * @fn InternalMap#getProperty
 * @brief Get method for an item's property.
 * @signature TValue getProperty(pm, d)
 * @param pm  The property map.
 * @param d   A vertex or edge descriptor that identifies the item in the property map.
 *            Types: @link VertexDescriptor @endlink, @link Graph#EdgeDescriptor @endlink
 * @return TValue Reference to the item in the property map of type @link GetValue @endlink.
 */

template <typename T1, typename T2, typename TEdgeDescriptor>
typename Value<InternalMap<Pair<T1, T2>, 1> >::Type
getProperty(InternalMap<Pair<T1, T2>, 1> const &,
            TEdgeDescriptor e)
{
    return (getCargo(e)).i1;
}

template <typename T1, typename T2, typename TEdgeDescriptor>
typename Value<InternalMap<Pair<T1, T2>, 1> >::Type
getProperty(InternalMap<Pair<T1, T2>, 1> &,
            TEdgeDescriptor e)
{
    return (getCargo(e)).i1;
}

template <typename T1, typename T2, typename TEdgeDescriptor>
typename Value<InternalMap<Pair<T1, T2>, 2> >::Type
getProperty(InternalMap<Pair<T1, T2>, 2> const &,
            TEdgeDescriptor e)
{
    return (getCargo(e)).i2;
}

template <typename T1, typename T2, typename TEdgeDescriptor>
typename Value<InternalMap<Pair<T1, T2>, 2> >::Type
getProperty(InternalMap<Pair<T1, T2>, 2> &,
            TEdgeDescriptor e)
{
    return (getCargo(e)).i2;
}

template <typename T, typename TEdgeDescriptor>
typename Value<InternalMap<T, 0> >::Type
getProperty(InternalMap<T, 0> const &,
            TEdgeDescriptor e)
{
    return getCargo(e);
}

template <typename T, typename TEdgeDescriptor>
typename Value<InternalMap<T, 0> >::Type
getProperty(InternalMap<T, 0> &,
            TEdgeDescriptor e)
{
    return getCargo(e);
}

// --------------------------------------------------------------------------
// Class InternalPointerMap
// --------------------------------------------------------------------------

/*!
 * @class InternalPointerMap
 * @brief An internal property map using pointer to members.
 * @signature template <typename TPropmap, TPropmap const Instance>
 *            class InternalPointerMap<TPropmap, Instance>;
 * tparam TPropmap  A pointer to member type.
 * @tparam Instance A pointer to a member of type TPropmap.
 * @section Remark
 *
 * Internal property maps are used to access internal edge cargos.
 */

template <typename TPropmap, TPropmap const Instance>
struct InternalPointerMap
{};

// --------------------------------------------------------------------------
// Metafunction InternalPointerMap#Value
// --------------------------------------------------------------------------

template <typename TClass, typename TValue, TValue TClass::* TPMember>
struct Value<InternalPointerMap<TValue TClass::*, TPMember> const>
{
    typedef TValue const Type;
};

template <typename TClass, typename TValue, TValue TClass::* TPMember>
struct Value<InternalPointerMap<TValue TClass::*, TPMember> >
{
    typedef TValue Type;
};

// --------------------------------------------------------------------------
// Function InternalPointerMap#resizeEdgeMap
// --------------------------------------------------------------------------

template <typename TPropMap, TPropMap const INSTANCE, typename TSpec>
void resizeEdgeMap(InternalPointerMap<TPropMap, INSTANCE> &, Graph<TSpec> const &)
{}

// --------------------------------------------------------------------------
// Function InternalPointerMap#assignProperty
// --------------------------------------------------------------------------

/*!
 * @fn InternalPointerMap#assignProperty:
 * @brief Assigns a property to an item in the property map.
 * @signature void assignProperty(pm, d, val)
 * @param pm  The property map
 * @param d   A vertex or edge descriptor that identifies the item in the property map.
 *            Types: @link VertexDescriptor @endlink, @link Graph#EdgeDescriptor @endlink
 * @param val The new value, where thg type of the new value must match the value type of the property map.
 */
template <typename TClass, typename TValue, TValue TClass::* TPMember, typename TEdgeDescriptor>
void assignProperty(InternalPointerMap<TValue TClass::*, TPMember> &,
                    TEdgeDescriptor const e,
                    TValue const val)
{
    (cargo(e)).*TPMember = val;
}

// --------------------------------------------------------------------------
// Function InternalPointerMap#property
// --------------------------------------------------------------------------

/*!
 * @fn InternalPointerMap#property:
 * @brief Accesses the property of an item in the property map.
 * @signature TRef property(pm, d)
 * @param pm  The property map.
 * @param d   A vertex or edge descriptor that identifies the item in the property map.
 *            Types: @link VertexDescriptor @endlink, @link Graph#EdgeDescriptor @endlink
 * @return TRef Reference to the item in the property map of type @link Reference @endlink.
 */

template <typename TClass, typename TValue, TValue TClass::* TPMember, typename TEdgeDescriptor>
typename Value<InternalPointerMap<TValue TClass::*, TPMember> >::Type &
property(InternalPointerMap<TValue TClass::*, TPMember>&,
         TEdgeDescriptor const e)
{
    return (cargo(e)).*TPMember;
}

template <typename TClass, typename TValue, TValue TClass::* TPMember, typename TEdgeDescriptor>
typename Value<InternalPointerMap<TValue TClass::*, TPMember> const>::Type &
property(InternalPointerMap<TValue TClass::*, TPMember> const &,
         TEdgeDescriptor const e)
{
    return (cargo(e)).*TPMember;
}

// --------------------------------------------------------------------------
// Function InternalPointerMap#getProperty
// --------------------------------------------------------------------------

/*!
 * @fn InternalPointerMap#getProperty
 * @brief Get method for an item's property.
 * @signature TValue getProperty(pm, d)
 * @param pm  The property map.
 * @param d   A vertex or edge descriptor that identifies the item in the property map.
 *            Types: @link VertexDescriptor @endlink, @link Graph#EdgeDescriptor @endlink
 * @return TValue Reference to the item in the property map of type @link GetValue @endlink.
 */

template <typename TClass, typename TValue, TValue TClass::* TPMember, typename TEdgeDescriptor>
typename Value<InternalPointerMap<TValue TClass::*, TPMember> >::Type
getProperty(InternalPointerMap<TValue TClass::*, TPMember> const &,
            TEdgeDescriptor const e)
{
    return (getCargo(e)).*TPMember;
}

template <typename TClass, typename TValue, TValue TClass::* TPMember, typename TEdgeDescriptor>
typename Value<InternalPointerMap<TValue TClass::*, TPMember> >::Type
getProperty(InternalPointerMap<TValue TClass::*, TPMember> &,
            TEdgeDescriptor const e)
{
    return (getCargo(e)).*TPMember;
}

// --------------------------------------------------------------------------
// Function Graph#assignVertexMap()
// --------------------------------------------------------------------------

/*!
 * @fn Graph#assignVertexMap
 * @brief Initializes a vertex map with values of an array.
 * @signature void assignVertexMap(g, pm, prop)
 *
 * @param[out] pm   An External Property Map. Types: @link ExternalPropertyMap @endlink
 * @param[in]  g    A Graph.
 * @param[in]  prop An array with properties that are to be assigned to the items in the property map.
 *
 * @section Remarks
 *
 * For every vertex descriptor there must be an entry in the array.
 */

template <typename TSpec, typename TPropertyMap, typename TProperties>
void assignVertexMap(TPropertyMap & pm,
                     Graph<TSpec> const & g,
                     TProperties const & prop)
{
    resize(pm, getIdUpperBound(_getVertexIdManager(g)), Generous());
    typedef Graph<TSpec> TGraph;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    for (; !atEnd(it); goNext(it))
        assignProperty(pm, getValue(it), getValue(prop, _getId(value(it))));
}

// --------------------------------------------------------------------------
// Function Graph#assignEdgeMap()
// --------------------------------------------------------------------------

/*!
 * @fn Graph#assignEdgeMap
 * @brief Initializes a vertex map with values of an array.
 * @signature void assignEdgeMap(g, pm, prop)
 *
 * @param[in]  pm   An External or Internal Property Map. Types: @link ExternalPropertyMap @endlink, @link InternalMap @endlink,
 *                  @link InternalPointerMap @endlink.
 * @param[out] prop An array with properties that are to be assigned to the items in the property map.
 * @param[in]  g    A Graph.
 *
 * @section Remarks
 *
 * For every edge id there must be an entry in the array.
 */
template <typename TSpec, typename TPropertyMap, typename TProperties>
void assignEdgeMap(TPropertyMap & pm,
                   Graph<TSpec> const & g,
                   TProperties const & prop)
{
    resize(pm, getIdUpperBound(_getEdgeIdManager(g)), Generous());
    typedef Graph<TSpec> TGraph;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    TEdgeIterator it(g);
    for (; !atEnd(it); goNext(it))
        assignProperty(pm, *it, prop[_getId(*it)]);
}

}  // namespace seqan

#endif //#ifndef SEQAN_HEADER_...
