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
// Implements a dynamic prefix sum data structure.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_MISC_MISC_DYNAMIC_PREFIX_SUM_H_
#define EXTRAS_INCLUDE_SEQAN_MISC_MISC_DYNAMIC_PREFIX_SUM_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <unsigned B_FACTOR>
struct DpsHelper_
{
    static const unsigned MAX_CHILD_NUMBER;
    static const unsigned MAX_KEY_NUMBER;

    static const unsigned MIN_CHILD_NUMBER;
    static const unsigned MIN_KEY_NUMBER;
};

template <unsigned B_FACTOR>
const unsigned DpsHelper_<B_FACTOR>::MAX_CHILD_NUMBER = B_FACTOR << 1;
template <unsigned B_FACTOR>
const unsigned DpsHelper_<B_FACTOR>::MAX_KEY_NUMBER = (B_FACTOR << 1) - 1;

template <unsigned B_FACTOR>
const unsigned DpsHelper_<B_FACTOR>::MIN_CHILD_NUMBER = B_FACTOR;
template <unsigned B_FACTOR>
const unsigned DpsHelper_<B_FACTOR>::MIN_KEY_NUMBER = B_FACTOR - 1;

// ----------------------------------------------------------------------------
// Class DpsTreeNode_
// ----------------------------------------------------------------------------

template <typename TValue, typename TCargo, unsigned B_FACTOR>
class DpsTreeNode_
{
public:

    typedef String<TValue>        TKeyTab;
    typedef String<TCargo>        TCargoTab;
    typedef String<DpsTreeNode_*> TChildTab;

    TKeyTab   _keyTab;
    TCargoTab _cargoTab;
    TChildTab _childTab;
    bool      _isLeaf;

    DpsTreeNode_() : _keyTab(), _cargoTab(), _childTab(), _isLeaf(false)
    {
        // Set the capacities.
        reserve(_keyTab, DpsHelper_<B_FACTOR>::MAX_KEY_NUMBER, Exact());
        reserve(_cargoTab, DpsHelper_<B_FACTOR>::MAX_KEY_NUMBER, Exact());
        reserve(_childTab, DpsHelper_<B_FACTOR>::MAX_CHILD_NUMBER, Exact());
    }

    // Copy Constructor
    DpsTreeNode_(DpsTreeNode_ const & other) : _keyTab(other._keyTab),
                                               _cargoTab(other._cargoTab),
                                               _childTab(other._childTab),
                                               _isLeaf(other._isLeaf)
    {}

    // Assignment Operator
    DpsTreeNode_ &
    operator=(DpsTreeNode_ const & other)
    {
        if (this != &other)
        {
            _keyTab   = other._keyTab;
            _cargoTab = other._cargoTab;
            _childTab = other._childTab;
            _isLeaf   = other._isLeaf;
        }
        return *this;
    }
};

// ----------------------------------------------------------------------------
// Class DynamicPrefixSumTree
// ----------------------------------------------------------------------------

template <typename TKeyValue, typename TCargo, unsigned B_FACTOR = 64, typename TSpec = Default>
class DynamicPrefixSumTree;

template <typename TKeyValue, typename TCargo, unsigned B_FACTOR>
class DynamicPrefixSumTree<TKeyValue, TCargo, B_FACTOR, Default>
{
public:

    typedef typename Value<DynamicPrefixSumTree>::Type TNode;
    typedef Allocator<SinglePool<sizeof(TNode)> > TNodeAllocator;

    static const TNode* NIL;

    TNodeAllocator _nodeAllocator;  // Allocator used to allocate/deallocate nodes.

    TNode * rootPtr;

    DynamicPrefixSumTree() : rootPtr(NULL)
    {
        rootPtr = _createNode(*this);
        _makeLeaf(*rootPtr);
    }

    // Copy Constructor
    DynamicPrefixSumTree(DynamicPrefixSumTree const & other) : rootPtr(other.rootPtr)
    {}

    // Assignment Operator
    DynamicPrefixSumTree &
    operator=(DynamicPrefixSumTree const & other)
    {
        if (this != &other)
            rootPtr = other.rootPtr;
        return *this;
    }
};

template <typename TKeyValue, typename TCargo, unsigned B_FACTOR>
const typename DynamicPrefixSumTree<TKeyValue, TCargo, B_FACTOR, Default>::TNode* DynamicPrefixSumTree<TKeyValue, TCargo, B_FACTOR, Default>::NIL = (typename DynamicPrefixSumTree<TKeyValue, TCargo, B_FACTOR, Default>::TNode *) 0;

// ----------------------------------------------------------------------------
// Class PrefixSumFunctor_
// ----------------------------------------------------------------------------

template <typename TNode>
struct PrefixSumFunctor_{};

template <typename TKey, typename TCargo, unsigned B_FACTOR>
struct PrefixSumFunctor_<DpsTreeNode_<TKey, TCargo, B_FACTOR> >
{
    typedef DpsTreeNode_<TKey, TCargo, B_FACTOR> TNode;

    TCargo currSum;

    PrefixSumFunctor_() : currSum(0)
    {}

    template <typename TPos>
    inline void
    operator()(TNode const & node, TPos idx)
    {
        currSum += getCargo(node, idx);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TCargo, unsigned B_FACTOR, typename TSpec>
struct Value<DynamicPrefixSumTree<TValue, TCargo, B_FACTOR, TSpec> >
{
    typedef DpsTreeNode_<TValue, TCargo, B_FACTOR> Type;
};

template <typename TValue, typename TCargo, unsigned B_FACTOR, typename TSpec>
struct Value<DynamicPrefixSumTree<TValue, TCargo, B_FACTOR, TSpec> const>
{
    typedef DpsTreeNode_<TValue, TCargo, B_FACTOR> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Cargo
// ----------------------------------------------------------------------------

template <typename TKeyValue, typename TCargo, unsigned B_FACTOR, typename TSpec>
struct Cargo<DynamicPrefixSumTree<TKeyValue, TCargo, B_FACTOR, TSpec> >
{
    typedef TCargo Type;
};

template <typename TKeyValue, typename TCargo, unsigned B_FACTOR, typename TSpec>
struct Cargo<DynamicPrefixSumTree<TKeyValue, TCargo, B_FACTOR, TSpec> const>
{
    typedef TCargo const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _createNode()
// ----------------------------------------------------------------------------

// Allocates new memory for an internal node and default constructs it.
// Returns pointer to the newly constructed value.
template <typename TKeyValue, typename TCargo, unsigned B_FACTOR, typename TSpec>
inline typename Value<DynamicPrefixSumTree<TKeyValue, TCargo, B_FACTOR, TSpec> >::Type *
_createNode(DynamicPrefixSumTree<TKeyValue, TCargo, B_FACTOR, TSpec> & tree)
{
    typedef DynamicPrefixSumTree<TKeyValue, TCargo, B_FACTOR, TSpec> TTree;
    typedef typename Value<TTree>::Type TNode;

    TNode *tmp;
    allocate(tree._nodeAllocator, tmp, 1);
    return new (tmp) TNode();
}

// ----------------------------------------------------------------------------
// Function _leafState()
// ----------------------------------------------------------------------------

// Returns the isLeaf flag.
template <typename TValue, typename TCargo, unsigned B_FACTOR>
inline bool &
_leafState(DpsTreeNode_<TValue, TCargo, B_FACTOR> & node)
{
    return node._isLeaf;
}

// Returns the isLeaf flag.
template <typename TValue, typename TCargo, unsigned B_FACTOR>
inline bool const &
_leafState(DpsTreeNode_<TValue, TCargo, B_FACTOR> const & node)
{
    return node._isLeaf;
}

// ----------------------------------------------------------------------------
// Function _makeLeaf()
// ----------------------------------------------------------------------------

// Sets the leaf flag of the node to true.
template <typename TValue, typename TCargo, unsigned B_FACTOR>
inline void
_makeLeaf(DpsTreeNode_<TValue, TCargo, B_FACTOR> & node)
{
    node._isLeaf = true;
}

// ----------------------------------------------------------------------------
// Function _makeInnerNode()
// ----------------------------------------------------------------------------

// Sets the leaf flag of the node to false.
template <typename TValue, typename TCargo, unsigned B_FACTOR>
inline void
_makeInnerNode(DpsTreeNode_<TValue, TCargo, B_FACTOR> & node)
{
    node._isLeaf = false;
}

// ----------------------------------------------------------------------------
// Function _getRoot()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR, typename TSpec>
inline typename Value<DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> >::Type &
_getRoot(DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> & tree)
{
    return *tree.rootPtr;
}

template <typename TKey, typename TCargo, unsigned B_FACTOR, typename TSpec>
inline typename Value<DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> const>::Type &
_getRoot(DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> const & tree)
{
    return *tree.rootPtr;
}

// ----------------------------------------------------------------------------
// Function _isFull()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR>
inline bool
_isFull(DpsTreeNode_<TKey, TCargo, B_FACTOR> const & node)
{
    return length(_getKeyTable(node)) == DpsHelper_<B_FACTOR>::MAX_KEY_NUMBER;
}

// ----------------------------------------------------------------------------
// Function _getKeyTable()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR>
inline String<TKey> &
_getKeyTable(DpsTreeNode_<TKey, TCargo, B_FACTOR> & node)
{
    return node._keyTab;
}

template <typename TKey, typename TCargo, unsigned B_FACTOR>
inline String<TKey> const &
_getKeyTable(DpsTreeNode_<TKey, TCargo, B_FACTOR> const & node)
{
    return node._keyTab;
}

// ----------------------------------------------------------------------------
// Function _getCargoTable()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR>
inline String<TCargo> &
_getCargoTable(DpsTreeNode_<TKey, TCargo, B_FACTOR> & node)
{
    return node._cargoTab;
}

template <typename TKey, typename TCargo, unsigned B_FACTOR>
inline String<TCargo> const &
_getCargoTable(DpsTreeNode_<TKey, TCargo, B_FACTOR> const & node)
{
    return node._cargoTab;
}

// ----------------------------------------------------------------------------
// Function _getChildTable()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR>
inline String<DpsTreeNode_<TKey, TCargo, B_FACTOR> *>  &
_getChildTable(DpsTreeNode_<TKey, TCargo, B_FACTOR> & node)
{
    return node._childTab;
}

template <typename TKey, typename TCargo, unsigned B_FACTOR>
inline String<DpsTreeNode_<TKey, TCargo, B_FACTOR> *> const &
_getChildTable(DpsTreeNode_<TKey, TCargo, B_FACTOR> const & node)
{
    return node._childTab;
}

// ----------------------------------------------------------------------------
// Function _splitChild()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR, typename TSpec, typename TPos>
inline void
_splitChild(DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> & tree,
            DpsTreeNode_<TKey, TCargo, B_FACTOR> & target,
            DpsTreeNode_<TKey, TCargo, B_FACTOR> & parent,
            TPos const & childIdx)
{
    typedef DpsTreeNode_<TKey, TCargo, B_FACTOR>                TNode;
    typedef typename TNode::TCargoTab                           TCargoTab;
    typedef typename Iterator<TCargoTab, Standard>::Type        TCargoIter;
    typedef typename Size<TCargoTab>::Type                      TSize;

    TNode* tmpPtr = _createNode(tree);  // Create new node in memory pool.
    _leafState(*tmpPtr) = _leafState(target);

    // We need to set the partial sum first.

    // Move the right half of the target node to the beginning of the new split node.
    resize(_getKeyTable(*tmpPtr), DpsHelper_<B_FACTOR>::MIN_KEY_NUMBER);
    resize(_getCargoTable(*tmpPtr), DpsHelper_<B_FACTOR>::MIN_KEY_NUMBER);
    arrayMoveForward(begin(_getKeyTable(target), Standard()) + B_FACTOR, end(_getKeyTable(target), Standard()),
                     begin(_getKeyTable(*tmpPtr), Standard()));
    arrayMoveForward(begin(_getCargoTable(target), Standard()) + B_FACTOR, end(_getCargoTable(target), Standard()),
                     begin(_getCargoTable(*tmpPtr), Standard()));

    // Update the prefix sums for the split node.
    TCargo cargoBeforeSplit = _getCargoTable(target)[DpsHelper_<B_FACTOR>::MIN_KEY_NUMBER];
    TCargoIter cargoItSplit = begin(_getCargoTable(*tmpPtr), Standard());
    TCargoIter cargoItSplitEnd = end(_getCargoTable(*tmpPtr), Standard());

    for (; cargoItSplit != cargoItSplitEnd; ++cargoItSplit)
        *cargoItSplit -= cargoBeforeSplit;

    // Move the right half of the target node's child table to the left half of the split node's child table.
    if (!isLeaf(target))
    {
        resize(_getChildTable(*tmpPtr), DpsHelper_<B_FACTOR>::MIN_CHILD_NUMBER);
        arrayMoveForward(begin(_getChildTable(target), Standard()) + B_FACTOR, end(_getChildTable(target), Standard()),
                         begin(_getChildTable(*tmpPtr), Standard()));
        resize(_getChildTable(target), DpsHelper_<B_FACTOR>::MIN_CHILD_NUMBER);
    }

    // Move the pivot element to the parent and update tables accordingly.
    TSize oldLength = length(_getChildTable(parent));
    resize(_getChildTable(parent), oldLength + 1);  // Make space for the new element.
    arrayMoveBackward(begin(_getChildTable(parent), Standard()) + childIdx + 1,
                      begin(_getChildTable(parent), Standard()) + oldLength,
                      begin(_getChildTable(parent), Standard()) + childIdx + 2);

    resize(_getKeyTable(parent), oldLength);  // Make space for the new element.
    arrayMoveBackward(begin(_getKeyTable(parent), Standard()) + childIdx,
                      begin(_getKeyTable(parent), Standard()) + oldLength - 1,
                      begin(_getKeyTable(parent), Standard()) + childIdx + 1);

    resize(_getCargoTable(parent), oldLength);  // Make space for the new element.
    arrayMoveBackward(begin(_getCargoTable(parent), Standard()) + childIdx ,
                      begin(_getCargoTable(parent), Standard()) + oldLength - 1,
                      begin(_getCargoTable(parent), Standard()) + childIdx + 1);

    _getKeyTable(parent)[childIdx] = _getKeyTable(target)[B_FACTOR - 1];
    _getCargoTable(parent)[childIdx] = _getCargoTable(target)[B_FACTOR - 1];
    _getChildTable(parent)[childIdx + 1] = tmpPtr;

    // Shrink the elements of target.
    resize(_getKeyTable(target), DpsHelper_<B_FACTOR>::MIN_KEY_NUMBER);
    resize(_getCargoTable(target), DpsHelper_<B_FACTOR>::MIN_KEY_NUMBER);
}

// ----------------------------------------------------------------------------
// Function _updateCargoRight()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR, typename TPos, typename TCargo2>
inline void
_updateCargoRight(DpsTreeNode_<TKey, TCargo, B_FACTOR> & node,
                  TPos const & idx,
                  TCargo2 const & newCargo)
{
    typedef DpsTreeNode_<TKey, TCargo, B_FACTOR>         TNode;
    typedef typename TNode::TCargoTab                    TCargoTab;
    typedef typename Iterator<TCargoTab, Standard>::Type TCargoTabIter;

    // Update the cargos of all right values in O(B_FACTOR) time.
    TCargoTabIter itCargo = begin(_getCargoTable(node), Standard()) + idx;
    for (; itCargo != end(_getCargoTable(node), Standard()); ++itCargo)
        *itCargo += newCargo;
}

// ----------------------------------------------------------------------------
// Function _insertNonFull()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR, typename TSpec, typename TKey2, typename TCargo2>
inline bool
_insertNonFull(DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> & tree,
               DpsTreeNode_<TKey, TCargo, B_FACTOR> & node,
               TKey2 const & newKey,
               TCargo2 const & newCargo)
{
    typedef DpsTreeNode_<TKey, TCargo, B_FACTOR>       TNode;
    typedef typename TNode::TKeyTab                    TKeyTab;
    typedef typename Size<TKeyTab>::Type               TSize;
    typedef typename Iterator<TKeyTab, Standard>::Type TKeyTabIter;

    // Found lower bound.
    TKeyTabIter it = std::lower_bound(begin(_getKeyTable(node), Standard()), end(_getKeyTable(node), Standard()), newKey);

    // Key exists already.
    if (*it == static_cast<TKey>(newKey))
        return false;

    TSize foundPos = it - begin(_getKeyTable(node), Standard());

    if (isLeaf(node))
    {
        _updateCargoRight(node, foundPos, newCargo);
        insertValue(_getKeyTable(node), foundPos, newKey);  // Insert the new key in leaf.
        if (foundPos == 0)
            insertValue(_getCargoTable(node), foundPos, newCargo);  // Insert the new cargo.
        else
            insertValue(_getCargoTable(node), foundPos, newCargo + getCargo(node, foundPos - 1));  // Update with left neighbor before insertion.
        return true;
    }

    TNode& child = getChild(node, foundPos);  // Get the child node affected by the insertion.
    if (_isFull(child))  // If already full, we split the child.
    {
        _splitChild(tree, child, node, foundPos);
        // We need to update the cargo here.
        if (foundPos > 0)
            _getCargoTable(node)[foundPos] += getCargo(node, foundPos -1);
        if (getKey(node, foundPos) < static_cast<TKey>(newKey))
            ++foundPos;
        _updateCargoRight(node, foundPos, newCargo);
        return _insertNonFull(tree, getChild(node, foundPos), newKey, newCargo);  // Go into subtree left of split.
    }
    _updateCargoRight(node, foundPos, newCargo);
    return _insertNonFull(tree, getChild(node, foundPos), newKey, newCargo);  // Go into subtree left of split.
}

template <typename TFunctor, typename TKey, typename TCargo, unsigned B_FACTOR, typename TSpec, typename TKey2>
inline bool
_find(TFunctor & functor,
      DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> const & tree,
      TKey2 const & key)
{
    typedef DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> const TTree;
    typedef typename Value<TTree>::Type                               TNode;
    typedef typename TNode::TKeyTab                                   TKeyTab;
    typedef typename Iterator<TKeyTab const, Standard>::Type          TKeyIterator;
    typedef typename Size<TKeyTab>::Type                              TSize SEQAN_TYPEDEF_FOR_DEBUG;

    SEQAN_ASSERT_NEQ(tree.rootPtr, (TNode*) 0);

    TNode* currNodePtr = &_getRoot(tree);
    TKeyIterator it;

    while(true)
    {
        SEQAN_ASSERT_GT(length(_getKeyTable(*currNodePtr)) ,0u);
        // Binary search on the keys.
        it = std::upper_bound(begin(_getKeyTable(*currNodePtr), Standard()), end(_getKeyTable(*currNodePtr), Standard()), key);
        if (it != begin(_getKeyTable(*currNodePtr), Standard()))
        {
            --it;
            // Now either it compares equal to the value or it must be in the subtree right of the current it.
            functor(*currNodePtr, it - begin(_getKeyTable(*currNodePtr), Standard()));
            if (*it == key)
                return true;

            if (isLeaf(*currNodePtr))
                return false;
            SEQAN_ASSERT_LT(static_cast<TSize>((it - begin(_getKeyTable(*currNodePtr), Standard())) + 1), length(_getChildTable(*currNodePtr)));
            currNodePtr = &getChild(*currNodePtr, (it - begin(_getKeyTable(*currNodePtr), Standard())) + 1);
        }
        else
        {
            if (isLeaf(*currNodePtr))
                return false;
            currNodePtr = &getChild(*currNodePtr, 0);
        }
    }
}


// ----------------------------------------------------------------------------
// Function getKey()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR, typename TPos>
inline TKey const &
getKey(DpsTreeNode_<TKey, TCargo, B_FACTOR> const & node,
       TPos const & idx)
{
    SEQAN_ASSERT_GEQ(idx, static_cast<TPos>(0));
    SEQAN_ASSERT_LT(idx, static_cast<TPos>(length(node._keyTab)));

    return node._keyTab[idx];
}

// ----------------------------------------------------------------------------
// Function getCargo()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR, typename TPos>
inline TCargo const &
getCargo(DpsTreeNode_<TKey, TCargo, B_FACTOR> const & node,
         TPos const & idx)
{
    SEQAN_ASSERT_GEQ(idx, static_cast<TPos>(0));
    SEQAN_ASSERT_LT(idx, static_cast<TPos>(length(node._cargoTab)));

    return node._cargoTab[idx];
}

// ----------------------------------------------------------------------------
// Function getChild()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR, typename TPos>
inline DpsTreeNode_<TKey, TCargo, B_FACTOR> &
getChild(DpsTreeNode_<TKey, TCargo, B_FACTOR>  & node,
         TPos const & idx)
{
    SEQAN_ASSERT_GEQ(idx, static_cast<TPos>(0));
    SEQAN_ASSERT_LT(idx, static_cast<TPos>(length(node._childTab)));

    return *node._childTab[idx];
}

template <typename TKey, typename TCargo, unsigned B_FACTOR, typename TPos>
inline DpsTreeNode_<TKey, TCargo, B_FACTOR> const &
getChild(DpsTreeNode_<TKey, TCargo, B_FACTOR> const & node,
         TPos const & idx)
{
    SEQAN_ASSERT_GEQ(idx, static_cast<TPos>(0));
    SEQAN_ASSERT_LT(idx, static_cast<TPos>(length(node._childTab)));

    return *node._childTab[idx];
}

// ----------------------------------------------------------------------------
// Function isLeaf()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCargo, unsigned B_FACTOR>
inline bool
isLeaf(DpsTreeNode_<TValue, TCargo, B_FACTOR> const & node)
{
    return node._isLeaf;
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR, typename TSpec>
inline bool
empty(DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> const & tree)
{
    typedef DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> TTree SEQAN_TYPEDEF_FOR_DEBUG;
    typedef typename Value<TTree const>::Type                   TNode SEQAN_TYPEDEF_FOR_DEBUG;

    SEQAN_ASSERT_NEQ(tree.rootPtr, (TNode *) 0);

    return empty(tree.rootPtr->_keyTab) && empty(tree.rootPtr->_childTab);
}

// ----------------------------------------------------------------------------
// Function insert()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR, typename TSpec, typename TKey2, typename TCargo2>
inline bool
insert(DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> & tree,
       TKey2 const & newKey,
       TCargo2 const & newCargo)
{
    typedef DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> TTree;
    typedef typename Value<TTree>::Type TNode;

    if (empty(tree))
    {
        SEQAN_ASSERT_NEQ(tree.rootPtr, (TNode*) 0);

        appendValue(_getKeyTable(_getRoot(tree)), newKey);
        appendValue(_getCargoTable(_getRoot(tree)), newCargo);
        return true;
    }

    if (_isFull(_getRoot(tree)))  // Increase tree by one level.
    {
        TNode* tmp = _createNode(tree);  // tmp is no leaf, and initialized with zero.
        std::swap(tree.rootPtr, tmp);    // Swap the root and the new generated node.
        appendValue(_getChildTable(_getRoot(tree)), tmp);  // Store pointer to child in child tab.
        _splitChild(tree, getChild(_getRoot(tree), 0), _getRoot(tree), 0);  // Balance the tree.
    }
    return _insertNonFull(tree, _getRoot(tree), newKey, newCargo);
}

// ----------------------------------------------------------------------------
// Function prefixSum()
// ----------------------------------------------------------------------------

template <typename TKey, typename TCargo, unsigned B_FACTOR, typename TSpec, typename TKey2>
inline TCargo
prefixSum(DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> const & tree,
          TKey2 const & newKey)
{
    typedef DynamicPrefixSumTree<TKey, TCargo, B_FACTOR, TSpec> TTree;
    typedef typename Value<TTree>::Type TNode;

    PrefixSumFunctor_<TNode> func;

    _find(func, tree, newKey);
    return func.currSum;
}

}

#endif // EXTRAS_INCLUDE_SEQAN_MISC_MISC_DYNAMIC_PREFIX_SUM_H_
