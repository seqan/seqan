TOC

Indices
-------

A substring index is a datatype which allows one to seek efficiently for all
occurrences of a pattern in a string or a set of strings. Substring
indices are very efficient for the exact string matching problem, i.e.
finding all exact occurrences of a pattern in a text or a text
collection. Instead of searching through the text in O(n) like
online-search algorithms do, a substring index looks up the pattern in
sublinear time o(n). Substring indices are full-text indices, i.e. they
handle all substrings of a text in contrast to inverted files or
signature files, which need word delimiters. SeqAn contains data
structures to create, hold and use substring indices. Based on a unified
concept, SeqAn offers the following concrete implementations defined as
specializations of seqan:Class.Index:

+-------------------------------------------------+-------------------------------------------------------------------------------------+
| **Specialization**                              | **Description**                                                                     |
+=================================================+=====================================================================================+
| seqan:Spec.IndexEsa                             | Abouelhoda et al., 2004]])                                                          |
+-------------------------------------------------+-------------------------------------------------------------------------------------+
| seqan:Spec.IndexWotd                            | Giegerich et al., 2003]])                                                           |
+-------------------------------------------------+-------------------------------------------------------------------------------------+
| seqan:Spec.IndexDfi                             | Weese, Schulz, 2008]])                                                              |
+-------------------------------------------------+-------------------------------------------------------------------------------------+
| seqan:Spec.IndexQGram                           | q-gram index                                                                        |
+-------------------------------------------------+-------------------------------------------------------------------------------------+
| [seqan:"Spec.Pizza & Chili Index" PizzaChili]   | An adapter for the `Pizza & Chili <http://pizzachili.dcc.uchile.cl/>`__ index API   |
+-------------------------------------------------+-------------------------------------------------------------------------------------+

Suffix Tree Interface
~~~~~~~~~~~~~~~~~~~~~

The unified concept allows the first three indices (seqan:Spec.IndexEsa,
seqan:Spec.IndexWotd, seqan:Spec.IndexDfi) to be accessed just like a
`suffix tree <Tutorial/Indices/SuffixTree>`__ independently of its
concrete implementation. To access this (virtual) suffix tree SeqAn
offers various [seqan:"Spec.VSTree Iterator" iterators].

Depth-First Search
^^^^^^^^^^^^^^^^^^

In SeqAn a suffix tree (see definition
`here <Tutorial/Indices/SuffixTree>`__) can be accessed with special
suffix tree iterators, which differ in the way the tree nodes are
traversed. For many sequence algorithms it is neccessary to do a full
depth-first search (dfs) over all suffix tree nodes beginning either in
the root (preorder dfs) or in a leaf node (postorder dfs). A preorder
traversal (Fig.1) halts in a node when visiting it for the first time
whereas a postorder traversal (Fig.2) halts when visiting a node for the
last time. The following two figures give an example in which order the
tree nodes are visited.

+---------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------+
| `Image(source:trunk/docs/img/streePreorder.png, 300px) <Image(source:trunk/docs/img/streePreorder.png, 300px)>`__   | `Image(source:trunk/docs/img/streePostorder.png, 300px) <Image(source:trunk/docs/img/streePostorder.png, 300px)>`__   |
+=====================================================================================================================+=======================================================================================================================+
| **Figure 1:** Preorder DFS                                                                                          | **Figure 2:** Postorder DFS                                                                                           |
+---------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------+

There are currently 2 iterators in SeqAn supporting a DFS search:

+----------------------------------------+----------------+-----------------+
| **Iterator**                           | **Preorder**   | **Postorder**   |
+========================================+================+=================+
| seqan:"Spec.BottomUp Iterator"         | -              | +               |
+----------------------------------------+----------------+-----------------+
| seqan:"Spec.TopDownHistory Iterator"   | +              | +               |
+----------------------------------------+----------------+-----------------+

If solely a postorder traversal is needed the seqan:"Spec.BottomUp
Iterator" should be preferred as it is more memory efficient. Please
note that the BottomUp Iterator is only applicable to
seqan:Spec.IndexEsa indices.

Example
^^^^^^^

We want to construct the suffix tree of the string "abracadabra" and
output the substrings represented by tree nodes in preorder dfs. All
index related data structures are defined in ``seqan/index.h``, so we
have to include this header (``seqan/sequence.h`` is included
implicitly).
.. includefrags:: demos/tutorial/index/index_preorder.cpp
   :fragment: includes

Next we create the string "abracadabra" and an index specialized with
the type of this string. The string can be given to the index
constructor.

.. includefrags:: demos/tutorial/index/index_preorder.cpp
   :fragment: initialization

The seqan:Metafunction.Iterator metafunction expects two arguments, the
type of the container to be iterated and a specialization tag, see the
seqan:"Spec.VSTree Iterator" hierarchy. In this example we chose a
seqan:"Spec.TopDownHistory Iterator" whose signature in the second
template argument is ``TopDown< ParentLinks<Preorder> >``. Each suffix
tree iterator constructor expects at least the index object and optional
problem-dependent parameters. As all DFS suffix tree iterators implement
the seqan:Concept.Iterator concept, they can be used via
seqan:Function.goNext, seqan:Function.atEnd, etc. The string that
represents the node the iterator points to is returned by
seqan:Function.representative.
.. includefrags:: demos/tutorial/index/index_preorder.cpp
   :fragment: iteration

Program output:

::

    #html
    <pre class="wiki" style="background-color:black;color:lightgray">

    a
    abra
    abracadabra
    acadabra
    adabra
    bra
    bracadabra
    cadabra
    dabra
    ra
    racadabra

.. raw:: html

   </pre>

**Note:** A relaxed suffix tree (see
`definition <Tutorial/Indices/SuffixTree>`__) is a suffix tree after
removing the $ characters and empty edges. For some bottom-up algorithms
it would be better not to remove empty edges and to have a one-to-one
relationship between leaves and suffices. In that cases you can use the
tags PreorderEmptyEdges or PostorderEmptyEdges instead of Preorder or
Postorder or EmptyEdges for the TopDown Iterator.

Assignments
^^^^^^^^^^^

| *``Task``
``1``*\ `` :: Write a program that constructs an index of the seqan:Class.StringSet "tobeornottobe", "thebeeonthecomb", "beingjohnmalkovich" and outputs the strings corresponding to suffix tree nodes in postorder DFS.``
| *``Difficulty``*\ `` :: 2``
| *``Solution``*\ `` :: can be found ``\ ```here`` <Tutorial/Indices/Assignment1>`__

| *``Task``
``2``*\ `` :: Write a program that outputs all maximal unique matches (MUMs) between "CDFGHC" and "CDEFGAHC".``
| *``Difficulty``*\ `` :: 2``
| *``Solution``*\ `` ::  can be found ``\ ```here`` <Tutorial/Indices/Assignment2>`__

Top-Down Iteration
^^^^^^^^^^^^^^^^^^

For index based pattern search or algorithms traversing only the upper
parts of the suffix tree the seqan:"Spec.TopDown Iterator" or
seqan:"Spec.TopDownHistory Iterator" is the best solution. Both provide
the functions seqan:Function.goDown and seqan:Function.goRight to go
down to the first child node or go to the next sibling. The
seqan:"Spec.TopDownHistory Iterator" additionally provides
seqan:Function.goUp to go back to the parent node. The child nodes in
seqan:Spec.IndexEsa indices are lexicographically sorted from first to
last. For seqan:Spec.IndexWotd and seqan:Spec.IndexDfi indices this
holds for all children except the first.

Example
^^^^^^^

In the next example we want to use the seqan:"Spec.TopDown Iterator" to
efficiently search a text for exact matches of a pattern. We therefore
want to use seqan:Function.goDown which has an overload to go down an
edge beginning with a specific character. First we create an index of
the text "How many wood would a woodchuck chuck."
.. includefrags:: demos/tutorial/index/index_search.cpp
   :fragment: initialization

The main search can then be implemented as follows. The algorithm
descends the suffix tree along edges beginning with the corresponding
pattern character. In each step the unseen edge characters have to be
verified.
.. includefrags:: demos/tutorial/index/index_search.cpp
   :fragment: iteration

If all pattern characters could successfully be compared we end in the
topmost node pattern is a prefix of. Thus, the suffixes represented by
this node are the occurrences of our pattern.
.. includefrags:: demos/tutorial/index/index_search.cpp
   :fragment: output

Program output:

::

    #html
    <pre class="wiki" style="background-color:black;color:lightgray">
    w
    wo
    wood
    9
    22

.. raw:: html

   </pre>

Alternatively, we could have used seqan:Function.goDown to go down the
path of a pattern instead single characters:
.. includefrags:: demos/tutorial/index/index_search2.cpp
   :fragment: output

::

    #html
    <pre class="wiki" style="background-color:black;color:lightgray">
    9
    22

.. raw:: html

   </pre>

Assignments
^^^^^^^^^^^

| *``Task``
``3``*\ `` ::  Write a program that iterates over all nodes of the suffix tree of the string "tobeornottobe" in preorder DFS. Use seqan:Function.goDown, seqan:Function.goRight and seqan:Function.goUp to iterate instead of seqan:Function.goNext or the operator++. Output the representatives.``
| *``Difficulty``*\ `` :: 4``
| *``Solution``*\ `` :: can be found ``\ ```here`` <Tutorial/Indices/Assignment3>`__

| *``Task``
``4``*\ `` ::  Modify the program to efficiently skip nodes with representatives longer than 3. Move the whole program into a template function whose argument specifies the index type and call this function twice, once for the seqan:Spec.IndexEsa and once for the seqan:Spec.IndexWotd index.``
| *``Difficulty``*\ `` :: 5``
| *``Solution``*\ `` ::  can be found ``\ ```here`` <Tutorial/Indices/Assignment4>`__

Access Suffix Tree Nodes
^^^^^^^^^^^^^^^^^^^^^^^^

In the previous subsection we have seen how to walk through a suffix
tree. We now want to know what can be done with a suffix tree iterator.
As all iterators are specializations of the general VSTree Iterator
class, they inherit all of its functions. There are various functions to
access the node the iterator points at, so we concentrate on the most
important ones.

+-----------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| **Function**                                                                      | **Description**                                                                                                                              |
+===================================================================================+==============================================================================================================================================+
| seqan:Function.representative                                                     | returns the substring that represents the current node, i.e. the concatenation of substrings on the path from the root to the current node   |
+-----------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| seqan:Function.getOccurrence                                                      | returns a position where the representative occurs in the text                                                                               |
+-----------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| seqan:Function.getOccurrences                                                     | returns a string of all positions where the representative occurs in the text                                                                |
+-----------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| seqan:Function.isRightTerminal                                                    | suffix tree]] figures)                                                                                                                       |
+-----------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| `isLeaf <http://www.seqan.de/dddoc/html_devel/FUNCTION_Index_23is_Leaf.html>`__   | tests if the current node is a tree leaf                                                                                                     |
+-----------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| seqan:Function.parentEdgeLabel                                                    | returns the substring that represents the edge from the current node to its parent (only TopDownHistory Iterator)                            |
+-----------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+

**Note:** There is a difference between the functions isLeaf and
isRightTerminal. In a relaxed suffix tree (see
`definition <Tutorial/Indices/SuffixTree>`__) a leaf is always a suffix,
but not vice versa, as there can be internal nodes a suffix ends in. For
them isLeaf returns false and isRightTerminal returns true.

Property Maps
^^^^^^^^^^^^^

Some algorithms require to store auxiliary information (e.g. weights,
scores) to the nodes of a suffix tree. To attain this goal SeqAn
provides so-called property maps, simple Strings of a property type.
Before storing a property value, these strings must first be resized
with seqan:Function.resizeVertexMap. The property value can then be
assigned or retrieved via seqan:Function.assignProperty or
seqan:Function.getProperty, seqan:Function.property. It is recommended
to call seqan:Function.resizeVertexMap prior to every call of
seqan:Function.assignProperty to ensure that the property map has
sufficient size. The following example iterates over all nodes in
preorder dfs and recursively assigns the node depth to each node. First
we create a seqan:Class.String of ``int`` to store the node depth for
each suffix tree node.
.. includefrags:: demos/tutorial/index/index_property_maps.cpp
   :fragment: initialization
The main loop iterates over all nodes in preorder DFS, i.e. parents are
visited prior children. The node depth for the root node is 0 and for
all other nodes it is the parent node depth increased by 1. The
functions seqan:Function.assignProperty, seqan:Function.getProperty and
seqan:Function.property must be called with a
seqan:Metafunction.VertexDescriptor. The vertex descriptor of the
iterator node is returned by seqan:Function.value and the descriptor of
the parent node is returned by seqan:Function.nodeUp.
.. includefrags:: demos/tutorial/index/index_property_maps.cpp
   :fragment: iteration
At the end we again iterate over all nodes and output the calculated
node depth.
.. includefrags:: demos/tutorial/index/index_property_maps.cpp
   :fragment: output
Program output:

::

    #html
    <pre class="wiki" style="background-color:black;color:lightgray">
    0
    1       a
    2       abra
    3       abracadabra
    2       acadabra
    2       adabra
    1       bra
    2       bracadabra
    1       cadabra
    1       dabra
    1       ra
    2       racadabra

.. raw:: html

   </pre>

*``Hint``*\ `` :: In SeqAn there is already a function seqan:Function.nodeDepth defined to return the node depth.``

Additional iterators
^^^^^^^^^^^^^^^^^^^^

By now, we know the following iterators (n=text size, σ=alphabet size,
d=tree depth):

+----------------------------------------+------------------------------------------+-------------+---------------------+
| **Iterator specialization**            | **Description**                          | **Space**   | **Index tables**    |
+========================================+==========================================+=============+=====================+
| seqan:"Spec.BottomUp Iterator"         | postorder dfs                            | O(d)        | SA, LCP             |
+----------------------------------------+------------------------------------------+-------------+---------------------+
| seqan:"Spec.TopDown Iterator"          | can go down and go right                 | O(1)        | SA, Lcp, Childtab   |
+----------------------------------------+------------------------------------------+-------------+---------------------+
| seqan:"Spec.TopDownHistory Iterator"   | can also go up, preorder/postorder dfs   | O(d)        | SA, Lcp, Childtab   |
+----------------------------------------+------------------------------------------+-------------+---------------------+

Besides the iterators described above, there are some
application-specific iterators in SeqAn:

+---------------------------------------------+-----------------------------------------------------------+-------------+--------------------------+
| **Iterator specialization**                 | **Description**                                           | **Space**   | **Index tables**         |
+=============================================+===========================================================+=============+==========================+
| seqan:"Spec.MaxRepeats Iterator"            | maximal repeats                                           | O(n)        | SA, Lcp, Bwt             |
+---------------------------------------------+-----------------------------------------------------------+-------------+--------------------------+
| seqan:"Spec.SuperMaxRepeats Iterator"       | supermaximal repeats                                      | O(d+σ)      | SA, Lcp, Childtab, Bwt   |
+---------------------------------------------+-----------------------------------------------------------+-------------+--------------------------+
| seqan:"Spec.SuperMaxRepeatsFast Iterator"   | supermaximal repeats (optimized for enh. suffix arrays)   | O(σ)        | SA, Lcp, Bwt             |
+---------------------------------------------+-----------------------------------------------------------+-------------+--------------------------+
| seqan:"Spec.MUMs Iterator"                  | maximal unique matches                                    | O(d)        | SA, Lcp, Bwt             |
+---------------------------------------------+-----------------------------------------------------------+-------------+--------------------------+
| seqan:"Spec.MultiMEMs Iterator"             | multiple maximal exact matches (w.i.p.)                   | O(n)        | SA, Lcp, Bwt             |
+---------------------------------------------+-----------------------------------------------------------+-------------+--------------------------+

Given a string s a repeat is a substring r that occurs at 2 different
positions i and j in s. The repeat can also be identified by the triple
(i,j,\|r\|). A maximal repeat is a repeat that cannot be extended to the
left or to the right, i.e. s[i-1]≠s[j-1] and s[i+\|r\|]≠s[j+\|r\|]. A
supermaximal repeat r is a maximal repeat that is not part of another
repeat. Given a set of strings s1, ..., sm a MultiMEM (multiple maximal
exact match) is a substring r that occurs in each sequence si at least
once and cannot be extended to the left or to the right. A MUM (maximal
unique match) is a MultiMEM that occurs exactly once in each sequence.
The following examples demonstrate the usage of these iterators:

+---------------------------------------+
| **Example**                           |
+=======================================+
| seqan:"Demo.Maximal Unique Matches"   |
+---------------------------------------+
| seqan:"Demo.Supermaximal Repeats"     |
+---------------------------------------+
| seqan:"Demo.Maximal Repeats"          |
+---------------------------------------+

q-gram Index
~~~~~~~~~~~~

A q-gram index can be used to efficiently retrieve all occurrences of a
certain q-gram in the text. It consists of various tables, called fibres
(see `HowTo <HowTo/AccessIndexFibres>`__), to retrieve q-gram positions,
q-gram counts, etc. However, it has no support for suffix tree
iterators. A q-gram index must be specialized with a seqan:Class.Shape
type. A seqan:Class.Shape defines q, the number of characters in a
q-gram and possibly gaps between these characters. There are different
specializations of seqan:Class.Shape available:

+-----------------------------+--------------------+----------------------+
| **Specialization**          | **Modifiable\***   | **Number of Gaps**   |
+=============================+====================+======================+
| seqan:Spec.UngappedShape    | -                  | 0                    |
+-----------------------------+--------------------+----------------------+
| seqan:Spec.SimpleShape      | +                  | 0                    |
+-----------------------------+--------------------+----------------------+
| seqan:Spec.OneGappedShape   | +                  | 0/1                  |
+-----------------------------+--------------------+----------------------+
| seqan:Spec.GappedShape      | -                  | any                  |
+-----------------------------+--------------------+----------------------+
| seqan:Spec.GenericShape     | +                  | any                  |
+-----------------------------+--------------------+----------------------+

-  - *fixed at compile time*, + *can be changed at runtime*

Each shape evaluates a gapped or ungapped sequence of q characters to a
hash value by the Functions seqan:Function.hash,
seqan:Function.hashNext, etc. For example, the shape 1101 represents a
3-gram with one gap of length 1. This shape overlayed with the
seqan:Spec.Dna text "GATTACA" at the third position corresponds to
"TT-C". The function seqan:Function.hash converts this 3-gram into
61=((\ **3**\ \*4+\ **3**)\*4+\ **1**. 4 is the alphabet size in this
example (see seqan:Metafunction.ValueSize).

The q-gram index offers different function to search or count
occurrences of q-grams in an indexed text, see
seqan:Function.getOccurrences, seqan:Function.countOccurrences. A q-gram
index over a seqan:Class.StringSet stores occurrence positions in the
same way as the ESA index and in the same fibre (Fibre\_SA). If only the
number of q-grams per sequence are needed the QGram\_Counts and
QGram\_CountsDir fibres can be used. They store pairs
``(seqNo, count)``, ``count``>0, for each q-gram that occurs ``counts``
times in sequence number ``seqNo``.

To efficiently retrieve all occurrence positions or all pairs
``(seqNo, count)`` for a given q-gram, these positions or pairs are
stored in contiguous blocks (in QGram\_SA, QGram\_Counts fibres), called
buckets. The begin position of bucket i is stored in directory fibres
(QGram\_Dir, QGram\_CountsDir) at position i, the end position is the
begin positions of the bucket i+1. The default implementation of the
seqan:Spec.IndexQGram index maps q-gram hash values 1-to-1 to bucket
numbers. For large q or large alphabets the seqan:Spec.OpenAddressing
index can be more appropriate as its directories are additionally bound
by the text length. This is realized by a non-trivial mapping from
q-gram hashes to bucket numbers that requires an additional fibre
(QGram\_BucketMap).

For more details on q-gram index fibres see the
`HowTo <HowTo/AccessIndexFibres>`__ or seqan:"Tag.QGram Index Fibres".

Example
^^^^^^^

We want to construct the q-gram index of the string "CATGATTACATA" and
output the occurrences of the ungapped 3-gram "CAT". As 3 is fixed at
compile-time and the shape has no gaps we can use a
seqan:Spec.UngappedShape which is the first template argument of
seqan:Spec.IndexQGram, the second template argument of
seqan:Class.Index. Next we create the string "CATGATTACATA" and
specialize the first index template argument with the type of this
string. The string can be given to the index constructor.
.. includefrags:: demos/tutorial/index/index_qgram.cpp
   :fragment: initialization

To get all occurrences of a q-gram, we first have to hash it with a
shape of the same type as the index shape (we can even use the index
shape returned by seqan:Function.indexShape). The hash value returned by
seqan:Function.hash or seqan:Function.hashNext is also stored in the
shape and is used by the function seqan:Function.getOccurrences to
retrieve all occurrences of our 3-gram.
.. includefrags:: demos/tutorial/index/index_qgram.cpp
   :fragment: output

Program output:

::

    #html
    <pre class="wiki" style="background-color:black;color:lightgray">
    0
    8

.. raw:: html

   </pre>

Assignments
^^^^^^^^^^^

| *``Task``
``5``*\ `` ::  Write a program that outputs all occurrences of the gapped q-gram "AT-A" in "CATGATTACATA".``
| *``Difficulty``*\ `` :: 3``
| *``Solution``*\ `` ::  can be found ``\ ```here`` <Tutorial/Indices/Assignment5>`__

| *``Task``
``6``*\ `` :: Create and output a matrix M where M(i,j) is the number of common ungapped 5-grams between sequence i and sequence j for 3 random seqan:Spec.Dna sequences, each not longer than 200 characters. Optional: Run the matrix calculation twice, once for an seqan:Spec.IndexQGram and once for an seqan:Spec.OpenAddressing index and output the directory sizes (QGram_Dir, QGram_CountsDir fibre).``
| *``Difficulty``*\ `` :: 5``
| *``Hint``*\ `` :: A common g-gram that occurs a times in one and b times in the other sequence counts for min(a,b).``
| *``Solution``*\ `` ::  can be found ``\ ```here`` <Tutorial/Indices/Assignment6>`__

Handling Multiple Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The previous sections briefly described how an index of a set of strings
can be instantiated. Instead of creating an seqan:Class.Index of a
seqan:Class.String you create one of a seqan:Class.StringSet. A
character position of this string set can be one of the following:

#. A local position (default), i.e. seqan:Class.Pair (seqNo, seqOfs)
   where seqNo identifies the string within the stringset and the seqOfs
   identifies the position within this string.

``2. A global position, i.e. single integer value between 0 and the sum of string lengths minus 1 (global position). This integer is the position in the gapless concatenation of all strings in the seqan:Class.StringSet to a single string.``
``The meta-function seqan:Metafunction.SAValue determines, which position type (local or global) will be used for internal index tables (suffix array, q-gram array) and what type of position is returned by functions like seqan:Function.getOccurrence or seqan:Function.position of a seqan:Class.Finder. ``
``seqan:Metafunction.SAValue returns a seqan:Class.Pair = local position by default, but could be specialized to return an integer type = global position for some applications.``
``If you want to write algorithms for both variants you should use the functions seqan:Function.posLocalize, seqan:Function.posGlobalize, seqan:Function.getSeqNo and seqan:Function.getSeqOffset.``

Submit a comment
^^^^^^^^^^^^^^^^

If you found a mistake, or have suggestions about an improvement of this
page press:
[/newticket?component=Documentation&description=Tutorial+Enhancement+for+page+http://trac.seqan.de/wiki/Tutorial/Indices&type=enhancement
submit your comment]

.. raw:: mediawiki

   {{TracNotice|{{PAGENAME}}}}
