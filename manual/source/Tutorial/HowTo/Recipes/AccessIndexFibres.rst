.. sidebar:: ToC

    .. contents::

.. _how-to-recipes-access-index-fibres:

Accessing Index Fibres
======================

Overview
--------

Basically each index consists of a set of tables, called fibres.
The set of available fibres of an index ``Index<TText, TSpec>`` depends on the index specialization ``TSpec``.

+------------------------+---------------------------------------------------------------------------------------------------------------------+
| :dox:`IndexEsa Fibres` | Description                                                                                                         |
+========================+=====================================================================================================================+
| EsaText                | The original text the index should be based on.                                                                     |
|                        | Can be either a sequence or a :dox:`StringSet`.                                                                     |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| EsaSA                  | The suffix array stores the begin positions of all suffixes in lexicographical order.                               |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| EsaLcp                 | The lcp table stores at position i the length of the longest common prefix between suffix with rank i and rank i+1. |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| EsaChildtab            | See :cite:`Abouelhoda2004`                                                                                          |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| EsaBwt                 | The Burrows-Wheeler table stores at position i the character left of the suffix with rank i.                        |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| EsaRawText             | Virtually concatenates all strings of the EsaText fibre.                                                            |
+------------------------+---------------------------------------------------------------------------------------------------------------------+

+------------------------+----------------------------------------------------------------------------------------------+
| :dox:`WOTDIndexFibres` | Description                                                                                  |
+========================+==============================================================================================+
| WotdText               | The original text the index should be based on.                                              |
+------------------------+----------------------------------------------------------------------------------------------+
| WotdSA                 | The suffix array stores the begin positions of all suffixes in lexicographical order.        |
+------------------------+----------------------------------------------------------------------------------------------+
| WotdDir                | :cite:`Giegerich2003`                                                                        |
+------------------------+----------------------------------------------------------------------------------------------+
| WotdRawText            | Virtually concatenates all strings of the WotdText fibre.                                    |
+------------------------+----------------------------------------------------------------------------------------------+

+-----------------------+--------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+
| :dox:`DfiIndexFibres` | Description                                                                                      | Type                                                                                         |
+=======================+==================================================================================================+==============================================================================================+
| DfiText               | The original text the index should be based on.                                                  | First template argument of the :dox:`Index`. Can be either a sequence or a :dox:`StringSet`. |
+-----------------------+--------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+
| DfiSA                 | The suffix array stores the begin positions of all suffixes in lexicographical order.            | String over the :dox:`SAValue` type of the index.                                            |
+-----------------------+--------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+
| DfiDir                | See :cite:`Giegerich2003`.                                                                       | String over the :dox:`Size` type of the index.                                               |
+-----------------------+--------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+
| DfiRawText            | Virtually concatenates all strings of the DfiText fibre.                                         | :dox:`ContainerConcept` over the alphabet of the text.                                       |
+-----------------------+--------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+

+-------------------------+------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+
| :dox:`QGramIndexFibres` | Description                                                                                                | Type                                                                                         |
+=========================+============================================================================================================+==============================================================================================+
| QGramText               | The original text the index should be based on.                                                            | First template argument of the :dox:`Index`. Can be either a sequence or a :dox:`StringSet`. |
+-------------------------+------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+
| QGramShape              | The q-gram :dox:`Shape`.                                                                                   | Specified by the first template argument of :dox:`IndexQGram`.                               |
+-------------------------+------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+
| QGramSA                 | The suffix array stores the begin positions of all suffixes in lexicographical order.                      | String over the :dox:`SAValue` type of the index.                                            |
+-------------------------+------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+
| QGramDir                | The directory maps q-gram hash values to start indices in the QGramSA fibre.                               | String over the :dox:`Index#Size` type of the index.                                         |
+-------------------------+------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+
| QGramCounts             | Stores numbers of occurrences per sequence of each q-gram in pairs (seqNo,count), count>0.                 | String over :dox:`Pair` of the :dox:`Index#Size` type of the index.                          |
+-------------------------+------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+
| QGramCountsDir          | The counts directory maps q-gram hash values to start indices in the QGramCounts fibre.                    | String over the :dox:`Index#Size` type of the index.                                         |
+-------------------------+------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+
| QGramBucketMap          | Used by the :dox:`OpenAddressingQGramIndex` index to store the hash value occupancy in the QGramDir fibre. | String over the :dox:`Value` type of the shape.                                              |
+-------------------------+------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+
| QGramRawText            | Virtually concatenates all strings of the QGramText fibre.                                                 | :dox:`ContainerConcept` over the alphabet of the text.                                       |
+-------------------------+------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------+

The first column in each table above contains the tags to select the corresponding fibre.
Most of these tags are aliases for the same tag, e.g. ``EsaSA``, ``QGramSA``, ... are aliases for ``FibreSA``.
If you write an algorithm that is generic in the type of index, you can use ``FibreText`` to specify the fibre that stores the index text.

Creation
--------

In most cases you don't need to create the fibres of an index by hand.
Most algorithms and data structures create them automatically, e.g. :dox:`Finder` or :dox:`VSTreeIterator`.
If you want to specify a certain index construction algorithm, have to recreate a fibre or manually access a fibre you can recreate or create on-demand a fibre by :dox:`Index#indexCreate` and :dox:`Index#indexRequire`.
If your algorithm should behave differently depending on the presence or absence of a fibre (and the fibre should then not be created), you can test for presence by :dox:`Index#indexSupplied`.

Access
------

The type of each fibre can be determined by the metafunction :dox:`Fibre`.
To access a fibre you can use the function :dox:`Index#getFibre` whose return type is the result of :dox:`Fibre`.
The second argument of both functions is a tag to select a specific fibre.
See the first column in the tables above.
One fibre in every index is the text to be indexed itself.
This fibre can be assigned during the construction.
For the ease of use, there exist shortcuts to access frequently used fibres:

+--------------------------------------------------------+---------------------------------------------------------+
| Shortcut                                               | Expands To ...                                          |
+========================================================+=========================================================+
| :dox:`IndexQGram#indexBucketMap indexBucketMap(index)` | :dox:`Index#getFibre getFibre(index, FibreBucketMap())` |
+--------------------------------------------------------+---------------------------------------------------------+
| :dox:`IndexEsa#indexBwt indexBwt(index)`               | :dox:`Index#getFibre getFibre(index, FibreBwt())`       |
+--------------------------------------------------------+---------------------------------------------------------+
| :dox:`IndexEsa#indexChildtab indexChildtab(index)`     | :dox:`Index#getFibre getFibre(index, FibreChildtab())`  |
+--------------------------------------------------------+---------------------------------------------------------+
| :dox:`IndexQGram#indexCounts indexCounts(index)`       | :dox:`Index#getFibre getFibre(index, FibreCounts())`    |
+--------------------------------------------------------+---------------------------------------------------------+
| :dox:`IndexQGram#indexCountsDir indexCountsDir(index)` | :dox:`Index#getFibre getFibre(index, FibreCountsDir())` |
+--------------------------------------------------------+---------------------------------------------------------+
| :dox:`IndexEsa#indexLcp indexLcp(index)`               | :dox:`Index#getFibre getFibre(index, FibreLcp())`       |
+--------------------------------------------------------+---------------------------------------------------------+
| :dox:`Index#indexRawText indexRawText(index)`          | :dox:`Index#getFibre getFibre(index, FibreRawText())`   |
+--------------------------------------------------------+---------------------------------------------------------+
| :dox:`IndexEsa#indexSA indexSA(index)`                 | :dox:`Index#getFibre getFibre(index, FibreSA())`        |
+--------------------------------------------------------+---------------------------------------------------------+
| :dox:`IndexQGram#indexShape indexShape(index)`         | :dox:`Index#getFibre getFibre(index, FibreShape())`     |
+--------------------------------------------------------+---------------------------------------------------------+
| :dox:`Index#indexText indexText(index)`                | :dox:`Index#getFibre getFibre(index, FibreText())`      |
+--------------------------------------------------------+---------------------------------------------------------+

and to access a single value:

+----------------------------------------------+---------------------------------------------------------+
| Shortcut                                     | Expands To ...                                          |
+==============================================+=========================================================+
| :dox:`IndexEsa#bwtAt bwtAt(pos, index)`      | :dox:`IndexEsa#indexBwt indexBwt(index)[pos]`           |
+----------------------------------------------+---------------------------------------------------------+
| :dox:`IndexEsa#childAt childAt(pos, index)`  | :dox:`IndexEsa#indexChildtab indexChildtab(index)[pos]` |
+----------------------------------------------+---------------------------------------------------------+
| :dox:`IndexQGram#dirAt dirAt(pos, index)`    | :dox:`IndexQGram#indexDir indexDir(index)[pos]`         |
+----------------------------------------------+---------------------------------------------------------+
| :dox:`IndexEsa#lcpAt lcpAt(pos, index)`      | :dox:`IndexEsa#indexLcp indexLcp(index)[pos]`           |
+----------------------------------------------+---------------------------------------------------------+
| :dox:`Index#rawtextAt rawtextAt(pos, index)` | :dox:`Index#indexRawText indexRawText(index)[pos]`      |
+----------------------------------------------+---------------------------------------------------------+
| :dox:`IndexEsa#saAt saAt(pos, index)`        | :dox:`IndexEsa#indexSA indexSA(index)[pos]`             |
+----------------------------------------------+---------------------------------------------------------+
| :dox:`IndexEsa#textAt textAt(pos, index)`    | :dox:`Index#indexText indexText(index)[pos]`            |
+----------------------------------------------+---------------------------------------------------------+

Please note that :dox:`IndexEsa#textAt` can also be used if the index text is a :dox:`StringSet`.
``pos`` can then be a :dox:`SAValue`.
