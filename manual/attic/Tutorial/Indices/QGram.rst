`Index\_QGram <http://www.seqan.de/dddoc/html/Spec_Index_QGram.html>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

--------------

''' I've instantiated an Index and delivered a StringSet. Is the q-gram
index thereby already created? '''

::

    Index< StringSet<String<Dna>>, Index_QGram< FixedShape<4> > > TQGramIndex;
    TQGramIndex myIndex(myReads);

No, not yet. Although the text to be indexed is given, the index tables
(called fibres) will be created on demand. If you explicitly want to
access the list of q-gram occurrences (SuffixArray) or the directory
mapping q-grams to the first occurrence entry (Directory), you have to
invoke the following function that creates both tables (if they haven't
been created already):

::

    indexRequire(myIndex, QGram_SADir());

--------------

''' How can I create a q-gram index with q given by a variable? '''

`FixedShape <http://www.seqan.de/dddoc/html/Spec_FixedShape.html>`__
(ungapped q-gram) or
`FixedGappedShape <http://www.seqan.de/dddoc/html/Spec_FixedGappedShape.html>`__
(gapped q-gram) can only be used if q (or the gapped shape) is fixed at
compile time. If q or the q-gram shape is variable and not known at
compile time, you have to specialize Index\_QGram using a
`SimpleShape <http://www.seqan.de/dddoc/html/Spec_SimpleShape.html>`__
(ungapped) or
`GappedShape <http://www.seqan.de/dddoc/html/Spec_GappedShape.html>`__
(gapped). Directly after the instantiation of the index, you should set
the q or the gapped shaped with resize(indexShape(myIndex), q) or
`shapeToString <http://www.seqan.de/dddoc/html/Function_shapeToString.html>`__\ (indexShape(myIndex),
bitString):

::

    int q = 5;
    Index< StringSet<String<Dna>>, Index_QGram< SimpleShape > > TQGramIndex;
    TQGramIndex myIndex(myReads);
    resize(indexShape(myIndex), q);

--------------

**Is there a function to save an index to disk to load it again without
recreating?**

Storing and loading an index can be done with:

::

    const char *fileName = "/home/user/myindex";
    save(index, fileName);

or

::

    const char *fileName = "/home/user/myindex";
    open(index, fileName);

If you have built your q-gram index with variable shapes (i.e.
`SimpleShape <http://www.seqan.de/dddoc/html/Spec_SimpleShape.html>`__
(ungapped) or
`GappedShape <http://www.seqan.de/dddoc/html/Spec_GappedShape.html>`__
(gapped)), you have to mind that q or the shape is not stored or loaded.
This must be done manually directly before or after loading with resize
(SimpleShape) oder stringToShape (GappedShape).

**Hint:** A newly instantiated index is initially empty. If you assign a
text to be indexed, solely the text fibre is set (fibre=a certain index
table, see
`Fibres <http://www.seqan.de/dddoc/html/Tag_QGram%20Index%20Fibres.html>`__).
All other fibres are empty and created on demand. Normally, a full
created index should be saved to disk. Therefore, you have to create the
required fibres explicitly by hand:

::

    const char *fileName = "/home/user/myindex";
    indexRequire(index, QGram_SADir());
    save(index, fileName);

For the `ESA index <Tutorial/Indices/ESA>`__ you would do:

::

    const char *fileName = "/home/user/myindex";
    indexRequire(index, ESA_SA());
    indexRequire(index, ESA_LCP());
    indexRequire(index, ESA_ChildTab());  // for TopDown iterators
    indexRequire(index, ESA_BWT());       // for (Super-)MaxRepeats iterators
    save(index, fileName);

--------------

''' I've created a q-gram index of a set of sequences (using a
StringSet). How can I determine all the sequences and positions a
certain q-gram occurs at? '''

Suppose, an index of an ungapped 5-gram was created. The sequence
numbers and positions within of occurrences of a certain 5-gram can be
dumped as follows:

::

    String<Dna> kmer = "agaag";
    unsigned hashValue = hash(indexShape(index), begin(kmer));
    for(unsigned i = dirAt(hashValue, index); i < dirAt(hashValue +1, index); ++i)
      ::std::cout << getSeqNo(saAt(i,index)) << ":" << getSeqOffset(saAt(i,index)) << ::std::endl;

The **suffix array** contains all occurrences of q-grams, s.t. the
occurrences of a single q-gram are a stored in a contiguous block. The
**directory** stores for any possible q-gram the index of the first
occurrence in the suffix array (see `QGram Index
Fibres <http://www.seqan.de/dddoc/html/Tag_QGram%20Index%20Fibres.html>`__).
If the hash value of a g-gram is determined, the q-gram occurrences are
stored from position ``dirAt(hash,index)`` to ``dirAt(hash+1,index)-1``
in the suffix array. An occurrence is stored either as a single integer
(index was built over a single string) or as a Pair (index of a
StringSet). The sequence number can be obtained with ``getSeqNo`` and
the position within the sequence with ``getSeqOffset``.

Alternatively, the Finder class can be used to determine all occurrences
of a certain q-gram. See
`Example <http://www.seqan.de/dddoc/html/Demo_Index%20Finder.html>`__.

--------------

.. raw:: mediawiki

   {{TracNotice|{{PAGENAME}}}}
