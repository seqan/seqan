Tutorial FMIndex
----------------

TOC

| ``Learning Objective ::``
| `` This tutorial will teach you how to configure and use SeqAn's FMIndex.``
| ``Difficulty ::``
| `` Average``
| ``Duration ::``
| `` 1h``
| ``Prerequisites ::``
| `` ``\ ```Tutorial/Indices`` <Tutorial/Indices>`__

The FM Index
~~~~~~~~~~~~

The FM index is a Burrows-Wheeler transform (BWT) based index developed
by `Ferragina and Manzini <Bibliography#FerraginaManzini2001>`__.

The index consists of two major fibres (seqan:Class.Fibre) which are a
LF table (storing all necessary information for the LF mapping) and a
compressed suffix array (used to retrieve positions of hits of a pattern
in a text). Per default the types of the involved fibres a chosen such
that text of arbitrary length can be processed. However, if your text
length is smaller than 4.29 billion characters you might want to adjust
some data types as shown next.

Reducing the memory consumption
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several ways to reduce the memory consumption.

#. Reduce the default sample rate of the compressed suffix array and

``2. Change the type of the compressed suffix array such that each entry requires less memory.``

You can make use of the first option by specifying the sample rate in
the constructor of the FM index (seqan:Spec.FMIndex) like this:

::

    #c++
    String<Dna> text = "...";

    Index<String<Dna>>, FMIndex> fmIndex(text, 1000);
    ...

Here the sample rate was increased to 1000, meaning that only 1 in 1000
suffix array entries of the complete suffix array is stored in the
compressed suffix array. The major drawback is an increase in running
time if one is interested in locating hits of queries in a text.

The second option has no drawback concerning running time but one has to
make sure that the text to index does not exceed 4.29 billion
characters. The critical observation is that each suffix array entry
consumes 64 bit of memory per default where 32 bit would be sufficient
if the text size is appropriate. In order to change the size type of the
suffix array entry we simply have to overload the metafunction
``SAValue``.

::

    #c++

    template<>
    struct SAValue<String<Dna> >
    {
        typedef unsigned Type;
    }

If your text is a seqan:Class.StringSet than ``SAValue`` will return a
seqan:Class.Pair that can be overloaded in the same way.

::

    #c++
    template<>
    struct SAValue<StringSet<String<Dna> > >
    {
        typedef Pair<unsigned, unsigned> Type;
    }

The first type of the pair is used as the type for the index of a string
in the string set. So if you only have a few strings you could save even
more memory like this:

::

    #c++
    template<>
    struct SAValue<StringSet<String<Dna> > >
    {
        typedef Pair<unsigned char, unsigned> Type;
    }

.. raw:: mediawiki

   {{TracNotice|{{PAGENAME}}}}
