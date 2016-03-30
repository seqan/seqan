.. _tutorial-datastructures-sequences:

Sequences
=========

..  toctree::
    :hidden:
    :titlesonly:

    StringsAndSegments
    StringSets
    Alphabets

Tasks, such as computing an alignment, searching for patterns online or indexing a genome or protein database are required in many bioinformatics applications.
For these tasks, SeqAn provides fundamental data structures to work efficiently with biological sequences.
SeqAn implements special alphabets, such as DNA, RNA, AminoAcid, Iupac and more.
The alphabets available in SeqAn can be reviewed :ref:`here <tutorial-datastructures-sequences-alphabets>`.
You will also find some more information about using the alphabet types.

Besides the alphabets SeqAn also implements a string class.
SeqAn strings are generic containers in which characters of any alphabet are stored continuously in memory.
The default string class implementation is equivalent to the STL `vector <http://en.cppreference.com/w/cpp/container/vector>`_ class.
However, the memory mangement of the SeqAn string class is optimized for working with SeqAn's alphabets.
Apart of the default string class implementation SeqAn provides many useful specializations of this class, which are very useful in the bioinformatics context.
The tutorial about :ref:`tutorial-datastructures-sequences-strings-and-segments` gives you a more detailed overview over the string class and it's functionality.

Another generic container data structure used very often is the string set.
The string set is a special container to represent a collection of strings, which in addition provides many helpful functions to make the work even easier and more efficient.
The :ref:`StringSet tutorial <tutorial-datastructures-sequences-string-sets>` introduces you to this class and how to work with it.
