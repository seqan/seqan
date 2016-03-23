.. _tutorial-datastructures-sequences:

Sequences
=========

..  toctree::
    :hidden:
    :titlesonly:

    Alphabets
    StringsAndSegments
    StringSets
    JournaledSet

Tasks, such as computing an alignment, searching for patterns online or indexing a genome or protein database are required in many bioinformatics applications.
For these tasks, SeqAn provides fundamental data structures to work efficiently with biological sequences.
SeqAn implements special alphabets, such as DNA, RNA, AminoAcid, Iupac and more.
The alphabets available in SeqAn can be reviewed :ref:`here <tutorial-datastructures-sequences-alphabets>`.
You will also find some more information about usage of alphabet types.

Besides the alphabets SeqAn also implements a string class.
SeqAn Strings are generic containers in which characters of any alphabet are stored continuously in memory.
The default String class implementation is equivalent to the STL `vector <http://en.cppreference.com/w/cpp/container/vector>`_ class.
However, the memory mangement of the SeqAn's String class is optimized for working with SeqAn's alphabets.
Also are many more specializations of the String class, that are needed in bioinformatic applications frequently.
The tutorial about :ref:`tutorial-datastructures-sequences-strings-and-segments` gives you a more detailed overview over the String class and it's functionality.

Another generic container data structure used very often is the StringSet.
The StringSet is a special container to represent a collection of Strings which provides many helpful functions to make the work even easier and more efficient.
The :ref:`StringSet tutorial <tutorial-datastructures-sequences-string-sets>` introduces you to this class and its functionalities.
