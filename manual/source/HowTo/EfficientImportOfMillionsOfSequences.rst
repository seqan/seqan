.. sidebar:: ToC

   .. contents::


.. _how-to-efficiently-import-millions-of-sequences:

Efficiently Importing Millions Of Sequences
===========================================

Memory Mapped Files
-------------------

The fastest way to import tons of sequences in Fasta/Fastq/GSeq/... file format is to avoid slow C++ I/O streams and instead map the file directly into memory.
This can be done via the :dox:`MMapString` which uses memory mapping of the operating system or via the :dox:`ExternalString` which emulates memory mapping by doing the paging by-hand.
Most commonly used file formats concatenate sequences separated by a delimiter, e.g. ``>``, ``@``, line-break, that marks the begin of each sequence.
In SeqAn there is also a data structure that represents multiple sequences using one concatenation string and the begin positions of each sequence, the :dox:`ConcatDirectStringSet`.
We therefore defined the type :dox:`MultiSeqFile` as an alias for a :dox:`ConcatDirectStringSet` using a single :dox:`MMapString`.

In the next example we are going to open a sequence file, recognize its format, split the file into sequence fractions and import each sequence, its quality values and id.

.. includefrags:: core/demos/howto/efficiently_import_sequences.cpp
   :fragment: includes

First we associate our sequence file with the memory mapped string underlying the :dox:`ConcatDirectStringSet` using :dox:`MMapString#open`.

.. includefrags:: core/demos/howto/efficiently_import_sequences.cpp
   :fragment: open_file

Next we guess the file format of the single concatenation string and store the result in a :dox:`AutoSeqFormat` object, which is used subsequently to select the right import function.
:dox:`split` expects a :dox:`ConcatDirectStringSet` and divides the underlying string into sequence fragments separated by a file format specific delimiter.

.. includefrags:: core/demos/howto/efficiently_import_sequences.cpp
   :fragment: guess_and_split

After calling :dox:`split` the ``multiSeqFile`` StringSet represents the sequence fragments and can be used to reserve memory for the StringSets that store sequences and ids.

.. includefrags:: core/demos/howto/efficiently_import_sequences.cpp
   :fragment: reserve

The main loop iterates over each sequence fragment and uses the functions :dox:`assignSeq`, :dox:`assignQual` and :dox:`assignSeqId` to extract sequence data, qualities and id.
The quality values are encoded in ASCII and have to be converted into integer values between 0 and 62 before assigning it to a :dox:`Dna5Q` character via :dox:`AlphabetWithQualitiesConcept#assignQualityValue`.

.. includefrags:: core/demos/howto/efficiently_import_sequences.cpp
   :fragment: read_sequences

Finally we output the number of imported sequences, the overall runtime and the first 10 sequences in Fasta format.

.. includefrags:: core/demos/howto/efficiently_import_sequences.cpp
   :fragment: output

Program Output
--------------

.. code-block:: console

   $ cd build/Release
   $ make efficiently_import_sequences
   [...]
   $ ./core/demos/howto/efficiently_import_sequences reads.fq
   Loading 1000000 sequences took 4.82109 seconds

   >HWI-EAS299_3_30MAPAAXX:6:1:1561:1481/1
   GTTTATTTCACCTCCTTTACTTGTAGTCCAGGCGGTA
   >HWI-EAS299_3_30MAPAAXX:6:1:1561:1481/2
   AAAGAATTTAAATATTTCCTTAATAAGGCACGCCGTT
   >HWI-EAS299_3_30MAPAAXX:6:1:1703:1976/1
   GTTTTGATGTACAACGCCGTTACAGGTATAGTGAGAG
   >HWI-EAS299_3_30MAPAAXX:6:1:1703:1976/2
   TTCTAAATTAAAACCTCCAGAATAAGGAACATAAGAG
   >HWI-EAS299_3_30MAPAAXX:6:1:1638:1932/1
   GAAATTTTTGAGGTTATTCGCTCTTGCAACACTTTTC
   >HWI-EAS299_3_30MAPAAXX:6:1:1638:1932/2
   CACCCATACTATTAAAGCAAGCATCGGGAAAAGTAAT
   >HWI-EAS299_3_30MAPAAXX:6:1:1726:1928/1
   GCATAATGCAAAGGGTTAGTATATGATTTTTAGTATG
   >HWI-EAS299_3_30MAPAAXX:6:1:1726:1928/2
   GAGACGACAACTCCCTCCGGGAACTAAACGTGCGTAT
   >HWI-EAS299_3_30MAPAAXX:6:1:720:1208/1
   GCATATTCTATAAATGCTAAGCATAAAAATAATTTTC
   >HWI-EAS299_3_30MAPAAXX:6:1:720:1208/2
   TGCCTGTTTACCATTTAGACAGGGTTCACAAATTTCA

Remarks
-------

* We intentionally use :dox:`ContainerConcept#appendValue` to fill the StringSets as for some applications it is more memory efficient to use a :dox:`ConcatDirectStringSet` to store imported sequences and ids.
  The :dox:`ConcatDirectStringSet` consists of only one :dox:`String` concatenating all sequences and a String containing the begin positions which induce less overhead compared to storing millions of single Strings separately on heap with their own begin, end and capacity information.
* Although not visible in the example, the import functions can of course also import large sequences spanning multiple lines in various formats.

Fragment Store
--------------

The whole program above is condensed into the function :dox:`FragmentStore#loadReads` working on a :dox:`FragmentStore`.
An example for this function is given in :ref:`how-to-filter-similar-sequences`.
