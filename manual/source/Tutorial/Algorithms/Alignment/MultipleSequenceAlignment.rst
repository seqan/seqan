.. sidebar:: ToC

    .. contents::

.. _tutorial-algorithms-alignment-multiple-sequence-alignment:

Multiple Sequence Alignment
===========================

Learning Objective
  You will learn how to compute a multiple sequence alignment (MSA) using SeqAn's alignment data structures and algorithms.

Difficulty
  Basic

Duration
  30 min

Prerequisites
  :ref:`tutorial-datastructures-sequences`, :ref:`tutorial-datastructures-alignment`

Alignments are at the core of biological sequence analysis and part of the "bread and butter" tasks in this area.
As you have learned in the :ref:`pairwise alignment tutorial <tutorial-algorithms-alignment-pairwise-sequence-alignment>`, SeqAn offers powerful and flexible functionality for coputing such pairwise alignments.
This tutorial shows how to compute multiple sequence alignments (MSAs) using SeqAn.
First, some background on MSA will be given and the tutorial will then explain how to create multiple sequence alignments.

Note that this tutorial focuses on the ``<seqan/graph_msa.h>`` module whose purpose is the computation of **global** MSAs, i.e. similar to SeqAn::T-Coffe :cite:`Rausch2008` or ClustalW :cite:`Thompson1994`.
If you are interested in computing consensus sequences of multiple overlapping sequences (e.g. NGS reads), similar to assembly after the layouting step, then have a look at the :ref:`tutorial-algorithms-consensus-alignment` tutorial.

While the pairwise alignment of sequences can be computed exactly in quadratic time usind dynamic programming, the computation of exact MSAs is harder.
Given :math:`n` sequences of length :math:`\ell`, the exact computation of an MSA is only feasible in time :math:`\mathcal{O}(\ell^n)`.
Thus, global MSAs are usually computed using a heuristic called **progressive alignment**.
For an introduction to MSAs, see the `Wikipedia Article on Multiple Sequence Aligment <http://en.wikipedia.org/wiki/Multiple_sequence_alignment>`_.

Computing MSAs with SeqAn
-------------------------

The SeqAn library gives you access to the engine of SeqAn::T-Coffee :cite:`Rausch2008`, a powerful and efficient MSA algorithm based on the progressive alignment strategy.
The easiest way to compute multiple sequence alignments is using the function :dox:`globalMsaAlignment`.
The following example shows how to compute a global multiple sequence alignment of proteins using the :dox:`Blosum62` scoring matrix with gap extension penalty ``-11`` and gap open penalty ``-1``.

First, we include the necessary headers and begin the ``main`` function by declaring our strings as a char array.

.. includefrags:: demos/tutorial/multiple_sequence_alignment/msa.cpp
   :fragment: main

Next, we build a :dox:`Align` object with underling :dox:`String SeqAn Strings` over the :dox:`AminoAcid` alphabet.
We create four rows and assign the previously defined amino acid strings into the rows.

.. includefrags:: demos/tutorial/multiple_sequence_alignment/msa.cpp
   :fragment: init

Finally, we call :dox:`globalMsaAlignment` and print ``align`` to the standard output.
We use the :dox:`Blosum62` score matrix with the penalties from above.

.. includefrags:: demos/tutorial/multiple_sequence_alignment/msa.cpp
   :fragment: alignment

The output of the program look as follows.

.. includefrags:: demos/tutorial/multiple_sequence_alignment/msa.cpp.stdout

Note that we stored the MSA in an :dox:`Align` object which allows easy access to the individual rows of the MSA as :dox:`Gaps` objects.
:dox:`globalMsaAlignment` also allows storing the alignment as an :dox:`AlignmentGraph`.
While this data structure makes other operations easier, it is less intuitive than the tabular represention of the :dox:`Align` class.

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     Compute a multiple sequence alignments between the four protein sequences from above using a :dox:`Align` object and the :dox:`Blosum80` score matrix.

   Solution
     .. container:: foldable

        The solution looks as follows.

        .. includefrags:: demos/tutorial/multiple_sequence_alignment/assignment1.cpp

        And here is the program's output.

        .. includefrags:: demos/tutorial/multiple_sequence_alignment/assignment1.cpp.stdout

Computing Consensus Sequences
-----------------------------

One common task following the computation of a global MSA for DNA sequences is the computation of a consensus sequence.
The type :dox:`ProfileChar` can be used for storing counts for a profile's individual characters.
It is used by creating a :dox:`String` over :dox:`ProfileChar` as the alphabet.

The following program first computes a global MSA of four variants of exon1 of the gene SHH.
First, we compute the alignment as in the example above.

.. includefrags:: demos/tutorial/multiple_sequence_alignment/consensus.cpp
   :fragment: align

Then, we create the profile string with the length of the MSA.
We then count the number of characters (and gap pseudo-characters which have an ``ordValue`` of ``4`` for :dox:`Gaps` over :dox:`Dna`) at each position.

.. includefrags:: demos/tutorial/multiple_sequence_alignment/consensus.cpp
   :fragment: profile-computation

Finally, we compute the consensus and print it to the standard output.
At each position, the consensus is called as the character with the highest count.
Note that ``getMaxIndex`` breaks ties by the ordinal value of the caracters, i.e. ``A`` would be preferred over ``C``, ``C`` over ``G`` and so on.

.. includefrags:: demos/tutorial/multiple_sequence_alignment/consensus.cpp
   :fragment: consensus-calling

The output of the program is as follows.

.. includefrags:: demos/tutorial/multiple_sequence_alignment/consensus.cpp.stdout
