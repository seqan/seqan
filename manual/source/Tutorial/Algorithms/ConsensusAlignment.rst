.. sidebar:: ToC

    .. contents::

.. _tutorial-algorithms-consensus-alignment:

Consensus Alignment
===================

Learning Objective
  You will learn how to perform a consensus alignment of sequences (e.g. NGS reads) stored in a FragmentStore.
  After completing this tutorial, you will be able to perform a consensus alignment of reads with and without using alignment information.
  This is useful for the consensus step in sequence assembly.

Difficulty
  Advanced

Duration
  1 h

Prerequisites
  :ref:`tutorial-datastructures-store-fragment-store`

The SeqAn module ``<seqan/consensus.h>`` allows the computation of consensus alignments based on the method by Rausch et al. :cite:`Rausch2009`.
It can be used for the consensus step in Overlap-Layout-Consensus assemblers.

Consensus with Approximate Positions
------------------------------------

The consensus module has two modes.
The first one is applicable when approximate positions of the reads are known.
The following program demonstrates this functionality.

First, we include the necessary headers.

.. includefrags:: demos/tutorial/consensus_alignment/with_positions.cpp
   :fragment: includes

Next, the fragment store is filled with reads and approximate positions.
The true alignment is shown in the comments.

.. includefrags:: demos/tutorial/consensus_alignment/with_positions.cpp
   :fragment: fill_store

This is followed by computing the consensus alignment using the function :dox:`consensusAlignment`.

.. includefrags:: demos/tutorial/consensus_alignment/with_positions.cpp
   :fragment: compute_consensus

Finally, the alignment is printed using an :dox:`AlignedReadLayout` object.

.. includefrags:: demos/tutorial/consensus_alignment/with_positions.cpp
   :fragment: print_layout

Here is the program's output:

.. includefrags:: demos/tutorial/consensus_alignment/with_positions.cpp.stdout

Consensus without Approximate Positions
---------------------------------------

When setting the ``useContigID`` member of the :dox:`ConsensusAlignmentOptions` object to ``false`` then we can also omit adding approximate positions for the reads.
In this case, the consensus step performs an all-to-all alignment of all reads and then computes a consensus multi-read alignment for all of them.
This is demonstrated by the following program.


.. includefrags:: demos/tutorial/consensus_alignment/without_positions.cpp

Here is this modified programs' output:

.. includefrags:: demos/tutorial/consensus_alignment/without_positions.cpp.stdout
