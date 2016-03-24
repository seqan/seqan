.. sidebar:: ToC

    .. contents::

.. _tutorial-algorithms-realignment:

Realignment
===========

.. TODO should be greatly fleshed out!

Learning Objective
  In this tutorial, you will learn how to refine multi-sequence alignments in a fragment store.
  This can be useful for refining multi-read alignments around indels prior to small variant calling.
  After completing the tutorial, you will be able to load reads into a fragment store and compute a realignment.

Difficulty
  Advanced

Duration
  30 min

Prerequisites
  :ref:`tutorial-datastructures-store-fragment-store`

A common task in NGS data analysis is small variant calling (SNVs or indels with a length of up to 10 bp) after the read mapping step.
Usually, one considers the "pileup" of the reads and looks for variant signatures (e.g. a certain number of non-reference characters in the aligned reads).
Usually, read mappers compute pairwise alignments of each read and the reference and store them in a SAM or BAM file.
In the absence of indels, such pairwise alignments can be converted to a multi-read alignment without problems.
However, there can be an undesired multi-read alignment around indels (Figure 1).

.. 
   Commented out because of missing picture in source dir.
   .. figure:: raw_alignment.png
       :alt: MSA as interpolated from pairwise alignments.

The task of improving such an alignment is called **realignment** and there is a small number of algorithms and tools available for realignment.
This tutorial describes the ``<seqan/realign.h>`` module which implements a variant of the ReAligner algorithm by Anson and Myers :cite:`Anson1997`.

Getting Started
---------------

Consider the following program.
It creates a fragment store and then reads a small reference (with a length of 2kb) from a FASTA file and also a SAM file with reads spanning a complex indel region at 1060 ~ 1140. 
Finally, it prints the multi-read alignment around this position using :dox:`AlignedReadLayout`.

.. includefrags:: demos/tutorial/realignment/step1.cpp

The output of the program is as follows:

.. includefrags:: demos/tutorial/realignment/step1.cpp.stdout

**Figure 1:** An example of a multi-read alignment from pairwise alignments

Performing the Realignment
--------------------------

We can now use the function :dox:`reAlignment` for performing a realignment of the reads in the fragment store.

contigID
  The numeric ID of the contig to realign.
realignmentMethod
  Whether to use linear (0) or affine gap costs (1).
  It is recommended to use affine gap costs.
bandwidth
  The bandwidth to use in the realignment step.
includeReference
  Whether or not to include the reference as a pseudo read.

The algorithm works as follows:
A profile is computed, with a count of each base and the gap character at each position in the multi-read alignment.
Each read is taken and aligned against this profile.
This is repeated until convergence.
Finally, the consensus of the multi-read alignment is written into ``store.contigStore[contigID].seq``.

The parameter ``bandwidth`` controls the bandwidth of the banded alignment used in the alignment of reads against the profile.
If ``includeReference`` is ``true`` then the reference is added as a pseudo-read (a new read at the end of the read store).
This can be used for computing alignments of the reads agains the original reference.

.. includefrags:: demos/tutorial/realignment/step2.cpp

Here is the program's output.
The reference pseudo-read is here shown as the first read (second row) below the reference (first row).

.. includefrags:: demos/tutorial/realignment/step2.cpp.stdout
