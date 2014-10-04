.. sidebar:: ToC

   .. contents::


.. _tutorial-multiple-sequence-alignment:

Multiple Sequence Alignment
---------------------------

Learning Objective
 You will learn how to compute a multiple sequence alignment using SeqAn's alignment data structures and algorithms.

Difficulty
  Basic

Duration
  15 min

Prerequisites
  :ref:`tutorial-first-steps-in-seqan`, :ref:`tutorial-sequences`, :ref:`tutorial-alphabets`, :ref:`tutorial-alignment-representation`

Apart from pairwise alignments, also multiple sequence alignments can be computed in SeqAn.
The easiest way to do this is by using the function :dox:`globalMsaAlignment`.
This function computes a heuristic alignment based on a consistency-based progressive alignment strategy as described in `SeqAn::TCoffee <http://bioinformatics.oxfordjournals.org/cgi/content/abstract/24/16/i187>`_ paper.

In the following example, we align four amino acid sequences using the :dox:`AlignmentGraph` data structure and the :dox:`Blosum62` scoring matrix with gap extension penalty -11 and gap open penalty -1.
The required header for multiple sequence alignments is ``<seqan/graph_msa.h>``.

.. includefrags:: core/demos/tutorial/alignments/alignment_msa.cpp
   :fragment: main

First, the sequence type ``TSequence`` is defined and a :dox:`StringSet` is declared.
The four sequences to be aligned are appended to the StringSet ``seq``.

.. includefrags:: core/demos/tutorial/alignments/alignment_msa.cpp
   :fragment: init

Now we initialize our :dox:`AlignmentGraph` with the sequences.
The graph and the :dox:`Blosum62` scoring matrix are handed to the function :dox:`globalMsaAlignment` which computes the desired alignment.

.. includefrags:: core/demos/tutorial/alignments/alignment_msa.cpp
   :fragment: alignment

And here is the output of this example program:

.. code-block:: console

   Alignment matrix:
         0     .    :    .    :    .    :    .    :    .    :
           DPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSA
             | |   |          |       |      || ||     ||
           RVKRP---MNAFIVWSRDQRRKMALENP--RMRNSEISKQLGYQWKMLTE
             | |              | | |   |   | |    | |    | | |
           FPKKP---LTPYFRFFMEKRAKYAKLHP--EMSNLDLTKILSKKYKELPE
             |||   |        | ||                  ||      |
           HIKKP---LNAFMLYMKEMRANVVAEST--LKESAAINQILGRRWHALSR

        50     .    :    .    :    .    :    .    :
           KEKGKFEDMAKADKARYEREMKTY----------IPPKGE
            ||  |   |    |        |
           AEKWPFFQEAQKLQAMHREKYPNYKYRP---RRKAKMLPK
             |    |  |                            |
           KKKMKYIQDFQREKQEFERNLARFREDH---PDLIQNAKK
               ||      | |                        |
           EEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     Compute a multiple sequence alignments between the four protein sequences

     * ``DPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE``
     * ``RVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK``
     * ``FPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK``
     * ``HIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK``

     using a :dox:`Align` object and the :dox:`Blosum80` score matrix.

     Repeat the above example using the Align data structure and the Blosum80 scoring matrix.

   Solution
     .. container:: foldable

        After the usual includes, the :dox:`Align` object `align` is initialized and the four sequences are appended as rows.

        .. includefrags:: core/demos/tutorial/alignments/alignment_msa_assignment1.cpp
           :fragment: main

        Now the MSA is computed, using the :dox:`Blosum80` matrix for scoring.

        .. includefrags:: core/demos/tutorial/alignments/alignment_msa_assignment1.cpp
           :fragment: alignment

        And here is the output:

        .. code-block:: console

              0     .    :    .    :    .    :    .    :    .    :
                DPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSA
                  | |   |          |       |    |  | ||     ||
                RVKRP---MNAFIVWSRDQRRKMALENPRMR-NS-EISKQLGYQWKMLTE
                  | |              | | |   | |  |     | |    | | |
                FPKKP---LTPYFRFFMEKRAKYAKLHPEMS-NL-DLTKILSKKYKELPE
                  |||   |        | ||                  ||      |
                HIKKP---LNAFMLYMKEMRANVVAESTLKE-SA-AINQILGRRWHALSR

             50     .    :    .    :    .    :    .    :    .    :
                KEKGKFEDMAKADKARYEREMKTY---------------IP--PKG---E
                 ||  |   |    |   || |                  |
                AEKWPFFQEAQKLQAMH-RE-K-----YP------NYKYRPRRKAKMLPK
                  |    |  |       |         |               ||   |
                KKKMKYIQDFQREKQEFERNLARFREDHP------DL--IQ--NAK---K
                    ||      | |             |                    |
                EEQAKYYELARKERQLH-MQ-L-----YPGWSARDNYGKKKKRKRE---K
