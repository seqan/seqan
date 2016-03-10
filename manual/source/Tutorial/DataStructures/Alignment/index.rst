.. _tutorial-datastructures-alignment:

Alignment
========================

..  toctree::
    :hidden:
    :titlesonly:

    ScoringSchemes
    AlignmentGaps
    AlignmentGraph

Alignment Algorithms ( e.g.
:ref:`Pairwise <tutorial-algorithms-alignment-pairwise-sequence-alignment>` and
:ref:`Multiple <tutorial-algorithms-alignment-multiple-sequence-alignment>` )
are one of the core algorithms in SeqAn. In this section you can learn how SeqAn
represents alignments as C++ Objects and how you could use those data structures
for your own alignment algorithm. Furthermore, you can learn which different
kinds of :ref:`tutorial-datastructures-alignment-scoringschemes` exist,
i.e. which combinations of Match/Mismatch Evaluation (e.g. Simple Score,
Substitutional Matrices Score) and Insertion/Deletion Evaluation (e.g. Linear
Gap Model, Affine Gap Model and Dynamic Gap Model) are possible and how you can
define your own Scoring Matrices.
