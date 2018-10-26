.. sidebar:: ToC

    .. contents::

.. _how-to-use-cases-simple-genome-alignment:

Simple Genome Alignment
=======================

Learning Objective
 You will learn how to write a simple genome aligner for large-scale sequences.

Difficulty
  Medium

Duration
  4h

Prerequisites
  :ref:`tutorial-getting-started-parsing-command-line-arguments`, :ref:`tutorial-algorithms-alignment-pairwise-sequence-alignment`, :ref:`tutorial-algorithms-seed-extension`, :ref:`tutorial-datastrucures-indices-q-gram-index`, :ref:`tutorial-io-sequence-io`

Introduction
""""""""""""

Pairwise sequence alignment is the tool of choice if one wants to compare two biological sequences.
However, for genome sized sequences the standard approach is to inefficient as it requires quadratic time and space to compute an alignment.
Alternatively, one could compute only a small band but it is unclear whether evolutionary events such as insertions or deletions would shift the band away from the main diagonal.
Thus the alignment would not reflect the actual alignment of the sequences.
To cope with this issue different methods have been developed that can efficiently deal with large-scale sequences.
One of these applications was the LAGAN (Limited Area Global Alignment of Nucleotides) [1.].
It is an iterative algorithm with the following three steps.

    #. Generate local alignments with CHAOS chaining [2.].
    #. Construct a rough global map
    #. Compute final alignment along the global map

The CHAOS chaining finds local alignments depending on three parameters: the seed size (of the q-grams), the distance and the gap parameter.
The distance parameter limits the distance between two endpoints of two neighbouring seeds along the diagonal.
The band limits the shift of two seeds in vertical and horizontal direction of the matrix.
In order to find a good tradeoff between speed and sensitivity the algorithm performs step 1 and 2 recursively until two neighbouring anchors are less than a threshold apart.
In the first iteration the parameters are set more restrictive (large seed size, smaller distance and gap size.)
For every gap that is bigger than the given threshold, the parameters are set more permissive (smaller seed size and larger distance and gap sizes).

Finally, after the global map has been constructed, the algorithm performs a limited alignment around the anchors and connects the anchors with a standard alignment algorithm.

Main method and ArgumentParser
""""""""""""""""""""""""""""""

In this tutorial we want to implement a simplified version of the LAGAN approach using containing only of the first iteration.
At first we will look at the basic main method and setup the argument parser.

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: initial_main

We create a ``LaganOption`` class where we store the arguments passed to our tool:

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: lagan_option

After that, we parse given command line arguments using SeqAn's :dox:`ArgumentParser`.
Firstly, we include the source code for the argument parser by adding the following include directory:

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: include_arg_parse

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: parse_arguments

.. hint::

    If you want to learn more about parsing arguments with SeqAn read the :ref:`tutorial-getting-started-parsing-command-line-arguments` tutorial.

With this, we have set up our initial tool. Let's start implementing the algorithm.
To do so, we need to first load the sequences using the class :dox:`SeqFileIn`.
We can get access to the data structures and methods by including the following modules:

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: include_seq_io

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

   Type
     Application

   Objective
     Use the code template below (click **more...**) and implement the function ``loadSequence`` to load a single sequence file from the specified path.
     Use the file paths given in the options object and report an error if the files could not be opened.

     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
            :fragment: load_sequence_template

   Hint
     * :dox:`SeqFileIn` constructor accepts a c-style string.
     * Use `string::c_str <http://www.cplusplus.com/reference/string/string/c_str>`_ to convert the option strings into C-style strings.
     * The function :dox:`SeqFileIn#readRecord` expects the input file, a sequence, e.g. :dox:`Dna5String` and an id, e.g. :dox:`CharString`.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
           :fragment: load_sequence_solution

Finally, we can update our main method and use our ``loadSequence`` function to load sequence 1 ...

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: load_seq_1

and sequence 2 ...

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: load_seq_2

Filtering seeds
"""""""""""""""

After we read the sequences from the command line it is time to write our actual algorithm.
A naive algorithm would scan for every q-gram of sequence 2 the complete sequence 1 to find all possible positions.
But this approach of course is to slow for large-scale sequences and we need a better strategy.
SeqAn provides for this kind of task indexes which can be queried efficiently to find all occurrences of a pattern in an indexed text.

.. hint::

    There are several index implementations and we recommend to read :ref:`tutorial-datastructures-indices` to learn more about the available index data structures.

In this tutorial we are going to use a :dox:`IndexQGram` (see :ref:`tutorial-datastrucures-indices-q-gram-index` for more information), which can be included with the following module:

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: include_index

This index type will create a directory with all distinct q-grams and stores the positions of the indexed sequence, where a specific q-gram occurs.
It therefor will generate a suffix array, which is sorted by the first ```q`` symbols of every suffix.

The following line declares our q-gram index type:

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: declare_index

The ``Index`` type is a template class which requires two type template parameters.
The first type template parameter names the sequence type that this index is constructed for.
In our case this will be a `Dna5String`.
The second type template parameter is a tag or also known a policy, that defines the type of index to use.
In our case we use the :dox:`IndexQGram` policy, which itself can be further specialized through two type template parameters.
We need to select the policy used for the q-gram shape and the policy for managing the q-grams.
In our example we will need a :dox:`SimpleShape`, which is a variable length ungapped shape.
Thus, we are able to change the size of the q-gram at runtime.
Note, there also constant length shapes, whose sizes are fixed at compile time.
And as a storage policy we use :dox:`OpenAdressingTags#OpenAdressing`, which allows us to use longer values for our q-gram.
Alternatively, we could leave the parameter unspecified and would therefor enable the default behavior which is direct addressing.
Direct addressing, however, will create a lookup table for every possible q-mer (:math:`\Sigma^{q}`), which can become quite large for small q already.

In the next step we are going to initialize the index.

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: init_index

First the index is constructed with the sequence it should be created for.
Note, that this will not yet create the index.
The creation will be triggered in a lazy manner, which means it will be first created when an access to the index is requested.
Before the index can be created we need to give the index the size of the q-gram shape.
This is done in the second line of the above snippet.
The method :dox:`IndexQGram#indexShape` returns a reference to the shape stored within the index.
We resize this shape to the requested length.
The last line initializes the index shape with the first q-gram of the second sequence.

Assignment 2
^^^^^^^^^^^^

.. container:: assignment

   Type
     Application

   Objective
     Write a loop over the second sequence and write the number of occurrences per q-gram to the console.

   Hint
     * Use the function :dox:`Shape#hashNext` to update the shape for the current q-gram.
     * Use the function :dox:`IndexQGram#getOccurrences` to get a list of hits.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
           :fragment: solution_assignment2

Local chaining
""""""""""""""

Now we can stream over the second sequence and can extract all locations of a given q-gram in the indexed sequence.
To implement the second step of the LAGAN algorithm, we need to apply local chaining to the filtered q-grams and extend them to longer anchors.
SeqAn offers a data structure called :dox:`SeedSet` for this.
The following snippet shows the declaration of a ``SeedSet``:

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: include_seeds

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: declare_seed_set

The first line declares the type of the seed we want to use.
We can define the chaining policy as type template parameter.
In our case we use the :dox:`ChainedSeed` policy, which enables us to locally chain the seeds.
In addition we define a :dox:`SeedSet` with our ``ChainedSeed`` as type template parameter.

Now we create an instance of the seed set and of a scoring scheme, which we will need to score the local chain.

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: init_seed_set

Assignment 3
^^^^^^^^^^^^

.. container:: assignment

   Type
     Application

   Objective
     Update the loop from assignment 2 and fill the previously created seed set.
     Use the CHAOS chaining method to chain seeds locally, using the scoring scheme, the gap and distance criteria
     and the current position of the matching q-grams.

   Hint
     * The method :dox:`SeedSet#addSeed` has different overloads for various chaining policies.
     * If the seed could not be combined to any other in the set it must be added as a single seed to set.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
           :fragment: solution_seeding

Global chaining
"""""""""""""""

After scanning the second sequence and filling the seed set the highest scoring global chain must be computed,
which gives a map of good matching anchors.

Assignment 4
^^^^^^^^^^^^
.. container:: assignment

   Type
     Application

   Objective
     Read the documentation to :dox:`chainSeedsGlobally` and build a global chain over
     the local anchors stored in the current seed set.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
           :fragment: chain_seeds

Final alignment
"""""""""""""""

In the original algorithm, the steps from above would be repeated for the gaps between the anchors selected by the global chaining algorithm.
In this tutorial we skip the iterative step and directly compute the final alignment along the global map produced by the chaining algorithm.
SeqAn already offers an alignment function for filling the gaps and connecting them with the anchors, which is available in the ``seeds`` module.

Assignment 5
^^^^^^^^^^^^
.. container:: assignment

   Type
     Application

   Objective
     Compute an alignment around the global anchors identified by the chaining algorithm.

   Hint
     * Use the :dox:`bandedChainAlignment` function to compute the alignment.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
           :fragment: build_alignment

Finally we can output the alignment in the specified output file.

.. includefrags:: demos/tutorial/simple_genome_alignment/lagan.cpp
    :fragment: output_alignment

Congratulation, you wrote a simple genome aligner!!!
You can use the code to add iterative steps to make it more sensitive for large indels.

References
""""""""""

#. Brudno M. and Morgenstern, B., 2002. Fast and sensitive alignment of large genomic sequences. In Proceeding of the IEEE Computer Society Bioinformatics Conference (CSB).
#. Brudno, M., Do, C. B., Cooper, G. M., Kim, M. F., Davydov, E., Green, E. D., Sidow, A., Batzoglou, S., and and NISC Comparative Sequencing Program. (2003). LAGAN and Multi-LAGAN: efficient tools for large-scale multiple alignment of genomic DNA. Genome research, 13(4), 721-731.
