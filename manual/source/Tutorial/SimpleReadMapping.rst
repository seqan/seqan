.. sidebar:: ToC

   .. contents::


.. _tutorial-simple-read-mapping:

Simple Read Mapping
===================

Learning Objective
  You will be able to write read mappers using SeqAn.

Difficulty
  Hard

Duration
  2 h

Prerequisites
  :ref:`tutorial-indices`, :ref:`tutorial-fragment-store`

In this tutorial, we will walk you through the code of a simple read mapper **minimapper** that uses the SWIFT filter and uses approximate string search for verification.
There are severe limitations to its capabilities but it's a read mapper in **12 effective lines of code** (ignoring includes, comments, typedefs, I/O and lines with closing braces).

Try It
------

You can find the source code in the directory ``core/demos/tutorial/read_mapping/core/demos/tutorial/read_mapping``.
Copy over the FASTA files into your build directory and test it:

.. code-block:: console

   $ cp .../core/demos/tutorial/read_mapping/*.fasta .
   $ make demo_tutorial_minimapper
   ...
   $ ./core/demos/tutorial/read_mapping/demo_tutorial_minimapper
   Invalid number of arguments.
   USAGE: minimapper GENOME.fasta READS.fasta OUT.sam
   $ ./core/demos/tutorial/read_mapping/tutorial_minimapper nc_001454.fasta reads_hamming.fasta out.sam
   $ cat out.sam
   @HD VN:1.0
   @SQ SN:gi|9626553|ref|NC_001454.1|  LN:34214
   @PG ID:SeqAn
   nc_001454.fasta.fasta.000000005 8   gi|9626553|ref|NC_001454.1| 1396    255 36M *   0   0   TAGTGTTAGTTTATTCTGATGGAGTTGTGGAGTGAG    ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
   nc_001454.fasta.fasta.000000003 8   gi|9626553|ref|NC_001454.1| 20574   255 36M *   0   0   CCGGCGGCGTACACTGGCTGGCCCTNGCCTGGAACC    ]]]]]]]]]]]]]]]]]]]]]]]]]!]]]]]]]]]]
   nc_001454.fasta.fasta.000000007 8   gi|9626553|ref|NC_001454.1| 23191   255 36M *   0   0   GTGACAACGCGCGTTTGGCCGTACTCAAACGCACCA    ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]

Code Walkthrough
----------------

First, include the headers of the SeqAn modules we will use.

.. includefrags:: core/demos/tutorial/read_mapping/minimapper.cpp
   :fragment: includes

We will now use some typedefs for the FragmentStore and SWIFT filter to get shortcuts to types used below.

.. includefrags:: core/demos/tutorial/read_mapping/minimapper.cpp
   :fragment: typedefs

We define the global constant ``EPSILON`` (:math:`\vareps`), the allowed error rate.

.. includefrags:: core/demos/tutorial/read_mapping/minimapper.cpp
   :fragment: global-constants

Evaluate the arguments from the command line.
Use the functions :dox:`FragmentStore#loadContigs` and :dox:`FragmentStore#loadReads` to load the reference sequence (possibly more than one if the FASTA file contains more than one sequence) and reads into the FragmentStore.
Note that these functions will automatically guess the file type for you.

.. includefrags:: core/demos/tutorial/read_mapping/minimapper.cpp
   :fragment: main-input

Initialize :dox:`Finder` and :dox:`Pattern` for the q-gram index used by the swift filter.

.. includefrags:: core/demos/tutorial/read_mapping/minimapper.cpp
   :fragment: pattern-finder

Now, iterate over all input sequence contigs and enumerate all SWIFT hits.
These hits will contain all possible matches of the reads in the FragmentStore with up to :math:`\varepsilon \cdot \ell` (with :math:`\ell =` :dox:`ContainerConcept#length length(read)`) errors.
Mismatches and indels are taken into consideration.

.. includefrags:: core/demos/tutorial/read_mapping/minimapper.cpp
   :fragment: swift

Now, verify each possible match using a :dox:`HammingSimplePattern`.
The verified matches will have Hamming distance :math:`< \lfloor \varepsilon \cdot \ell \rfloor`, edit distance is not considered.

.. includefrags:: core/demos/tutorial/read_mapping/minimapper.cpp
   :fragment: verification

Finally, write out the resulting multiple read alignment to the SAM file with the file name on the command line.

.. includefrags:: core/demos/tutorial/read_mapping/minimapper.cpp
   :fragment: main-output

Hands On!
---------

Programming can only be learned by programming, so let's get started.
We create a new sandbox and a new app for the minimapper.
If you already have a sandbox, then you can skip the first step

.. code-block:: console

   $ ./util/bin/skel.py repository sandbox/my_sandbox
   $ ./util/bin/skel.py app minimapper sandbox/my_sandbox

Now, we copy over the code from the original location into our new app and build it.

.. code-block:: console

   $ cp core/demos/tutorial/read_mapping/minimapper.cpp sandbox/my_sandbox/apps/minimapper/minimapper.cpp
   $ cd build/Debug
   $ cmake .
   $ make minimapper
   $ ./sandbox/my_sandbox/apps/minimapper/minimapper
   Invalid number of arguments.
   USAGE: minimapper GENOME.fasta READS.fasta OUT.sam

Now, play around with the source code.
Here are some examples for things to try out.
There are no solutions, and they are merely thought to get you started playing...

Task 1: Use the ArgumentParser
""""""""""""""""""""""""""""""

Global constants are kind of inflexible.
Instead of the global constant *EPSILON*, create an *Options* struct with a member variable *epsilon*, initialize it to 0.8 in the constructor and use an *Option* struct in the main program.
Make the value for configurable using the class :dox:`ArgumentParser` described in the :ref:`tutorial-parsing-command-line-arguments` Tutorial.

Task 2: Allow Edit Distance for Verification
""""""""""""""""""""""""""""""""""""""""""""

Currently, the read mapper can only find reads with mismatches but not
with indels. The SWIFT filter will already create hits for positions
with indels so you only have to adjust the verification step.

Hint
  Use the :dox:`MyersPattern Myers Pattern` for the approximate search.
  Don't forget to call :dox:`Finder#findBegin` using the score (:dox:`MyersPattern#getScore`) of the last hit as the find begin score.
  You can use one Myers Pattern object per read sequence to only perform the precomputation once.
  If you reuse your finder object, don't forget to call :dox:`Finder#clear`.

Task 3: Find Matches On Reverse-Complement
""""""""""""""""""""""""""""""""""""""""""

Another limitation is that only reads from the forward strand will be found.
Either reverse-complement all reads or the contigs to find reads from the reverse strand.

Maybe add options to limit searching to the forward or reverse strand.

Hint
  Reverse-complementing the contigs will be faster in practice:
  First, an index is built over the reads which would have to be built twice if the reads were complemented.
  Second, there will usually be more reads data than genome data if the coverage is greater than 1.

Task 4: Allow Other Output Formats
""""""""""""""""""""""""""""""""""

Read the documentation on the function :dox:`FragmentStore#write` of the class :dox:`FragmentStore`.
