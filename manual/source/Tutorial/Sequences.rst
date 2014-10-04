.. sidebar:: ToC

   .. contents::


.. _tutorial-sequences:

Sequences
---------

Learning Objective
  You will learn about the SeqAn sequence concept and its main class :dox:`String` as well as the class :dox:`Segment`.
  After completing this tutorial, you will be able to use important functionalities of sequences in SeqAn and you will be ready to continue with the more specific tutorials, e.g. :ref:`tutorial-iterators`, :ref:`tutorial-alignment-representation`, or :ref:`tutorial-pairwise-sequence-alignment`.

Difficulty
  Very basic

Duration
  45 min

Prerequisites
  Basic C or C++ knowledge, the :ref:`tutorial-first-steps-in-seqan` tutorial helps.

Sequences are the core concept of SeqAn.
A sequence is a container that stores an ordered list of values.
In SeqAn, there are three kinds of sequences: Strings, Sequence Adaptions and Segments.

The :dox:`String` class is one of the most fundamental classes in SeqAn.
It is designed as a generic data structure that can be instantiated for all kinds of values, both simple (e.g. ``char``, :dox:`Dna`, :dox:`AminoAcid`) and non-simple value types (e.g. :dox:`Tuple`, :dox:`String`).
With sequence adaptions, SeqAn offers an interface for accessing data types that are not part of SeqAn, namely standard library strings and c-style char arrays.
Thus those built-in types can be handled in a similar way as SeqAn strings, for example with the :dox:`ContainerConcept#length` function.
:dox:`Segment Segments` are contiguous subsequences that represent parts of other sequences.

This tutorial will deal with the SeqAn sequence classes :dox:`String` and :dox:`Segment`.

Strings
~~~~~~~

In this section, we will have a detailed look at the SeqAn class :dox:`String`.
You will learn how to build and expand strings as well as how to compare and convert them.

Building Strings
^^^^^^^^^^^^^^^^

Let's first have a look at an example on how to define a :dox:`String`.
The type of the contained value is specified by the first template argument, e.g. ``char`` or ``int``.

.. code-block:: cpp

   String<char>  myText;     // A string of characters.
   String<int>   myNumbers;  // A string of integers.

Any type that provides a default constructor, a copy constructor and an assignment operator can be used as the alphabet / contained type of a :dox:`String`.
This includes the C++ `POD types <http://www.parashift.com/c++-faq-lite/intrinsic-types.html#faq-26.7>`_, e.g. ``char``, ``int``, ``double`` etc., but you can use more complex types, e.g. :dox:`String Strings`, too.

.. code-block:: cpp

   String<String<char> >   myStringList;   // A string of character strings.

.. hint::

   Nested Sequences (aka "Strings of Strings")

   A set of sequences can either be stored in a sequence of sequences, for example in a ``String< String<char> >``, or in :dox:`StringSet`.
   See the tutorial :ref:`tutorial-string-sets` for more information about the class :dox:`StringSet`.

SeqAn also provides the following types that are useful in bioinformatics: :dox:`AminoAcid`, :dox:`Dna`, :dox:`Dna5`, :dox:`DnaQ`, :dox:`Dna5Q`, :dox:`Finite`, :dox:`Iupac`, :dox:`Rna`, :dox:`Rna5`.
You can find detailed information in the tutorial :ref:`tutorial-alphabets`.

.. code-block:: cpp

   String<Dna>        myGenome;   // A string of nucleotides.
   String<AminAcid>   myProtein;  // A string of amino acids.

For commonly used string parameterizations, SeqAn has a range of shortcuts implemented, e.g. :dox:`DnaString`, :dox:`RnaString` and :dox:`Peptide`.

.. code-block:: cpp

   // Instead of String<Dna> dnaSeq we can also write:
   DnaString dnaSeq = "TATA";

The user can specify the kind of string that should be used in an optional second template argument of :dox:`String`.
This is also known as selecting the specialization of a class in SeqAn.
The default string implementation is :dox:`AllocString Alloc String`, which the best choice for most cases.

.. code-block:: cpp

   String<Dna>              myGenome;   // A default string of nucleotides.
   String<Dna, Alloc<> >    myGenome;   // The same as above.

For some scenarios though, there might be other types more suitable.
One such example is when processing extremely large strings that are much larger than the available main memory.
In this case, using :dox:`ExternalString External Strings` is a good choice.

.. code-block:: cpp

   // Most of the string is stored on the disk.
   String<Dna, External<> > myLargeGenome;

More details about the different specializations you can find in the tutorial :ref:`tutorial-sequences-in-depth`.

.. tip::

   String Simplify Memory Management

   One advantage of using Strings is that the user does not need to reserve memory manually with **new** and does not need **delete** to free memory.
   Instead, those operations are automatically handeld by the :dox:`String` class.

   .. code-block:: cpp

      String<Dna> myGenome = "TATACGCG";

Functionality
^^^^^^^^^^^^^

SeqAn also provides the common C++ operators for strings. You can use
them like STL strings, for example:

.. code-block:: cpp

   String<Dna> dnaSeq = "TATA";
   dnaSeq += "CGCG";
   std::cout << dnaSeq << std::endl;

.. code-block:: console

   TATACGCG

Each sequence object has a capacity, i.e. the maximum length of a sequence that can be stored in this object.
While some sequence types have a fixed capacity, the capacity of other sequence classes like :dox:`AllocString Alloc String` or ``std::basic_string`` can be changed at runtime.
The capacity can be set explicitly by functions such as :dox:`String#reserve` or :dox:`SequenceConcept#resize`.
It can also be set implicitly by functions like :dox:`SequenceConcept#append` or :dox:`SequenceConcept#replace`, if the operation's result exceeds the length of the target string.

In the following example, a :dox:`String` of :dox:`Dna5String`, we first set the new length of the container with :dox:`SequenceConcept#resize` to two elements.
After assigning two elements we append one more element with :dox:`SequenceConcept#appendValue`.
In the last step the capacity is implicitly changed.

.. code-block:: cpp

   String<Dna5String> readList;
   resize(readList, 2);
   readList[0] = "GGTTTCGACG";
   readList[1] = "AAGATGTCGC";
   appendValue(readList, "TATGCATGAT");

Using the function :dox:`ContainerConcept#length`, we can now get the length of our strings, e.g.:

.. code-block:: cpp

   std::cout << length(readList) << std::endl;
   std::cout << length(readList[0]) << std::endl;

.. code-block:: console

   3
   10

To empty a :dox:`String`, the function :dox:`SequenceConcept#clear` resets the object.

.. code-block:: cpp

   clear(readList);

SeqAn offers a range of other functions for the work with the :dox:`String` class, e.g. :dox:`AssignableConcept#assign`, :dox:`RandomAccessContainerConcept#assignValue`, :dox:`RandomAccessContainerConcept#value`, :dox:`IteratorAssociatedTypesConcept#getValue`, :dox:`ContainerConcept#empty`, etc.
The full list of functions you can find in the documentation :dox:`String`.

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     In the following assignment, you will write a small function that builds the reverse complement of a given string.
     Copy the code below and add the following functionalities:

     #. Use the ``resize`` function to ``resize`` the ``revComplGenome`` variable.
     #. Using the ``getRevCompl`` function, get the reverse complement for every nucleotide ``genome`` and store it in reverse order ``revComplGenome``.
     #. Print out the original genome and the reverse complement.

     .. code-block:: cpp

        #include <seqan/sequence.h>
        #include <seqan/basic.h>
        #include <iostream>
        #include <seqan/file.h>
        #include <seqan/modifier.h>

        using namespace seqan;

        Dna getRevCompl(Dna const & nucleotide)
        {
            if (nucleotide == (Dna)'A')
                return (Dna)'T';
            if (nucleotide == (Dna)'T')
                return (Dna)'A';
            if (nucleotide == (Dna)'C')
                return (Dna)'G';
            return (Dna)'C';
        }

        int main()
        {
            DnaString genome = "TATATACGCGCGAGTCGT";
            DnaString revComplGenome;

            // Your code snippet

            // And to check if your output is correct,
            // use the given SeqAn function reverseComplement(),
            // which modifies the sequence in-place
            reverseComplement(genome);
            std::cout << genome << std::endl;
            return 0;
        }

   Hints
     Remember that the last element in ``genome`` is stored at position ``length(genome) - 1``.

   Solution
     Click *more...* to see the solution.

     .. container:: foldable

        .. code-block:: cpp

           #include <seqan/sequence.h>
           #include <seqan/basic.h>
           #include <iostream>
           #include <seqan/file.h>
           #include <seqan/modifier.h>

           using namespace seqan;

           Dna getRevCompl(Dna const & nucleotide)
           {
               if (nucleotide == (Dna)'A')
                   return (Dna)'T';
               if (nucleotide == (Dna)'T')
                   return (Dna)'A';
               if (nucleotide == (Dna)'C')
                   return (Dna)'G';
               return (Dna)'C';
           }

           int main()
           {
               DnaString genome = "TATATACGCGCGAGTCGT";
               DnaString revComplGenome;
               resize(revComplGenome, length(genome));

               for (unsigned i = 0; i < length(genome); ++i)
               {
                   revComplGenome[length(genome) - 1 - i] = getRevCompl(genome[i]);
               }

               std::cout << genome << std::endl;
               std::cout << revComplGenome << std::endl;
               reverseComplement(genome);
               std::cout << genome << std::endl;
               return 0;
           }

Workshop Assignment 1
^^^^^^^^^^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     In this assignment, you will do some simple string building tasks, and write a simple alignment of the given reads and chromosomes.
     Use the given code template to solve these subtasks:

   #. Assume we have mapped the reads to the positions 7, 100, 172, and 272 in 'chr1'.
      Store these positions in another string 'alignPosList'.
   #. Build another String bsChr1 as a copy of chr1, and exchange every 'C' with a 'T', as in a bisulfite treated genome.
   #. Print alignments of the reads and chr1 (or bschr1) using the function ``printAlign`` and the string ``alignPosList``.

   .. code-block:: cpp

      #include <iostream>
      #include <seqan/sequence.h>
      #include <seqan/file.h>

      using namespace seqan;
      // Function to print simple alignment between two sequences with the same length
      template <typename TText1, typename TText2>
      void printAlign(TText1 const & genomeFragment, TText2 const & read)
      {
              std::cout <<  "Alignment " << std::endl;
              std::cout << "  genome : " << genomeFragment << std::endl;
              std::cout << "  read   : " << read << std::endl;
      }

      int main(int, char const **)
      {
          // Build reads and genomes
          DnaString chr1 = "TATAATATTGCTATCGCGATATCGCTAGCTAGCTACGGATTATGCGCTCTG"
                           "CGATATATCGCGCTAGATGTGCAGCTCGATCGAATGCACGTGTGTGCGATC"
                           "GATTAGCGTCGATCATCGATCTATATTAGCGCGCGGTATCGGACGATCATA"
                           "TTAGCGGTCTAGCATTTAG";

          // Build List containing all reads
          typedef String<DnaString> TDnaList;
          TDnaList readList;
          resize(readList, 4);
          readList[0] = "TTGCTATCGCGATATCGCTAGCTAGCTACGGATTATGCGCTCTGCGATATATCGCGCT";
          readList[1] = "TCGATTAGCGTCGATCATCGATCTATATTAGCGCGCGGTATCGGACGATCATATTAGCGGTCTAGCATT";
          readList[2] = "AGCCTGCGTACGTTGCAGTGCGTGCGTAGACTGTTGCAAGCCGGGGGTTCATGTGCGCTGAAGCACACATGCACA";
          readList[3] = "CGTGCACTGCTGACGTCGTGGTTGTCACATCGTCGTGCGTGCGTACTGCTGCTGACA";

          // Append a second chromosome sequence fragment to chr1
          DnaString chr2 = "AGCCTGCGTACGTTGCAGTGCGTGCGTAGACTGTTGCAAGCCGGGGGTTCAT"
                           "GTGCGCTGAAGCACACATGCACACGTCTCTGTGTTCCGACGTGTGTCACGTG"
                           "CACTGCTGACGTCGTGGTTGTCACATCGTCGTGCGTGCGTACTGCTGCTGAC"
                           "ACATGCTGCTG";
          append(chr1, chr2);

          // Print readlist
          std::cout << " \n Read list: " << std::endl;
          for(unsigned i = 0; i < length(readList); ++i)
              std::cout << readList[i] << std::endl;

          // Assume we have mapped the 4 reads to chr1 (and chr2) and now have the mapping start positions (no gaps).
          // Store the start position in a String alignPosList: 7, 100, 172, 272

          // Optional
          // Bisulfite conversion
          // Assume chr1 is beeing bisulfate treated: Copy chr1 to a new genome bsChr1 and exchange every 'C' with a 'T'
          DnaString bsChr1;

          // Print alignments of the reads with chr1 (or bsChr1) sequence using the function printAlign
          // and the positions in alignPosList.
          // To do that, you have to create a copy of the fragment in chr1 (bsChr1) that is aligned to the read.
          std::cout << " \n Print alignment: " << std::endl;
          for(unsigned i = 0; i < length(readList); ++i)
          {
              // Temporary copy of begin position (beginPosition) from alignPosList
              // of a given alignment between the read and the genome

              // Genome fragment
              DnaString genomeFragment;
              // We have to create a copy of the corresponding fragment of the genome, where the read aligns to


              // Call of our function to print the simple alignment
              printAlign(genomeFragment, readList[i]);
          }
          return 1;
      }

   Hints
     You have to create a copy of the fragment in chr1 (bsChr1) that is aligned to the read.

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: core/demos/tutorial/sequences/string_assignment_1_solution.cpp

Comparisons
^^^^^^^^^^^

Two sequences can be lexicographically **compared** using standard operators such as ``<`` or ``>=``.

.. code-block:: cpp

   String<char> a = "beta";
   String<char> b = "alpha";

   std::cout << (a != b) << std::endl;
   std::cout << (a < b) << std::endl;
   std::cout << (a > b) << std::endl;

.. code-block:: console

   1
   0
   1

Each comparison involves a scan of the two sequences for searching the first mismatch between the strings.
This could be costly if the two sequences share a long common prefix.
Suppose we want to branch in a program depending on whether ``a < b``, ``a == b``, or ``a > b``.

.. code-block:: cpp

   if (a < b)      { /* code for case "a < b"  */ }
   else if (a > b) { /* code for case "a > b"  */ }
   else            { /* code for case "a == b" */ }

In this case, although only one scan would be enough to decide what case is to be applied, each operator ``>`` and ``<`` performs a new comparison.
SeqAn offers the class :dox:`Lexical` to avoid unnecessary sequence scans.
Lexicals can store the result of a comparison, for example:

.. code-block:: cpp

   // Compare a and b and store the result in comp
   Lexical<> comp(a, b);

   if (isLess(comp))         { /* code for case "a < b"  */ }
   else if (isGreater(comp)) { /* code for case "a > b"  */ }
   else                      { /* code for case "a == b" */ }

Conversions
^^^^^^^^^^^

A sequence of type A values can be converted into a sequence of type B values, if A can be converted into B.
SeqAn offers different conversion alternatives.

**Copy conversion.**
The source sequence is copied into the target sequence.
This can be done by assignment (``operator=``) or using the function :dox:`AssignableConcept#assign`.

.. code-block:: cpp

   String<Dna> source = "acgtgcat";
   String<char> target;
   assign(target, source);
   std::cout << target;

.. code-block:: console

   acgtgcat

**Move conversion.**
If the source sequence is not needed any more after the conversion, it is always advisable to use :dox:`AssignableConcept#move` instead of :dox:`AssignableConcept#assign`.
The function :dox:`AssignableConcept#move` does not make a copy but can reuse the source sequence storage.
In some cases, :dox:`AssignableConcept#move` can also perform an in-place conversion.

.. code-block:: cpp

   String<char> source = "acgtgcat";
   String<Dna> target;

   // The in-place move conversion.
   move(target, source);
   std::cout << target;

.. code-block:: console

   acgtgcat

Assignment 2
^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     In this assignment you will sort nucleotides.
     Copy the code below. Adjust the code such that all nucleotides, which are lexicographically smaller than a Dna5 ``'G'`` are stored in a list ``lesser``, while all nucleotides which are greater, should be stored in a list ``greater``.
     Print out the final lists.

     .. code-block:: cpp

        #include <iostream>
        #include <seqan/sequence.h>
        #include <seqan/file.h>

        using namespace seqan;

        int main()
        {
            String<Dna5> nucleotides = "AGTCGTGNNANCT";
            String<Dna5> selected;
            // Append all elements of nucleotides, apart of Gs,
            // to the list selected.
            for (unsigned i = 0; i < length(nucleotides); ++i){
                appendValue(selected, nucleotides[i]);
            }
            std::cout << "Selected nucleotides: " << selected << std::endl;
            return 0;
        }

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. code-block:: cpp

           #include <iostream>
           #include <seqan/sequence.h>
           #include <seqan/file.h>

           using namespace seqan;

           int main()
           {
               String<Dna5> nucleotides = "AGTCGTGNNANCT";
               String<Dna5> lesser;
               String<Dna5> greater;

               for (unsigned i = 0; i < length(nucleotides); ++i){
                   if (nucleotides[i] < 'G')
                       appendValue(lesser, nucleotides[i]);
                   else if (nucleotides[i] > 'G')
                       appendValue(greater, nucleotides[i]);
               }
               std::cout << "Lesser nucleotides: " << lesser << std::endl;
               std::cout << "Greater nucleotides: " << greater << std::endl;
           }

Assignment 3
^^^^^^^^^^^^

.. container:: assignment

   Type
     Transfer

   Objective
     In this task you will compare whole sequences.
     Reuse the code from above. Instead of a ``String<Dna5>`` we will now deal with a ``String<Dna5String>``.
     Build a string which contains the Dna5Strings "ATATANGCGT", "AAGCATGANT" and "TGAAANTGAC".
     Now check for all elements of the container, if they are lexicographically smaller or bigger than the  given reference sequence "GATGCATGAT" and append them to a appropriate list.
     Print out the final lists.

   Hints
     Try to avoid unnecessary sequence scans.

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: core/demos/tutorial/sequences/lexical_assignment_solution.cpp

Segments
~~~~~~~~

The following section will introduce you into the :dox:`Segment` class of SeqAn.

:dox:`Segment Segments` are contiguous subsequences that represent parts of other sequences.
Therefore, their functionality is similar to the :dox:`String` functionality.
In SeqAn, there are three kinds of segments: :dox:`InfixSegment`, :dox:`PrefixSegment`, and :dox:`SuffixSegment`.
The metafunctions :dox:`SegmentableConcept#Infix`, :dox:`SegmentableConcept#Prefix`, and :dox:`SegmentableConcept#Suffix`, respectively, return the appropriate segment data type for a given sequence type.

For prefixes, we use the function :dox:`SegmentableConcept#prefix` to build the prefix.
The first parameter is the sequence we build the prefix from, the second the **excluding** end position.
For :dox:`SegmentableConcept#infix`\ es, we have to provide both the including start and the excluding end position.
For :dox:`SegmentableConcept#suffix`\ es, the second parameter of the function denotes the including starting position of the suffix:

.. code-block:: cpp

   String<Dna> dnaSeq = "AGTTGGCATG";
   Prefix<String<Dna> >::Type pre = prefix(dnaSeq, 4);
   std::cout << "Prefix: " << pre << std::endl;

   Infix<String<Dna> >::Type inf = infix(dnaSeq, 4, 7);
   std::cout << "Infix: " << inf << std::endl;

   Suffix<String<Dna> >::Type suf = suffix(dnaSeq, 4);
   std::cout << "Suffix: " << suf << std::endl;

.. code-block:: console

   Prefix: AGTT
   Infix: GGC
   Suffix: GGCATG

Segments store a pointer on the underlying sequence object, the *host*, and an start and/or end position, depending on the type of segment.
The segment is *not* a copy of the sequence segment.

.. warning::

   Please note that it is not possible anymore to change the underlying sequence by changing the segment.
   If you want to change the host sequence, you have to explicilty modify this.
   If you want to modify only the segment, you have to explicitly make a copy of the string.

Assignment 4
^^^^^^^^^^^^

.. container:: assignment

   Type
     Application

   Objective

     In this task you will use a segment to pass over an infix of a given sequence to a function without copying the corresponding fragment.
     Use the code given below.
     Lets assume that we have given a ``genome`` and a ``read`` sequence as well as the begin position of a given alignment.
     In the main function a fragment of the Dna5String ``genome`` is copied and passed together with the Dna5String ``read`` to a ``print`` function.
     Adjust the code to use an infix of the genome, instead of copying the corresponding fragment.

     .. includefrags:: core/demos/tutorial/sequences/segments_assignment_hint.cpp

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: core/demos/tutorial/sequences/segments_assignment_solution.cpp

Workshop Assignment 2
^^^^^^^^^^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     Take the solution from the workshop assignment above and change it to use Segments for building the genome fragment.

   Hints
     Note that because ``printAlign`` uses templates, you don't have to change the function even though the type of ``genomeFragment`` is different.

   Solution
    Click **more...** to see the solution.

    .. container:: foldable

       .. includefrags:: core/demos/tutorial/sequences/segments_assignment_2_solution.cpp
