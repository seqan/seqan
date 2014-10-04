.. sidebar:: ToC

   .. contents::


.. _tutorial-mini-bowtie:

Mini Bowtie
===========

Learning Objective
  You will be able to write read mappers using SeqAn.

Difficulty
  Hard

Duration
  2h

Prerequisites
  :ref:`tutorial-indices`, :ref:`tutorial-index-iterators`, :ref:`tutorial-fragment-store`

Mini-Bowtie is a very basic read aligner that is inspired by the well known Bowtie program :cite:`Langmead2009`.
It serves as an example to show that you can write sophisticated programs with SeqAn using few lines of code.

Background
----------

Bowtie is based on the FM index :cite:`Ferragina2001` which uses the Burrows-Wheeler transform (BWT) to search for patterns in a text in :math:`\mathcal{O}(m)` time with ``m`` being the pattern length.
Furthermore the index only consumes minute space because it can be compressed using favorable features of the BWT.

Since this tutorial serves as an example of SeqAn's functionality, we will not completely re-implement Bowtie.
Instead we will concentrate on providing a simple yet functional read mapper.
We will start by reading the reference sequence and reads, then search approximately each read in the reference allowing at most one mismatch, and finally write the result into an informative SAM file.

We do make use of the FM index, but we use only one key feature of Bowtie, namely a two phase search.
Such two phase search speeds up approximate search on the FM index but requires two indices, i.e. one of the forward and one of the backward reference sequence.
Instead of searching each read in one phase, allowing one mismatch in every possible position, we search it twice, but in two phases having some restrictions.
We first search each read in the forward index by allowing one mismatch only in the second half of the read.
Then we reverse each read and search it in the backward index, again by allowing one mismatch only in the second half of the read.
The advantage of a two phase search will become clear in the following.

.. tip::

    The FM index implicitly represent a prefix trie of a text.
    Searching on the FM index is performed backwards, starting with the last character of the pattern and working its way to the first.

In the first search phase we search backwards the right half of the read without mismatches.
If we succeed, in the second phase we continue backwards allowing one mismatch in the left half.
The whole procedure is based on a depth-first traversal on the prefix trie of the reference.
Trie’s edges are labeled with single characters and each path following some edges spells a substring of the reference.
Thus, first we start in the root node and follow the path spelling the right half of the read.
If we succeed, we continue in the subtree rooted in the previously reached node by following each possible path that spells the left half of the read with one substitution.

.. figure:: mini_bowtie_search.png
   :width: 600px

The red arrows indicate the exact search of the right half of the read, the green arrows the approximate search of the left half.

From the figure above you can see the advantage of using a two phase approach.
By starting from the root node and performing exact search, we follow only one path and avoid visiting the top of the trie.
Approximate search is more involved, since we have to try every path until we find a second mismatch, but it is limited to a small subtree.

So far we have found all locations of a read having at most one mismatch in its left half.
In order to find also all locations having at most one mismatch in the right half, we simply reverse the reference and the read and repeat the procedure from above.
Note that we only reverse but do not complement the reference or the reads.

Assignment 1
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Copy the code below and rearrange the rows such that they give a functional order.

     .. code-block:: cpp

        // ==========================================================================
        //                                mini_bowtie
        // ==========================================================================
        // Copyright (c) 2006-2012, Knut Reinert, FU Berlin
        // All rights reserved.


        #include <iostream>
        #include <seqan/basic.h>
        #include <seqan/sequence.h>
        #include <seqan/file.h>
        #include <seqan/index.h>
        #include <seqan/store.h>

        using namespace seqan;

        void search() {};

        int main(int argc, char \*argv[])
        {
            // type definitions
            typedef String<Dna5> TString;
            typedef StringSet<TString> TStringSet;
            typedef Index<StringSet<TString>, FMIndex<> > TIndex;
            typedef Iterator<TIndex, TopDown<ParentLinks<> > >::Type TIter;

            // reversing the sequences for backward search
            // backward search
            // reading the command line arguments
            // declaration and initialization of the fragment store
            // forward search
            // combining the contigs of the reference into one string set
            appendValue(text, fragStore.contigStore[i].seq);
            std::cerr << "Invalid number of arguments." << std::endl
                          << "USAGE: mini_bowtie GENOME.fasta READS.fasta OUT.sam" << std::endl;
            }
            if (argc < 3) {
            if (loadContigs(fragStore, argv[1])) return 1;
            if (loadReads(fragStore, argv[2])) return 1;
            clear(fmIndex);
            clear(fmIndex);
            StringSet<TString> text;
            for (unsigned i = 0; i < length(fragStore.contigStore); ++i)
            fmIndex = TIndex(text);
            TIndex fmIndex(text);
            TIter it(fmIndex);
            search();
            search();
            clear(it);
            clear(it);
            reverse(text);
            reverse(fragStore.readSeqStore);
            it = TIter(fmIndex);
            FragmentStore<> fragStore;
            return 0;
            return 1;
        }

   Hint
     .. container:: foldable

       We make use of the :dox:`FragmentStore`.
       While we can access the pattern/reads as if using a :dox:`StringSet`, we need to create a :dox:`StringSet` of the contigs, because in the :dox:`FragmentStore` the contigs are not stored in a :dox:`StringSet`.

   Hint
     .. container:: foldable

       The correct order of the comments is:

       .. code-block:: cpp

           // reading the command line arguments
           // declaration and initialization of the fragment store
           // combining the contigs of the reference into one string set
           // forward search
           // reversing the sequences for backward search
           // backward search

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/mini_bowtie/solution1.cpp

Now that we have the backbone of our program we can start to implement the fundamental part, the search routine.
The ``search`` function requires two input arguments, namely the iterator used to traverse the FM index of the reference sequence and the string set containing the reads.

The ``search`` function iterates over the reads and searches them in the already mentioned two phase fashion.
In the first phase the right half of the pattern is searched exactly.
The second phase is more involved and will be discussed after the second assignment.

Assignment 2
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Expand the solution to the last assignment by implementing the backbone of the search routine.
     The backbone should consist of function definition, an outer loop traversing the pattern (using a standard iterator) and the first step of the search, namely the exact search of the right pattern half.

   Hint
     .. container:: foldable

        The template function header could look like this:

        .. code-block:: cpp

           template <typename TIter, typename TStringSet>
           void search(TIter & it, TStringSet const & pattern)


   Hint
     .. container:: foldable

        You can obtain the correct iterator type using the metafunction :dox:`ContainerConcept#Iterator`.

        .. code-block:: cpp

           typedef typename Iterator<TStringSet const, Standard>::Type TPatternIter;

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/mini_bowtie/solution2.cpp

At this point we need to implement the most critical part of the our program, which is the second step of the search phase.
In our case this step works as follows:

Assume that we have already found a path in the trie representing the pattern from position :math:`i` to :math:`(m - 1)` with m being the pattern length and :math:`i < m / 2`.
Then we substitute the character of the pattern at position :math:`i - 1` with every character of the alphabet and search for the remaining characters exact.
We repeat those two steps until we processed every character of the pattern.

The corresponding pseudo code could look like this:

.. code-block:: cpp

   unsigned startApproxSearch = length(pattern) / 2;
   for i = startApproxSearch to 0
   {
       for all c in the aphabet
       {
           if (goDown(it, c))
           {
               if (goDown(it, pattern[0..i - 1]))
               {
                    HIT
               }
               goBack to last start
           }
       }
       goDown correct character
   }

Assignment 3
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Include the pseudo code from above into the search function.

   Hint
     Make a copy of the iterator before following the path of the substituted character.
     Doing so saves time and keeps the code simple because you do not need to use :dox:`TopDownHistoryIterator#goUp`.

   Hint
     .. container:: foldable

       :dox:`TopDownIterator#goDown` returns a boolean indicating if a path exists or not.
       In addition, you do not need to go through the steps of the pseudo code if the second pattern half was not found!

   Hint
     .. container:: foldable

       :dox:`OrderedAlphabetConcept#MinValue` returns the lowest value of an alphabet, while :dox:`FiniteOrderedAlphabetConcept#ValueSize` returns the number of different values of a data type.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/mini_bowtie/solution3.cpp

So far so good.
But there is a slight mistake.
While substituting we also substitute the character in the pattern with itself.
Therefore we find locations of exact matches multiple times.

Assignment 4
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Adjust the code to cope with the problem mentioned above.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/mini_bowtie/solution4.cpp

So this is already the fundamental part of our program.
What's left to do is to write the result into a SAM file.
In order to do so, we make use of the :dox:`FragmentStore`.
Everything we need to do is to fill the :dox:`FragmentStore::alignedReadStore` which is a member of the :dox:`FragmentStore`.
This is very easy, because we only need to append a new value of type :dox:`AlignedReadStoreElement::AlignedReadStoreElement` specifying the match id, the pattern id, the id of the contig, as well as the begin and end position of the match in the reference.

An ``addMatchToStore`` function could look like this:

.. code-block:: cpp

   template <typename TStore, typename TIter, typename TPatternIt>
   void addMatchToStore(TStore & fragStore, TPatternIt const & patternIt, TIter const & localIt)
   {
       typedef FragmentStore<>::TAlignedReadStore TAlignedReadStore;
       typedef Value<TAlignedReadStore>::Type TAlignedRead;

       for (unsigned num = 0; num < countOccurrences(localIt); ++num)
       {
           unsigned pos = getOccurrences(localIt)[num].i2;
           TAlignedRead match(length(fragStore.alignedReadStore), position(patternIt), getOccurrences(localIt)[num].i1 ,
               pos,  pos + length(value(patternIt)));
           appendValue(fragStore.alignedReadStore, match);
       }
   }

One could think we are done now.
Unfortunately we are not.
There are two problems.
Recall, that in the second search phase we reverse the text and pattern, therefore we messed up start and end positions in the original reference.
Furthermore we found exact matches twice, once in the forward index and once in the reversed index.

Assignment 5
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Include the following two lines into your code:

     .. code-block:: cpp

        struct ForwardTag {};
        struct ReverseTag {};

        and write a second ``addMatchToStore`` function that is called when we search in the reversed reference.
        In addition, make all necessary changes to the code such that exact matches are only added once.``

   Hint
     .. container:: foldable

        The metafunction :dox:`IsSameType` can be used to determine whether two types are equal or not.

   Solution
     .. container:: foldable

       .. includefrags:: extras/demos/tutorial/mini_bowtie/mini_bowtie.cpp

Done? Not quite.

We need to copy the following four lines of code into our code in order to write the correct result in the SAM file.
Calling the reverse function is necessary because an alignment must be computed for every match to be written into the SAM file.

.. code-block:: cpp

   reverse(text);
   reverse(fragStore.readSeqStore);
   std::ofstream samFile(argv[3], std::ios_base::out);
   write(samFile, fragStore, Sam());
