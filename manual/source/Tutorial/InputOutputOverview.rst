.. sidebar:: ToC

   .. contents::


.. _tutorial-input-output-overview:

I/O Overview
============

Learning Objective
  This article will give you an overview of the I/O in SeqAn.

Difficulty
  Basic

Duration
  20 min

Prerequisites
  :ref:`tutorial-sequences`

Overview
--------

#. Adaptions of system library provided file and stream routines to the SeqAn :dox:`StreamConcept` concept.
#. Code for tokenization and parsing.
#. Conversion from textual number representations to numeric values (aka "lexical casting").


Files
-----

Most file formats in bioinformatics are structured as lists of records.
Often, they start out with a header that itself contains different header records.
For example, the SAM format starts with an optional header where users can specify the contigs of their reference sequence.

Streams
-------

In computer science, it is common to call the abstraction to such data sources **streams**.
In SeqAn, the concept :dox:`StreamConcept` provides an interface for such stream data types.

SeqAn provides adaptions from the standard C and C++ file interfaces to the :dox:`StreamConcept` concept.
Furthermore, SeqAn provides the :dox:`Stream` class and specializations for accessing ``char`` arrays and zlib and bzip compressed files as streams.

Error Handling
--------------

Exceptions.

Next Steps
----------

If you want, you can now have a look at the API documentation of the :dox:`StreamConcept` concept as well as the documentation of the :dox:`SmartFile` class.

There are two "tracks" in this section of the tutorials which you can follow.
First, you can now read the tutorials for **already supported file formats**.

* :ref:`tutorial-sequence-io`
* :ref:`tutorial-sam-bam-io`

Second, if you want to learn how to develop **support for new file formats** then read the following article.

* :ref:`tutorial-custom-io`
