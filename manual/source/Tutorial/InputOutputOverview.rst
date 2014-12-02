.. sidebar:: ToC

   .. contents::


.. _tutorial-input-output-overview:

File I/O Overview
=================

Learning Objective
  This article will give you an overview of the file I/O in SeqAn.

Difficulty
  Basic

Duration
  30 min

Prerequisites
  :ref:`tutorial-sequences`


Overview
--------

Most file formats in bioinformatics are structured as lists of records.
Often, they start out with a header that itself contains different header records.
For example, the Binary Sequence Alignment/Map (BAM) format starts with an header that lists all contigs of the reference sequence.
The BAM header is followed by a list of BAM alignment records that contain query sequences aligned to some reference contig.

SeqAn allows to read or write record-structured files through two types of classes: :dox:`FileIn` and :dox:`FileOut`.
Classes of type :dox:`FileIn` allow to read files, whereas classes of type :dox:`FileOut` allow to write files.
For example, the class :dox:`BamFileIn` allows to read a BAM file, whereas the class :dox:`BamFileOut` allows to write a BAM file.
Note how these types of classes **do not allow to read and write the same file at the same time**.

This tutorial shows the basic functionalities provided by :dox:`FileIn` and :dox:`FileOut` class types.
In particular, this tutorial adopts the classes :dox:`BamFileIn` and :dox:`BamFileOut` as concrete types.
Nonetheless, **these functionalities are independent from the particular file format** and thus valid for all record-based file formats supported by SeqAn.

Basic I/O
---------

The demo application shown here is a simple SAM to BAM converter.

Includes
""""""""

Support for a specific format comes by including a specific header file.
In this case, we include the BAM header file:

.. includefrags:: demos/tutorial/base_io/example1.cpp
   :fragment: include


Opening and Closing Files
"""""""""""""""""""""""""

Classes :dox:`FileIn` and :dox:`FileOut` allow to :dox:`open` and :dox:`close` files.

A file can be opened by passing the filename to the constructor:

.. includefrags:: demos/tutorial/base_io/example1.cpp
   :fragment: ctor

Alternatively, a file can be opened after construction by calling :dox:`open`:

.. includefrags:: demos/tutorial/base_io/example1.cpp
   :fragment: open

Noe that any file is closed *automatically* whenever the :dox:`FileIn` or :dox:`FileOut` object goes out of scope.
Eventually, a file can be closed *manually* by calling :dox:`close`.

Accessing the Header
""""""""""""""""""""

To access the header, we need an object representing the format-specific header.
In this case, we use an object of type :dox:`BamHeader`.
The content of this object can be ignored for now, it will be covered in the :ref:`tutorial-sam-bam-io` tutorial.

.. includefrags:: demos/tutorial/base_io/example1.cpp
   :fragment: header

Function :dox:`BamFileIn#readRecord` reads the header from the input SAM file and :dox:`BamFileOut#writeRecord` writes it to the BAM output file.

Accessing the Records
"""""""""""""""""""""

There are three use cases for reading or writing record-based files:

#. read or write the file **record by record**;
#. read or write a **batch of records**, e.g. 100k records at a time;
#. read or write **all records** from or to the file.

These use cases are supported respectively by the functions :dox:`readRecord` and :dox:`readRecords`, or :dox:`writeRecord` and :dox:`writeRecords`.

In this example, we are going to read and write the files record by record.
Again, to access each record, we need an object representing the format-specific record.
In this case, we use an object of type :dox:`BamAlignmentRecord`.
Each call to :dox:`BamFileIn#readRecord` reads one record from the SAM input file and moves the :dox:`BamFileIn` forward.
Each call to :dox:`BamFileOut#writeRecord` writes the record just read to the BAM output files.
We check the end of the input file by calling :dox:`BamFileIn#atEnd`.

.. includefrags:: demos/tutorial/base_io/example1.cpp
   :fragment: records

Our small SAM to BAM conversion demo is ready.
The tool still lacks error handling, reading from standard input and writing to standard output.
You are now going to add these features.

Error Handling
--------------

We distinguish between two types of errors: *low-level* file I/O errors and *high-level* file format errors.
Possible file I/O errors can affect both input and output files.
Example of errors are: the file permissions forbid a certain operations, the file does not exist, there is a disk reading error, a file being read gets deleted while we are reading from it, or there is a physical error in the hard disk.
Conversely, file format errors can only affect input files.
Such errors arise whenever the input file content is damaged or incorrect.

Error handling in SeqAn is implemented by means of exceptions.
Classes of types :dox:`FileIn` and :dox:`FileOut` throw exceptions of type :dox:`IOError` to signal *low-level* file I/O errors and exceptions of type :dox:`ParseError` to signal *high-level* input file format errors.

I/O Errors
""""""""""

All :dox:`FileIn` and :dox:`FileOut` constructors and functions throw :dox:`IOError` exceptions on failure.
Therefore, it is sufficient to catch these exceptions to handle any error properly.

There is only one exception to this rule.
Function :dox:`FileIn#open` returns a ``bool`` to indicate whether the file was opened successfully or not.


Assignment 1
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Improve the program above to handle I/O errors.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/base_io/solution1.cpp


Parsing Errors
""""""""""""""

Functions :dox:`FileIn#readRecord` and :dox:`FileIn#readRecords` throw :dox:`ParseError` exceptions on failure.


Assignment 2
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Improve the program above to handle parsing errors.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/base_io/solution2.cpp


File Formats
------------

.. warning::
    Describe file format detection of FileIn and FileOut.

Compressed Files
""""""""""""""""

All above examples and your solutions to the assignments **already have compression support built-in**, if the compression libraries are available!
For accessing compressed files, you need to have zlib installed for reading ``.gz`` files and libbz2 for reading ``.bz2`` files.

If you are using Linux or Mac Os X and you followed the :ref:`tutorial-getting-started` tutorial closely then you should have already installed the necessary libraries.
On Windows, you will need to follow :ref:`how-to-install-contribs-on-windows` to get the necessary libraries.

You can check whether you have installed the libraries to use zlib and libbz2 by running CMake again.
Simply call ``cmake .`` in your build directory.
At the end of the output, there will be a section "SeqAn Features".
If you can read ``ZLIB - FOUND`` and ``BZIP2 - FOUND`` then you can use zlib and libbz2 in your programs.



Streams
-------

The constructors of :dox:`FileIn` and :dox:`FileOut` accept not only filenames, but also standard C++ streams or any other class fulfilling the :dox:`StreamConcept` concept.
For instance, you can pass `std::cin` to any :dox:`FileIn::FileIn FileIn constructor` and `std::cout` to any :dox:`FileIn::FileOut FileOut constructor`.

Assignment 3
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Improve the program above to read from standard input.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/base_io/solution3.cpp


Assignment 4
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Improve the program above to write to standard output.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/base_io/solution4.cpp



Next Steps
----------

If you want, you can now have a look at the API documentation of the :dox:`StreamConcept` concept as well as the documentation of the :dox:`SmartFile` class.

There are two "tracks" in this section of the tutorials which you can follow.
First, you can now read the tutorials for **already supported file formats**.

* :ref:`tutorial-sequence-io`
* :ref:`tutorial-sam-bam-io`

Second, if you want to learn how to develop **support for new file formats** then read the following article.

* :ref:`tutorial-custom-io`
