.. sidebar:: ToC

    .. contents::

.. _tutorial-io-input-output-overview:

File I/O Overview
=================

Learning Objective
  This article will give you an overview of the formatted file I/O in SeqAn.

Difficulty
  Basic

Duration
  30 min

Prerequisites
  :ref:`tutorial-datastructures-sequences`


Overview
--------

Most file formats in bioinformatics are structured as lists of records.
Often, they start out with a header that itself contains different header records.
For example, the Binary Sequence Alignment/Map (SAM/BAM) format starts with an header that lists all contigs of the reference sequence.
The BAM header is followed by a list of BAM alignment records that contain query sequences aligned to some reference contig.

.. _tutorial-io-input-output-overview-formatted-files:

Formatted Files
^^^^^^^^^^^^^^^

SeqAn allows to read or write record-structured files through two types of classes: :dox:`FormattedFileIn` and :dox:`FormattedFileOut`.
Classes of type :dox:`FormattedFileIn` allow to read files, whereas classes of type :dox:`FormattedFileOut` allow to write files.
Note how these types of classes **do not allow to read and write the same file at the same time**.

These types of classes provide the following I/O operations on formatted files:

#. Open a file given its filename or attach to an existing stream like `std::cin` or `std::cout`.
#. Guess the file format from the file content or filename extension.
#. Access compressed or uncompressed files transparently.

SeqAn provides the following file formats:

* :dox:`SeqFileIn`, :dox:`SeqFileOut` (see Tutorial :ref:`tutorial-io-sequence-io`)
* :dox:`BamFileIn`, :dox:`BamFileOut` (see Tutorial :ref:`tutorial-io-sam-bam-io`)
* :dox:`BedFileIn`, :dox:`BedFileOut` (see Tutorial :ref:`tutorial-io-bed-io`)
* :dox:`VcfFileIn`, :dox:`VcfFileOut` (see Tutorial :ref:`tutorial-io-vcf-io`)
* :dox:`GffFileIn`, :dox:`GffFileOut` (see Tutorial :ref:`tutorial-io-gff-and-gtf-io`)
* :dox:`RoiFileIn`, :dox:`RoiFileOut`
* :dox:`SimpleIntervalsFileIn`, :dox:`SimpleIntervalsFileInOut`
* :dox:`UcscFileIn`, :dox:`UcscFileOut`


.. warning::

    Access to compressed files relies on external libraries.
    For instance, you need to have zlib installed for reading ``.gz`` files and libbz2 for reading ``.bz2`` files.
    If you are using Linux or OS X and you followed the :ref:`tutorial-getting-started` tutorial closely, then you should have already installed the necessary libraries.
    On Windows, you will need to follow :ref:`infra-use-install-dependencies` to get the necessary libraries.

    You can check whether you have installed these libraries by running CMake again.
    Simply call ``cmake .`` in your build directory.
    At the end of the output, there will be a section "SeqAn Features".
    If you can read ``ZLIB - FOUND`` and ``BZIP2 - FOUND`` then you can use zlib and libbz2 in your programs.


Basic I/O
---------

This tutorial shows the basic functionalities provided by any class of type :dox:`FormattedFileIn` or :dox:`FormattedFileOut`.
In particular, this tutorial adopts the classes :dox:`BamFileIn` and :dox:`BamFileOut` as concrete types.
The class :dox:`BamFileIn` allows to read files in SAM or BAM format, whereas the class :dox:`BamFileOut` allows to write them.
Nonetheless, **these functionalities are independent from the particular file format** and thus valid for all record-based file formats supported by SeqAn.

The demo application shown here is a simple BAM to SAM converter.

Includes
^^^^^^^^

Support for a specific format comes by including a specific header file.
In this case, we include the BAM header file:

.. includefrags:: demos/tutorial/file_io_overview/example1.cpp
   :fragment: include

Opening and Closing Files
^^^^^^^^^^^^^^^^^^^^^^^^^

Classes of type :dox:`FormattedFileIn` and :dox:`FormattedFileOut` allow to :dox:`FormattedFile#open` and :dox:`FormattedFile#close` files.

A file can be opened by passing the filename to the constructor:

.. includefrags:: demos/tutorial/file_io_overview/example1.cpp
   :fragment: ctor

Alternatively, a file can be opened after construction by calling :dox:`FormattedFile#open`:

.. includefrags:: demos/tutorial/file_io_overview/example1.cpp
   :fragment: open

Note that any file is closed *automatically* whenever the :dox:`FormattedFileIn` or :dox:`FormattedFileOut` object goes out of scope.
Eventually, a file can be closed *manually* by calling :dox:`FormattedFile#close`.

Accessing the Header
^^^^^^^^^^^^^^^^^^^^

To access the header, we need an object representing the format-specific header.
In this case, we use an object of type :dox:`BamHeader`.
The content of this object can be ignored for now, it will be covered in the :ref:`tutorial-io-sam-bam-io` tutorial.

.. includefrags:: demos/tutorial/file_io_overview/example1.cpp
   :fragment: header

The function :dox:`FormattedFileIn#readHeader` reads the header from the input BAM file and :dox:`FormattedFileOut#writeHeader` writes it to the SAM output file.

Accessing the Records
^^^^^^^^^^^^^^^^^^^^^

Again, to access records, we need an object representing format-specific information.
In this case, we use an object of type :dox:`BamAlignmentRecord`.
Each call to :dox:`FormattedFileIn#readRecord` reads one record from the BAM input file and moves the :dox:`BamFileIn` forward.
Each call to :dox:`FormattedFileOut#writeRecord` writes the record just read to the SAM output files.
We check the end of the input file by calling :dox:`FormattedFile#atEnd`.

.. includefrags:: demos/tutorial/file_io_overview/example1.cpp
   :fragment: records

Our small BAM to SAM conversion demo is ready.
The tool still lacks error handling, reading from standard input and writing to standard output.
You are now going to add these features.

Error Handling
--------------

We distinguish between two types of errors: *low-level* file I/O errors and *high-level* file format errors.
Possible file I/O errors can affect both input and output files.
Example of errors are: the file permissions forbid a certain operation, the file does not exist, there is a disk reading error, a file being read gets deleted while we are reading from it, or there is a physical error in the hard disk.
Conversely, file format errors can only affect input files: such errors arise whenever the content of the input file is incorrect or damaged.
Error handling in SeqAn is implemented by means of exceptions.

I/O Errors
^^^^^^^^^^

All :dox:`FormattedFile#FormattedFile FormattedFileIn` and :dox:`FormattedFile#FormattedFile FormattedFileOut` constructors and functions throw exceptions of type :dox:`IOError` to signal *low-level* file I/O errors.
Therefore, it is sufficient to catch these exceptions to handle I/O errors properly.

There is only one exception to this rule.
Function :dox:`FormattedFile#open` returns a ``bool`` to indicate whether the file was opened successfully or not.


Assignment 1
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Improve the program above to detect file I/O errors.

   Hint
     Use the :dox:`IOError` class.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/file_io_overview/solution1.cpp


Format Errors
^^^^^^^^^^^^^

Classes of types :dox:`FormattedFileIn` throw exceptions of type :dox:`ParseError` to signal *high-level* input file format errors.


Assignment 2
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Improve the program above to detect file format errors.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/file_io_overview/solution2.cpp


Streams
-------

The :dox:`FormattedFile#FormattedFile FormattedFileIn` and :dox:`FormattedFile#FormattedFile FormattedFileOut` constructors accept not only filenames, but also standard C++ streams, or any other class implementing the :dox:`StreamConcept Stream` concept.
For instance, you can pass `std::cin` to any :dox:`FormattedFile#FormattedFile FormattedFileIn constructor` and `std::cout` to any :dox:`FormattedFile#FormattedFile FormattedFileOut constructor`.

.. note::

    When writing to `std::cout`, classes of type :dox:`FormattedFileOut` cannot guess the file format from the filename extension.
    Therefore, the file format has to be specified explicitly by providing a tag, e.g. :dox:`FileFormats#Sam` or :dox:`FileFormats#Bam`.

Assignment 3
^^^^^^^^^^^^

.. container:: assignment

   Type
     Application

   Objective
     Improve the program above to write to standard output.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/file_io_overview/solution3.cpp

        Running this program results in the following output.
        
        .. includefrags:: demos/tutorial/file_io_overview/solution3.cpp.stdout


Next Steps
----------

If you want, you can now have a look at the API documentation of the :dox:`FormattedFile` class.

You can now read the tutorials for **already supported file formats**:

* :ref:`tutorial-io-sequence-io`
* :ref:`tutorial-io-sam-bam-io`
* :ref:`tutorial-io-vcf-io`
* :ref:`tutorial-io-bed-io`
* :ref:`tutorial-io-gff-and-gtf-io`

.. COMMENT or, if you want to learn how to develop **support for new file formats** then read the following article:
    * :ref:`tutorial-custom-io`
