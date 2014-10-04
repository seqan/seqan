.. sidebar:: ToC

   .. contents::


.. _tutorial-file-io:

File I/O
========

Learning Objective
  In this tutorial, you will learn about the new file I/O infrastructure in SeqAn.
  You will get an overview of the different layers in the library, an introduction on the :dox:`StreamConcept` concept, the :dox:`Stream` class, and :dox:`MMapString MMap-Strings`.

Difficulty
  Advanced

Duration
  60 min

Prerequisites
  :ref:`tutorial-input-output-overview`, :ref:`tutorial-indexed-fasta-io`, :ref:`tutorial-basic-sam-bam-io`

This tutorial introduces the low-level facilities of file I/O in SeqAn:

* There is a concept called :dox:`StreamConcept` in the SeqAn library that stream data types have to implement.
  There also is the class :dox:`Stream` that provides implementations of the concept together with its specializations.
  (If you want to provide your own Stream implementation, you should specialize the class :dox:`Stream`).
* Particularly, there are the specializations :dox:`GzFileStream` and :dox:`BZ2FileStream BZ2 FileStream` that provide access to compressed files.
* Furthermore, SeqAn allows to access memory mapped files using the :dox:`MMapString MMap String` specialization.

The target audience consists of developers (1) who want to learn how to use memory mapped files and compressed streams, or (2) who want to have raw, byte-wise read and write access to files, or (3) who want to get a deeper understanding of the I/O system in the SeqAn library.

Note that this tutorial has more of a presentational character with fewer tasks.

Streams
-------

The :ref:`tutorial-input-output-overview` tutorial has already given you a good overview of streams in SeqAn and how to open them for reading and writing.
As a reminder: Always open your streams in binary mode to circument problems with getting and setting positions within files on Windows.
How exactly you can open files in binary mode depends on the library you are using. Consult the documentation of the library you are using for I/O.

The Stream Concept
^^^^^^^^^^^^^^^^^^

The stream concept requires the following functions which work on already files (e.g. ``FILE *``, ``std::fstream``, or :dox:`Stream` objects).

+--------------------------------------------------+---------------------------------------------------------------------------+
| Function                                         | Summary                                                                   |
+==================================================+===========================================================================+
| :dox:`StreamConcept#streamEof`                   | Return whether stream is at end of file.                                  |
+--------------------------------------------------+---------------------------------------------------------------------------+
| :dox:`StreamConcept#streamError`                 | Return error code of stream.                                              |
+--------------------------------------------------+---------------------------------------------------------------------------+
| :dox:`StreamConcept#streamFlush`                 | Flush stream buffer.                                                      |
+--------------------------------------------------+---------------------------------------------------------------------------+
| :dox:`StreamConcept#streamPeek`                  | Get next character from stream without changing the position in the file. |
+--------------------------------------------------+---------------------------------------------------------------------------+
| :dox:`StreamConcept#streamPut`                   | Write a value to the output, converted to string.                         |
+--------------------------------------------------+---------------------------------------------------------------------------+
| :dox:`StreamConcept#streamReadBlock streamBlock` | Read a block of ``char`` values from the stream.                          |
+--------------------------------------------------+---------------------------------------------------------------------------+
| :dox:`StreamConcept#streamReadChar`              | Read one character from the stream.                                       |
+--------------------------------------------------+---------------------------------------------------------------------------+
| :dox:`StreamConcept#streamSeek`                  | Set stream's location.                                                    |
+--------------------------------------------------+---------------------------------------------------------------------------+
| :dox:`StreamConcept#streamTell`                  | Retrieve stream's location.                                               |
+--------------------------------------------------+---------------------------------------------------------------------------+
| :dox:`StreamConcept#streamWriteBlock`            | Write an array of ``char`` to the stream.                                 |
+--------------------------------------------------+---------------------------------------------------------------------------+

Not all functions might be available for all streams.
The metafunction :dox:`StreamConcept#HasStreamFeature` provides information about the stream types.

Stream Adaptions
^^^^^^^^^^^^^^^^

The following C/C++ I/O interfaces can be adapted to the :dox:`StreamConcept` concept.

+-------------------------------------------------------------------------+---------------------------------------+
| File Type                                                               | Description                           |
+=========================================================================+=======================================+
| ``FILE*``                                                               | C standard library files.             |
+-------------------------------------------------------------------------+---------------------------------------+
| ``std::fstream``, ``std::ifstream``, ``std::ofstream``                  | C++ iostream library file streams     |
+-------------------------------------------------------------------------+---------------------------------------+
| ``std::stringstream``, ``std::istringstream``, ``std::ostringstream``   | C++ iostream library string streams   |
+-------------------------------------------------------------------------+---------------------------------------+

This way, we can use the common C++ I/O types through a common interface.
Also, we could add adaptions of other file and stream data types to the :dox:`StreamConcept` concept.

The following example shows how to use the :dox:`StreamConcept` global interface functions to copy the contents of the file ``in.txt`` to the file ``out.txt``.

.. includefrags:: extras/demos/tutorial/file_io/example1.cpp

Assignment 1
""""""""""""

.. container:: assignment

   Reading / Writing

   Type
     Review

   Objective
     Write a program that accepts three parameters from the command line.
     The first one should identify the stream type to use (e.g. ``"file"`` for ``FILE*`` and ``"fstream"`` for ``std::fstream``).
     The second should be either ``'r'`` or '``w'`` for reading/writing.
     The third one should be a file name.
     The program should, depending on the parameters, open the given file name in read/write mode using the given file type.
     When reading, it should display the file contents on stdout.
     When writing, it should put the string ``"Hello world!\n"`` into the file.

   Hint
     You can encapsulate the reading and writing in their own function templates.
     This allows you to remove redundancy from the code.

   Solution ::
    .. container:: foldable

       .. includefrags:: extras/demos/tutorial/file_io/solution1.cpp

Char Arrays As Streams
^^^^^^^^^^^^^^^^^^^^^^

Sometimes it is useful to treat variables of type ``char *`` or ``char[]`` as streams, e.g., for parsing.
You can use the :dox:`CharArrayStream Char-Array Stream` specialization for this purpose.

.. code-block:: cpp

   char const * str = "me, myself and my pony";
   seqan::Stream<seqan::CharArray<char const *> > wrapper(str, str + strlen(str));
   // We can now read from wrapper as if it was a stream.

Compressed Streams
^^^^^^^^^^^^^^^^^^

For accessing ``.gz`` and ``.bz2`` files, the ``stream`` module contains specializations of the class :dox:`Stream`.
The main reason for being :dox:`Stream` specializations instead of adaptions is that zlib and bzlib use too generic data types, e.g., ``void*``, where global functions might have unwanted side effects.

Use the following :dox:`Stream` specializations to read and write zlib and bzlib compressed files.

+---------------------------------------+---------------------------------------------------------------------------+
| Stream Class                          | Description                                                               |
+=======================================+===========================================================================+
| :dox:`GzFileStream GZ File Stream`    | Wraps the `zlib <http://zlib.org>`_ functionality for ``.gz`` files.      |
+---------------------------------------+---------------------------------------------------------------------------+
| :dox:`BZ2FileStream BZ2 File Stream`  | Wraps the `bzlib <http://bzlib.net>`_ functionality for ``.bz2`` files.   |
+---------------------------------------+---------------------------------------------------------------------------+

zlib files have a decent compression ratio and support quite fast compression and decompression.
bz2 files are fairly slow to read and write, although the compression ratio is better.
For most bioinformatics applications, you will prefer zlib over bzlib.

If you are using SeqAn's build system, *zlib* and *libbz2* will be detected automatically.
On Linux and Mac Os X, these libraries are usually already installed.
If you are using Windows, then you can follow the instructions in :ref:`how-to-install-contribs-on-windows` for installing the libraries.
If you are using your own build system, see BuildManual/IntegrationWithYourOwnBuildSystem for the necessary configuration steps.

Both specializations can be constructed with an already open underlying compressed stream, e.g. you can pass the ``gzFile``/``BZFILE*``, that you want to work on, to the stream.
They are meant as very thin wrappers around the handle for the compressed stream.
This has the advantage that you have full access to the compression settings etc. and the wrappers only add error flags and so on when necessary.
For more convenience, you can also use the :dox:`File#open` function to open them.

The following example shows (1) how to conditionally enable zlib and bzlib support, (2) how to open ``gzFile`` and ``BZFILE*`` handles for reading and their corresponding wrappers and (3) the possibilities for error checking.

In the header of the program, we include the zlib and bzlib headers if the correct preprocessor symbols are set.
Also, we'll include the required SeqAn headers.

.. includefrags:: extras/demos/tutorial/file_io/stream_compression_formats.cpp
   :fragment: header

The first routine demonstrates how to open a ``.gz`` file and write its contents to stdout with full error handling.
Note that writing char-by-char is probably not the best idea in a real-world program.

.. includefrags:: extras/demos/tutorial/file_io/stream_compression_formats.cpp
   :fragment: open-gz

The next routine demonstrates how to open a ``.bz2`` file and write its contents to stdout, again with full error handling.

.. includefrags:: extras/demos/tutorial/file_io/stream_compression_formats.cpp
   :fragment: open-bz2

And finally, the code that calls the functions from above.

.. includefrags:: extras/demos/tutorial/file_io/stream_compression_formats.cpp
   :fragment: main

Now, let's test the program.
We'll first create gzip and bzip2 compressed text files and an uncompressed text file.
Then, we'll run our demo program on these files.
Note that the :dox:`BZ2FileStream` fails when reading from the file, not when opening the file.

.. code-block:: console

   # echo 'foo' > test.txt
   # gzip test.txt
   # echo 'bar' > test.txt
   # bzip2 test.txt
   # echo 'bz' > test.txt
   # ./extras/demos/tutorial/stream/tutorial_stream_compression_formats test.txt
   ERROR: GZip file has the wrong format!
   ERROR: Reading byte from BZ2 file.
   # ./extras/demos/tutorial/stream/tutorial_stream_compression_formats test.txt.gz
   foo
   ERROR: Reading byte from BZ2 file.
   # ./extras/demos/tutorial/stream/tutorial_stream_compression_formats test.txt.bz2
   ERROR: GZip file has the wrong format!
   bar

Assignment 2
""""""""""""

.. container:: assignment

   Writing a File Compression/Decompression Tool

   Type
     Application

   Objective
     Write a file compression/decompression tool.
     The first argument should be the format to read/write, e.g. ``"gz"`` for gzip and ``"bz2"`` for bzip2.
     The second argument should be the direction, i.e. "c" for "compress", "x" for "extract".
     The third and fourth arguments should be the source/target files.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/file_io/solution2.cpp

Memory Mapped Files
-------------------

Memory mapped files allow very fast access to files since they enable you to read data with few, if any additional buffers.
Wikipedia has a `nice article on memory mapped files <http://en.wikipedia.org/wiki/Memory-mapped_file>`_.

In SeqAn, you access memory mapped files using the :dox:`MMapString` specialization.
After opening the mapped string using :dox:`File#open`, you can access its contents as if you were manipulating a normal :dox:`String`.
The following shows a simple example:

.. includefrags:: extras/demos/tutorial/file_io/mmap_string_example.cpp

An example execution of the program:

.. code-block:: console

   # echo 'foo' > test.txt
   # ./extras/demos/tutorial/stream/tutorial_mmap_string_example test.txt
   This is the first mapped string!
   foo

Next Steps
----------

* Read `Wikipedia's article on memory mapped files <http://en.wikipedia.org/wiki/Memory-mapped_file>`_.
* Read the :ref:`tutorial-lexical-casting` tutorial to learn how to read text from files that represent numbers (e.g. ``"100"``) into values of numeric types such as ``int``.
* Read the :ref:`tutorial-parsing` tutorial to learn how to write parsers for your own file formats.
* Continue with the :ref:`tutorial`.
