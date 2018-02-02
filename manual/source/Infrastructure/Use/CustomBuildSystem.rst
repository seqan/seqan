.. sidebar:: ToC

   .. contents::


.. _infra-use-custom:

Integration with your own Build System
======================================

This page gives an example of how to use SeqAn in your application based on your own Makefiles.
You should be able to adapt the descriptions to configure your build system and/or IDE.

.. tip::

   **SeqAn is a header-only library.**

   Simply adding its include folder to your include path or installing it globally makes it available to your program.

C++14 Standard
--------------

On GNU/Linux, BSD and macOS, always add ``-std=c++14`` (or a newer standard) when building on the command line.

For XCode on macOS you need to set this option in the project settings.

As of Visual Studio 2015 our subset of C++14 is already available in all supported compilers.

OpenMP
------

On GNU/Linux, BSD and macOS, add ``-fopenmp`` unless you are using Clang versions older than 3.8.0.

For XCode on macOS OpenMP is not yet available.

With Visual Studio OpenMP is switched on by default.

Compiler specifics
------------------

Intel Compiler
^^^^^^^^^^^^^^

The Intel Compiler does not ship a c++ standard library on its own and will use
the one pre-installed on the system (e.g., the one from g++). This can be a
problem [especially for cluster users through the use of a module system], if
the standard library by a default g++ installation is to old.

Please check with the following command which g++ version is being used and make
sure it matches the supported gcc versions.

.. code-block:: console

    # icpc -v
    icpc version 17.0.2 (gcc version 5.4.0 compatibility)

If you have multiple g++ installations, you can choose the standard library by
``icpc -gxx-name=g++-5.4.0 -gcc-name=gcc-5.4.0 …``.

You may have to add the path of the library to ```$LD_LIBRARY_PATH`` for the
linker.

Operating System specifics
--------------------------

GNU/Linux
^^^^^^^^^

**Libraries**

Add ``-lrt -lpthread`` to the compiler call.

Note static linking against pthread might cause issues on some linux distributions, especially with gcc-4.9.
In this case you need to explicitly link against the whole archive like: ``-Wl,--whole-archive -lpthread -Wl,--no-whole-archive``.
Or just use a more recent version of gcc!

BSD
^^^

**Libraries**

Add ``-lpthread -lexecinfo -lelf`` to the compiler call.

**Misc**

Also define ``-D_GLIBCXX_USE_C99=1`` if you are using gcc-4.9.

Warning levels
--------------

It is recommended to compile your programs with as many warnings enabled as possible.
This section explains which flags to set for different compilers.

GCC, Clang, ICC (unix)
^^^^^^^^^^^^^^^^^^^^^^

The following flags are recommended:

``-W -Wall -pedantic``
  Maximal sensitivity of compiler against possible problems.


Visual Studio
^^^^^^^^^^^^^

For Visual Studio, the following flags are recommended:

::

    /W2 /wd4996 -D_CRT_SECURE_NO_WARNINGS

Explanation:

``/W2``
  Warning level 2 is pretty verbose already.
  In the future, we will support level 3 without warnings in SeqAn code.

``/wd4996``
  Allows the use of some deprecated functions without warnings.

``-D_CRT_SECURE_NO_WARNINGS`` ::``
   Some C functions like ``sprintf`` are prone to incorrect usage and security holes.
   Replacing such calls does not have a high priority right now since SeqAn is usually not used on servers facing the outside world.

Preprocessor Defines Affecting SeqAn
------------------------------------

There are certain preprocessor symbols that affect the behaviour of SeqAn.

SEQAN_ENABLE_DEBUG
^^^^^^^^^^^^^^^^^^

possible value
  0, 1

default
  0

meaning
  If set to 1, assertions within SeqAn (``SEQAN_ASSERT...``) are enabled, they are disabled otherwise.
  Is forced to 1 if ``SEQAN_ENABLE_TESTING`` is true.
  This flag will internally always correspond to the inverse of ``NDEBUG``, i.e. setting it to 1 will force ``NDEBUG`` to be undefined and setting it to 0 will forcefully set ``NDEBUG``.

SEQAN_ENABLE_TESTING
^^^^^^^^^^^^^^^^^^^^

possible value
  0, 1

default
  0

meaning
 This makes the code very slow, and should only be used when running the SeqAn unit tests.
 Has to be set to 1 for tests to work.

SEQAN_ASYNC_IO
^^^^^^^^^^^^^^

possible value
  0, 1

default
  0 on ``FreeBSD/32Bit`` and ``OpenBSD/*``; 1 otherwise

meaning
 Whether asynchronous input/output is available.

SEQAN_HAS_EXECINFO
^^^^^^^^^^^^^^^^^^

possible value
  0, 1

default
  depends on platform / existance of ``<execinfo.h>``

meaning
 This should almost always be set to 1 on non-Windows platforms!

SEQAN_HAS_BZIP2
^^^^^^^^^^^^^^^

possible value
  0, 1

default
  0

meaning
 If set to 1 then libbzip2 is expected to be available.
 You have to link against the library (e.g. add ``-lbz2`` to your linke rflags) and ``bzlib.h`` must be in your include path.

SEQAN_HAS_ZLIB
^^^^^^^^^^^^^^

possible value
  0, 1

default
  0

meaning
 If set to 1 then zlib is expected to be available.
 You have to link against the library (e.g. add ``-lz`` to your linker flags) and ``zlib.h`` must be in your include path.

SEQAN_HAS_OPENMP
^^^^^^^^^^^^^^^^

possible value
  0, 1

default
  0

meaning
 If set to 1 then OpenMP is expected to be available.
 You might have to add ``-fopenmp`` and possibly ``-lgomp`` to your build. And OpenMP needs to be supported by your compiler.

SEQAN_VERSION_CHECK_OPT_IN
^^^^^^^^^^^^^^^^^^^^^^^^^^

meaning 
 If set then the version update feature is enabled but deactivated in the argument parser options and must be explicitly activated by the user.
 If not specified the update feature is enabled and activated by default.

usage
 Add compiler flag: ``-DSEQAN_VERSION_CHECK_OPT_IN``

SEQAN_DISABLE_VERSION_CHECK
^^^^^^^^^^^^^^^^^^^^^^^^^^^

meaning
 If set the version update feature is disabled in the argument parser.

usage
 add compiler flag: ``-DSEQAN_DISABLE_VERSION_CHECK`` 

SEQAN_BGZF_NUM_THREADS
^^^^^^^^^^^^^^^^^^^^^^^^^^^

possible value
 positive integer

default
 16

meaning
 Number of threads to use for BGZF I/O.

usage
 add compiler flag: ``-DSEQAN_BGZF_NUM_THREADS=value`` 

Settings Projects Using Seqan
-----------------------------

You normally want to have at least two build modes: one for debugging and one for optimized compiling.
The following settings have to be applied to your IDE project/Makefiles (below is an example for a Makefile based project).

Debug Builds
^^^^^^^^^^^^

Besides enabling debug symbols and disabling optimization, there are the
following SeqAn specific settings to be applied.

- Add SeqAn to your include path
- Define ``SEQAN_ENABLE_DEBUG`` to be ``1``.

This translates into the following GCC flags:

::

    -g -O0 -DSEQAN_ENABLE_DEBUG=1 -I${PATH_TO_SEQAN_INSTALL}/include

Release/Optimized Builds
^^^^^^^^^^^^^^^^^^^^^^^^

Besides disabling debug symbols, enabling optimization and disabling assertions in the standard library, there are the following SeqAn specific settings to be applied.

- Add SeqAn to your include path
- Define ``NDEBUG``. This will make sure that ``SEQAN_ENABLE_DEBUG`` is 0 and also other STL includes of your program are not slowed down.

This translates into the following GCC flags:

::

    -O3 -DNDEBUG -I${PATH_TO_SEQAN_INSTALL}/include

.. caution::

    While some guides tell you to not use ``-O3`` this is absolutely crucial for SeqAn based applications to perform well. Unoptimized builds are slower by multiple factors!

An Example Project Based On Makefiles
-------------------------------------

We will create a project with good old Makefiles and GCC.
The program will not do much but can serve as a minimal example on how to use SeqAn with your own build process.
You should be able to adapt this guide to your favourite build system or IDE.

The example project can be found in ``util/makefile_project``.
The project layout looks like this:

::

    .
    |-- Makefile.rules
    |-- Makefile
    |-- README
    |-- debug
    |   `-- Makefile
    |-- release
    |   `-- Makefile
    `-- src
        `-- main.cpp

main.cpp
^^^^^^^^

We have one directory ``src`` for source files.
The file ``main.cpp`` looks as follows:

.. includefrags:: util/makefile_project/src/main.cpp
   :language: cpp

It includes SeqAn headers just as you would within the SeqAn CMake framework.

Now, consider the contents of the Makefiles:

Makefile.rules
^^^^^^^^^^^^^^

Contains the necessary commands to build the object file for the program ``main.cpp`` and then make an executeable ``main`` from it and clean targets.
This file is included from the files ``release/Makefile`` and ``debug/Makefile``.

.. includefrags:: util/makefile_project/Makefile.rules
   :language: make

Makefile
^^^^^^^^

Allows to build both debug and release builds by calling ``make debug``, ``make release`` or ``make all`` from the project directory.
Removes all binaries with ``make clean``.

.. includefrags:: util/makefile_project/Makefile
   :language: make

debug/Makefile, release/Makefile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The file ``debug/Makefile`` looks as follows.

.. includefrags:: util/makefile_project/debug/Makefile
   :language: make

The file ``release/Makefile`` looks as follows.

.. includefrags:: util/makefile_project/release/Makefile
   :language: make

These Makefiles include the file ``Makefile.rules``.
They add build type specific arguments to the variables ``$(CXXFLAGS)``.
For debug builds, debug symbols are enabled, optimization level 0 is chosen, testing is enabled in SeqAn and debugging is disabled.
For release builds, debug symbols are not, optimization level 3 is chosen, testing and debugging are both disabled in SeqAn.
For good measure, we also disable assertions in the C library with ``-DNDEBUG``.

Notes
^^^^^

Above we added the include path to SeqAn's include directory manually.
By changing the include path, we can install the SeqAn library anywhere.
For example, we could create a directory ``include`` parallel to ``src``, copy the release version of SeqAn into it and then change the include path of the compiler to point to this directory (value ``../include``).

Short Version
-------------

+---------+---------------------+---------------------------------------------------------------------------------------------------------------------------------------------+
| OS      | Compiler            | Flags                                                                                                                                       |
+=========+=====================+=============================================================================================================================================+
| Linux   | GCC/Clang≥3.8/ICC   | ``-I /path/to/seqan/include -std=c++14 -O3 -DNDEBUG -W -Wall -pedantic -fopenmp -lpthread -lrt``                                            |
+---------+---------------------+---------------------------------------------------------------------------------------------------------------------------------------------+
| BSD     | GCC/Clang≥3.8/ICC   | ``-I /path/to/seqan/include -std=c++14 -O3 -DNDEBUG -W -Wall -pedantic -fopenmp -lpthread -lexecinfo -lelf -D_GLIBCXX_USE_C99=1``           |
+---------+---------------------+---------------------------------------------------------------------------------------------------------------------------------------------+
| macOS   | system's Clang      | ``-I /path/to/seqan/include -std=c++14 -O3 -DNDEBUG -W -Wall -pedantic``                                                                    |
+---------+---------------------+---------------------------------------------------------------------------------------------------------------------------------------------+
| Windows | Visual Studio MSVC  | ``/W2 /wd4996 -D_CRT_SECURE_NO_WARNINGS``                                                                                                   |
+---------+---------------------+---------------------------------------------------------------------------------------------------------------------------------------------+

Adapt the include path to the actual place of SeqAn's ``include`` folder!
