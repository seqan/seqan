.. sidebar:: ToC

   .. contents::


.. _build-manual-integration-with-your-own-build-system:

Integration with your own Build System
--------------------------------------

The CMake build system that SeqAn ships with is meant for people who want to build the applications, tests, and demos that SeqAn ships with.
It has the advantage that new such programs have only to be added by a certain convention and they get added to the Makefiles/project files on the next ``cmake`` call.
If you just want to use SeqAn in your own project, it might not be a good fit.
One of the disadvantages is that CMake will overwrite your project files on every call.
Another disadvantage is that the generated project files are huge and might take a long while to load.

This page gives an example of how to use SeqAn in your application based on your own Makefiles.
You should be able to adapt the descriptions to configuring your build system and/or IDE.

.. tip::

   SeqAn is a header library only.
   Simply add ``core/include`` and ``extras/include`` to your include path and you can use SeqAn, as seen in the `Short Version`_.
   See below how to enable using zlib for BAM access, for example.

Libraries on Linux
~~~~~~~~~~~~~~~~~~

On Linux, you have to link against ``librt``.
For GCC, add the flag ``-lrt`` to the ``g++`` compiler call.

Compiler Flags
~~~~~~~~~~~~~~

It is recommended to compile your programs with as many warnings enabled as possible.
This section explains which flags to set for different compilers.

GCC
^^^

For GCC, the following flags are recommended:

::

    -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros

Explanation:

``-W -Wall -pedantic``
  Maximal sensitivity of compiler against possible problems.

``-Wno-variadic-macros``
  The assertion macros are variadic.
  Variadic macros were standardized in C99 but are not part of C++98 so GCC warns against their usage.
  Disable these warnings.

``-Wno-long-long``
  64 bit integers (``long long``) are not supported in C++98, but GCC implements them nevertheless but warns against their usage in pedantic mode.
  We really want 64 bit integers, though.

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are certain preprocessor symbols that affect the behaviour of SeqAn.

SEQAN_ENABLE_DEBUG
^^^^^^^^^^^^^^^^^^

possible value
  0, 1

default
  1

meaning
  If set to 1, assertions within SeqAn (``SEQAN_ASSERT...``) are enabled, they are disabled otherwise.
  Is forced to 1 if ``SEQAN_ENABLE_TESTING`` is true.
  If not set, is set to 0 if ``NDEBUG`` is defined and set to 1 if undefind and ``NDEBUG`` is not defined.

SEQAN_ENABLE_TESTING
^^^^^^^^^^^^^^^^^^^^^^

possible value
  0, 1

default
  0

meaning
 If set to 1, checkpoints are enabled.
 This makes the code very slow, however, and should only be used when running the tests.
 Has to be set to 1 for tests to work.

SEQAN_HAS_BZIP2
^^^^^^^^^^^^^^^

possible value
  0, 1

default
  0

meaning
 If set to 1 then libbzip2 is available.``
 You have to link against the library (e.g. add ``-lbz2`` to your linke rflags) and ``bzlib.h`` must be in your include path.

SEQAN_HAS_ZLIB
^^^^^^^^^^^^^^

possible value
  0, 1

default
  0

meaning
 If set to 1 then zlib is available.
 You have to link against the library (e.g. add ``-lz`` to your linker flags) and ``zlib.h`` must be in your include path.

Settings Projects Using Seqan
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You normally want to have at least two build modes: one for debugging and one for optimized compiling.
The following settings have to be applied to your IDE project/Makefiles (below is an example for a Makefile based project).

Debug Builds
^^^^^^^^^^^^

Besides enabling debug symbols and disabling optimization, there are the
following SeqAn specific settings to be applied.

- Add the path to the directory ``seqan`` to your include path.
- Define ``SEQAN_ENABLE_DEBUG`` to be ``1``.
  Alternatively, you can leave ``SEQAN_ENABLE_DEBUG`` undefined and not define ``NDEBUG``.
- Define ``SEQAN_ENABLE_TESTING`` to be ``0``.

This translates into the following GCC flags:

::

    -g -O0 -DSEQAN_ENABLE_TESTING=0 -I${PATH_TO_CORE}/include \
      -I${PATH_TO_EXTRAS}/include

Release/Optimized Builds
^^^^^^^^^^^^^^^^^^^^^^^^

Besides disabling debug symbols, enabling optimization and disabling assertions in the standard library, there are the following SeqAn specific settings to be applied.

* Add the path to the directory ``seqan`` to your include path.
* Define ``NDEBUG``.
  This will make ``SEQAN_ENABLE_DEBUG`` be defined as ``0`` if you don't defined ``SEQAN_ENABLE_DEBUG`` otherwise.
* Define ``SEQAN_ENABLE_TESTING`` to be ``0``.

This translates into the following GCC flags:

::

    -O3 -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -I${PATH_TO_CORE}/include \
      -I${PATH_TO_EXTRAS}/include

An Example Project Based On Makefiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We will create a project with good old Makefiles and GCC.
The program will not do much but can serve as a minimal example on how to use SeqAn with your own build process.
You should be able to adapt this guide to your favourite build system or IDE.

The example project can be found in ``misc/makefile_project``.
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

.. includefrags:: misc/makefile_project/src/main.cpp
   :language: cpp

It includes SeqAn headers just as you would within the SeqAn CMake framework.

Now, consider the contents of the Makefiles:

Makefile.rules
^^^^^^^^^^^^^^

Contains the necessary commands to build the object file for the program ``main.cpp`` and then make an executeable ``main`` from it and clean targets.
This file is included from the files ``release/Makefile`` and ``debug/Makefile``.

.. includefrags:: misc/makefile_project/Makefile.rules
   :language: make

Makefile
^^^^^^^^

Allows to build both debug and release builds by calling ``make debug``, ``make release`` or ``make all`` from the project directory.
Removes all binaries with ``make clean``.

.. includefrags:: misc/makefile_project/Makefile
   :language: make

debug/Makefile, release/Makefile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The file ``debug/Makefile`` looks as follows.

.. includefrags:: misc/makefile_project/debug/Makefile
   :language: make

The file ``release/Makefile`` looks as follows.

.. includefrags:: misc/makefile_project/release/Makefile
   :language: make

These Makefiles include the file ``Makefile.rules``.
They add build type specific arguments to the variables ``$(CXXFLAGS)``.
For debug builds, debug symbols are enabled, optimization level 0 is chosen, testing is enabled in SeqAn and debugging is disabled.
For release builds, debug symbols are not, optimization level 3 is chosen, testing and debugging are both disabled in SeqAn.
For good measure, we also disable assertions in the C library with ``-DNDEBUG``.

Notes
^^^^^

Note we that added include path to the directory ``include`` that contains the directory ``seqan``.
By changing the include path, we can install the SeqAn library anywhere.
For example, we could create a directory ``include`` parallel to ``src``, copy the release version of SeqAn into it and then change the include path of the compiler to point to this directory (value ``../include``).

Short Version
~~~~~~~~~~~~~

* Add both ``core/include`` and ``extras/include`` to your include path (``-I``).
* Linux/GCC flags: ``-lrt`` (required) ``-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros`` (optional).
* Windows/MSVC flags: ``/W2 /wd4996 -D_CRT_SECURE_NO_WARNINGS`` (optional).
* Defines: ``NDEBUG`` to also disable SeqAn assertions in release mode.
