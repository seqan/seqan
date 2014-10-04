.. sidebar:: ToC

   .. contents::


.. _build-manual-using-the-find-seqan-cmake-module:

Using the FindSeqAn CMake Module
--------------------------------

Overview
~~~~~~~~

`CMake <http://cmake.org/>`_ is a cross-platform build system generator.
That is, you describe the different executables and binaries and their dependencies ``CMakeLists.txt`` files.
Then, CMake generates build systems from this, for example in the form of Makefiles or Visual Studio projects.

This article will not describe how to use CMake in general but only how to use SeqAn easily from within CMake projects.
In CMake projects, one uses `modules to find libraries <http://www.vtk.org/Wiki/CMake:How_To_Find_Libraries>`_ such as SeqAn.
SeqAn ships with such a module in ``util/cmake/FindSeqAn.cmake``.

This article describes how to use this module.

Input / Output of the FindSeqAn Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As with all other modules, you have to make the file ``FindSeqAn.cmake`` available as a CMake module, either by putting it into the same directory as the ``CMakeLists.txt`` that you are using it from or by adding the path to the file ``FindSeqAn.cmake`` to the variable ``CMAKE_MODULE_PATH``.

Then, you can use it as follows (the argument ``REQUIRED`` is optional):

.. code-block:: cmake

    find_package (SeqAn REQUIRED)

Input
^^^^^

SeqAn is somewhat special as a library since it has some optional dependencies.
Certain features in SeqAn can be enabled or disabled, depending on whether the dependencies could be found.

You can set the dependencies to search for with the variable ``SEQAN_FIND_DEPENDENCIES`` (which is a list).
For example:

.. code-block:: cmake

    set (SEQAN_FIND_DEPENDENCIES ZLIB BZip2)
    find_package (SeqAn)

Note that ``FindSeqAn.cmake`` itself will not search for its dependencies with the argument ``REQUIRED``. Rather, it will set the variables ``SEQAN_HAS_*`` and add corresponding definitions to ``SEQAN_DEFINIONS`` (see below).

Currently, you can specify the following dependencies:

``ALL``
  Enable all dependencies.

``DEFAULT``
  Enable default dependencies (zlib, OpenMP if available)

``NONE``
  Disable all dependencies.

``ZLIB``
  zlib compression library

``BZip2``
  libbz2 compression library

``OpenMP``
  OpenMP language extensions to C/C++

``CUDA``
  CUDA language extensions to C/C++

If you want ``FindSeqAn.cmake`` to expect the SeqAn build system layout then set the variable ``SEQAN_USE_SEQAN_BUILD_SYSTEM`` to ``TRUE``.
In this case, it will try to locate the library parts from ``core`` and ``extras``.

Output
~~~~~~

The call to ``find_package(SeqAn)`` will set the following variables:

``SEQAN_FOUND``
  Indicate whether SeqAn was found.``

Variables indicating whether dependencies were found:

``SEQAN_HAS_ZLIB``
  ``TRUE`` `` if zlib was found.``

``SEQAN_HAS_BZIP2``
  ``TRUE`` `` if libbz2 was found.``

``SEQAN_HAS_OPENMP``
  ``TRUE`` `` if OpenMP was found.``

``SEQAN_HAS_CUDA``
  ``TRUE`` `` if CUDA was found.``

Variables to be passed to ``include_directories()``, ``target_link_directories()``, and ``add_definitions()`` in your ``CMakeLists.txt``:

``SEQAN_INCLUDE_DIRS``
  A list of include directories.

``SEQAN_LIBRARIES``
  A list of libraries to link against.

``SEQAN_DEFINITIONS``
  A list of definitions to be passted to the compiler.

Required additions to C++ compiler flags are in the following variable:

``SEQAN_CXX_FLAGS``
  C++ Compiler flags to add.

The following variables give the version of the SeqAn library, its major, minor, and the patch version part of the version string.

``SEQAN_VERSION_STRING``
  Concatenated version string, `` ``${SEQAN_VERSION_MAJOR}.${SEQAN_VERSION_MINOR}.${SEQAN_VERSION_PATCH}`` ``.``

``SEQAN_VERSION_MAJOR``
  Major version.

``SEQAN_VERSION_MINOR``
  Minor version.

``SEQAN_VERSION_PATCH``
  Patch-level version.

The following flag defines whether this is a trunk version and the version given by the variables above is meant to be used as the previously released version.

``SEQAN_VERSION_DEVELOPMENT``
  Whether or not this is a pre-release version.

Example
~~~~~~~

Below you can find a minimal example ``CMakeLists.txt`` file that uses the ``FindSeqAn.cmake``.

.. code-block:: cmake

   cmake_minimum_required (VERSION 2.8.2)
   project (core_apps_dfi)

   # ----------------------------------------------------------------------------
   # Dependencies
   # ----------------------------------------------------------------------------

   # Only search for zlib as a dependency for SeqAn.
   set (SEQAN_FIND_DEPENDENCIES ZLIB)
   find_package (SeqAn REQUIRED)

   # ----------------------------------------------------------------------------
   # Build Setup
   # ----------------------------------------------------------------------------

   # Add include directories.
   include_directories (${SEQAN_INCLUDE_DIRS})

   # Add definitions set by find_package (SeqAn).
   add_definitions (${SEQAN_DEFINITIONS})

   # Add CXX flags found by find_package (SeqAn).
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

   # Add executable and link against SeqAn dependencies.
   add_executable (app app.cpp)
   target_link_libraries (dfi ${SEQAN_LIBRARIES})
