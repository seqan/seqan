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

A Running Example
~~~~~~~~~~~~~~~~~

In the directory ``util/raw_cmake_project``, you can find a small example project that uses the ``FindSeqAn.cmake`` module outside the SeqAn build system.
The directory sturcture looks as follows:

.. code-block:: text

   .
   |-- CMakeLists.txt
   |-- README
   `-- src
       |-- CMakeLists.txt
       `-- main.cpp

The Project's Contents
^^^^^^^^^^^^^^^^^^^^^^

The file ``src/main.cpp`` contains a minimal SeqAn program.

.. includefrags:: util/raw_cmake_project/src/main.cpp

The root ``CMakeLists.txt`` file just sets up the project name, defines a minimal CMake version, makes all binaries go to the ``bin`` subdirectory, and then includes ``src/CMakeLists.txt``.

.. includefrags:: util/raw_cmake_project/CMakeLists.txt

This included file calls ``find_package(SeqAn REQUIRED)``.
If the library could not be found, the ``REQUIRED`` parameter will make the ``find_package()`` call fail.
Before this, the variable ``SEQAN_FIND_DEPENDENCIES`` is set such that zlib and libbz2 are searched for the in ``find_package()`` call and enabled in the SeqAn library through compiler defines.

.. includefrags:: util/raw_cmake_project/src/CMakeLists.txt

This is followed by adding the include directory, definitions, and compiler flags required for compiling a program against the SeqAn library,
Finally, the source file ``main.cpp`` is compiled into a program called ``main`` and the libraries that SeqAn was configured with are linked to ``main``.
Note that SeqAn itself does not require a linking step but when using compression (e.g. for the BAM format), we have to link to ``zlib``.

Building The Project
^^^^^^^^^^^^^^^^^^^^

By default, the ``cmake`` program will look for ``FindSeqAn.cmake`` in its module directory.
Usually, this is located in ``/usr/share/cmake-2.8/Modules`` or a similar location that is available system-wide.
Installing ``FindSeqAn.cmake`` there is one option of making it available in your ``CMakeLists.txt``.
A better option might be to pass the path to the ``util/cmake`` directory of your SeqAn checkout to the ``CMAKE_MODULE_PATH`` CMake variable through the command line.

Also, CMake will look for the SeqAn include files in central location such as ``/usr/local/include``.
Instead of installing SeqAn there, you can pass additional directories to search in using the CMake variable ``SEQAN_INCLUDE_PATH``.

Putting this together, you can execute ``cmake`` for the example CMake project with the following command line:

.. code-block:: console

   # mkdir -p ~/tmp/cmake_example_build
   # cd ~/tmp/cmake_example_build
   # cmake path/to/raw_cmake_project \
       -DCMAKE_MODULE_PATH=~/seqan_checkout/util/cmake \
       -DSEQAN_INCLUDE_PATH=~/seqan_checkout/include
   [...]
   # make main && ./bin/main
   Hello SeqAn!

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
In this case, it will try to locate the library parts from root of the SeqAn source files.

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
   project (apps_dfi)

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
