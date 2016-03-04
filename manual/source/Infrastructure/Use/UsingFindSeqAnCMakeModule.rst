.. sidebar:: ToC

    .. contents::

.. _infra-use-cmake:

Using SeqAn in CMake-based projects
===================================

Overview
--------

`CMake <http://cmake.org/>`_ is a cross-platform build system generator.
That is, you describe the different executables and binaries and their dependencies ``CMakeLists.txt`` files.
Then, CMake generates build systems from this, for example in the form of Makefiles or Visual Studio projects.

This article will only describe the most basic things about CMake in general and instead focus on how to use SeqAn easily from within CMake projects.
In CMake projects, one uses `modules to find libraries <http://www.vtk.org/Wiki/CMake:How_To_Find_Libraries>`_ such as SeqAn.
SeqAn ships with such a module.

In the following we assume that you have installed CMake on your operating system, if you have not, yet, install it via the operating systems mechanisms (see also :ref:`Setting up SeqAn <infra-use-install>`) and/or `download from the CMake homepage <https://cmake.org/download/>`_.

You should also have a valid C++-Compiler installed, refer to the `GitHub-README <https://github.com/seqan/seqan>`_ to see which compilers are currently supported.

A Running Example
-----------------

Create a folder somewhere, e.g. ``~/devel/my_project`` and in it the following two files:

``my_project.cpp``

.. includefrags:: util/raw_cmake_project/src/main.cpp


``CMakeLists.txt``

.. code-block:: cmake

   # Minimum cmake version
   cmake_minimum_required (VERSION 3.0.0)

   # Name of project and that it is C++ only.
   project (my_project CXX)

   # ----------------------------------------------------------------------------
   # Dependencies
   # ----------------------------------------------------------------------------

   # Search for zlib as a dependency for SeqAn.
   find_package (ZLIB)

   # Load the SeqAn module and fail if not found.
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
   add_executable (my_project my_project.cpp)
   target_link_libraries (my_project ${SEQAN_LIBRARIES})


.. The Project's Contents
.. ^^^^^^^^^^^^^^^^^^^^^^
..
.. The file ``src/main.cpp`` contains a minimal SeqAn program.
..
.. .. includefrags:: util/raw_cmake_project/src/main.cpp
..
.. The root ``CMakeLists.txt`` file just sets up the project name, defines a minimal CMake version, makes all binaries go to the ``bin`` subdirectory, and then includes ``src/CMakeLists.txt``.
..
.. .. includefrags:: util/raw_cmake_project/CMakeLists.txt
..
.. This included file calls ``find_package(SeqAn REQUIRED)``.
.. If the library could not be found, the ``REQUIRED`` parameter will make the ``find_package()`` call fail.
.. Before this, multiple ``find_package()`` calls detect optional dependencies and enable them in the SeqAn library through compiler defines. Note that it is important that these packages be found **before** the SeqAn package is searched.
..
.. .. includefrags:: util/raw_cmake_project/src/CMakeLists.txt
..
.. This is followed by adding the include directory, definitions, and compiler flags required for compiling a program against the SeqAn library,
.. Finally, the source file ``main.cpp`` is compiled into a program called ``main`` and the libraries that SeqAn was configured with are linked to ``main``.
.. Note that SeqAn itself does not require a linking step but when using compression (e.g. for the BAM format), we have to link to ``zlib``.

Building The Project
^^^^^^^^^^^^^^^^^^^^

First you should create a build directory, i.e. for cmake-builds everything happens in a different directory, than in the source. In our case create ``~/devel/my_project-build`` and then in that the ``release`` directory. More on why we use two levels :ref:`here <how-to-recipes-use-parallel-build-directories>`.

.. code-block:: console

   # mkdir -p ~/devel/my_project-build/release
   # cd ~/devel/my_project-build/release

By default, the ``cmake`` program will look for ``FindSeqAn.cmake`` in its module directory.
Usually, this is located in ``/usr/share/cmake/Modules`` or a similar location that is available system-wide.
Depending on how you :ref:`installed SeqAn <infra-use-install>` it will be found by cmake automatically, or not. If not, you have to give the path to cmake via the ``CMAKE_MODULE_PATH`` argument on the command line.

Also, CMake will look for the SeqAn include files in central locations such as ``/usr/local/include``. Again, depending on your install  this will *just work*, but if not, you needto  specify the location via the ``SEQAN_INCLUDE_PATH`` argument.

When using operating system packages of SeqAn it might look like this:

.. code-block:: console

   # cmake ../../my_project

Or, if you did a full git checkout instead, it will look like this:

.. code-block:: console

   # cmake ../../my_project \
       -DCMAKE_MODULE_PATH=~/devel/seqan/util/cmake \
       -DSEQAN_INCLUDE_PATH=~/devel/seqan/include

Finally you can then build the application by calling

* on Makefile-based builds (Linux/Mac/BSD):

    .. code-block:: console

        # make

* Windows

    .. code-block:: console

        # cmake-build (TODO double-check)

You can then run the application in the usual way

* on Makefile-based builds (Linux/Mac/BSD):

    .. code-block:: console

        # ./my_project

* Windows

    .. code-block:: console

        # my_project

.. note:: Changing compilers (Makefile-based builds)

    To use e.g. ``g++-5`` instead of the default ``g++``, add ``-DCMAKE_CXX_COMPILER=g++5`` to your cmake call.

.. note:: Using XCode on Mac

    To use XCode on mac instead of a Makefile-based build, add ``-G Xcode`` to your cmake call.
    TODO explain how to open project file with xcode.

.. note:: Using different Visual Studio versions

    To change the version of Visual Studio you are building against, add ``-G "Visual Studio 10 2010"`` to your cmake call. TODO double-checl that this creates 64bit TODO explain how to open project file with visual studio.


Input / Output of the FindSeqAn Module
--------------------------------------

As with all other modules, you have to make the file ``FindSeqAn.cmake`` available as a CMake module, either by putting it into the same directory as the ``CMakeLists.txt`` that you are using it from or by adding the path to the file ``FindSeqAn.cmake`` to the variable ``CMAKE_MODULE_PATH``.

Then, you can use it as follows (the argument ``REQUIRED`` is optional):

.. code-block:: cmake

    find_package (SeqAn REQUIRED)

Input
^^^^^

SeqAn is somewhat special as a library since it has some optional dependencies.
Certain features in SeqAn can be enabled or disabled, depending on whether the dependencies could be found.

For example:

.. code-block:: cmake

    find_package (ZLIB)
    find_package (BZip2)
    find_package (SeqAn)

If these packages are found **before** SeqAn is searched certain ``SEQAN_HAS_*`` macros are defined and corresponding features become available.

Currently, the following dependencies enable optional features:

``ZLIB``
  zlib compression library

``BZip2``
  libbz2 compression library

``OpenMP``
  OpenMP language extensions to C/C++

If you want ``FindSeqAn.cmake`` to expect the SeqAn build system layout then set the variable ``SEQAN_USE_SEQAN_BUILD_SYSTEM`` to ``TRUE``.
In this case, it will try to locate the library parts from root of the SeqAn source files.

Output
------

The call to ``find_package(SeqAn)`` will set the following variables:

``SEQAN_FOUND``
  Indicate whether SeqAn was found.``

Also the following MACROS are passed to the code indicating whether dependencies were (searched and) found:

``SEQAN_HAS_ZLIB``
  ``TRUE`` `` if zlib was found.``

``SEQAN_HAS_BZIP2``
  ``TRUE`` `` if libbz2 was found.``

``SEQAN_HAS_OPENMP``
  ``TRUE`` `` if OpenMP was found.``

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


