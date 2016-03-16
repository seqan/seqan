.. sidebar:: ToC

    .. contents::

.. _infra-use-cmake:

Using SeqAn in CMake-based projects
===================================

Overview
--------

`CMake <http://cmake.org/>`_ is a cross-platform build system generator where you describe different executables, binaries and their dependencies in ``CMakeLists.txt`` files.
Then, CMake generates build systems such as Makefiles or Visual Studio projects from these files. This article describes only the most basic things about CMake in general and focuses on how to use SeqAn easily from within CMake projects.

In CMake projects, one uses `modules to find libraries <http://www.vtk.org/Wiki/CMake:How_To_Find_Libraries>`_ such as SeqAn.
SeqAn ships with such a module.

In the following we assume that you have installed CMake on your operating system. If you have not yet, install it via the operating systems mechanisms (see also :ref:`Setting up SeqAn <infra-use-install>`) and/or `download from the CMake homepage <https://cmake.org/download/>`_.

You should also have a valid C++-Compiler installed. Refer to the `GitHub-README <https://github.com/seqan/seqan>`_ to see which compilers are currently supported.

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

Building The Project
^^^^^^^^^^^^^^^^^^^^

First you should create a build directory, i.e. for cmake-builds everything happens in a different directory, than in the source directory. In our case create the directory ``~/devel/my_project-build`` and in there a folder ``release``. More on why we use two levels :ref:`here <infra-use-cmake-build-dirs>`.

.. code-block:: console

   # mkdir -p ~/devel/my_project-build/release
   # cd ~/devel/my_project-build/release

By default, the ``cmake`` program will look for ``FindSeqAn.cmake`` in its module directory.
Usually, this is located in ``/usr/share/cmake/Modules`` or a similar location that is available system-wide.
Depending on how you :ref:`installed SeqAn <infra-use-install>` it might be found by cmake automatically. If not, you have to give the path to cmake via the ``CMAKE_MODULE_PATH`` argument on the command line.

Also, CMake will look for the SeqAn include files in central locations such as ``/usr/local/include``. Again, depending on your installation this might *just work*. If not, you need to specify the location via the ``SEQAN_INCLUDE_PATH`` argument.

When using operating system packages of SeqAn and the default compiler it might look like this:

.. code-block:: console

   # cmake ../../my_project

If you instead did a full git checkout to your home-directory in the previous step, it might look like this:

.. code-block:: console

   # cmake ../../my_project \
       -DCMAKE_MODULE_PATH=~/devel/seqan/util/cmake \
       -DSEQAN_INCLUDE_PATH=~/devel/seqan/include

.. tip::

    Depending on your setup you might need to manually choose a more modern compiler and/or activate C++11 support! Please read :ref:`this page <infra-use-cmake-build-dirs>` for more information on configuring CMake builds.

Finally you can then build the application by calling

* on Makefile-based builds (Linux/Mac/BSD):

    .. code-block:: console

        # make

* Windows

    .. code-block:: console

        # cmake --build .

**The above step is the only step you need to repeat when changing your source code.** You only have to run CMake again, if you have changed the ``CMakeLists.txt``.

You can then execute the application in the usual way

* on Makefile-based builds (Linux/Mac/BSD):

    .. code-block:: console

        # ./my_project

* Windows

    .. code-block:: console

        # my_project

Using IDEs
^^^^^^^^^^

On Linux and BSD many IDEs directly support cmake, just open/import the ``CMakeLists.txt`` with e.g. `KDevelop <https://www.kdevelop.org>`_ or `QtCreator <http://www.qt.io/ide/>`_.

To use XCode on Mac with your CMake-based project, add ``-G Xcode`` to the cmake call above and then run ``open TODO``.

On Windows a Visual Studio generator is used by default and you will find a ``.vcxproj`` in the source directory that you can open with Visual Studio.

See :ref:`this page <infra-use-cmake-build-dirs>` for more details.


Details of the FindSeqAn Module
-------------------------------

As mentioned above, this line is the important line for including SeqAn:

.. code-block:: cmake

    find_package (SeqAn REQUIRED)

If SeqAn is only an optional dependency of your program, you can omit the ``REQUIRED`` keyword. In this case you should check the contents of the ``SEQAN_FOUND`` CMake-variable and depending on that configure your build, e.g. with custom Macros.

You can also check for the definition of SeqAn's version macros from within your code:

``SEQAN_VERSION_STRING``
  Concatenated version string, ``${SEQAN_VERSION_MAJOR}.${SEQAN_VERSION_MINOR}.${SEQAN_VERSION_PATCH}``

``SEQAN_VERSION_MAJOR``
  Major version.

``SEQAN_VERSION_MINOR``
  Minor version.

``SEQAN_VERSION_PATCH``
  Patch-level version.

Dependencies
^^^^^^^^^^^^

SeqAn itself has some optional dependencies.
Certain features in SeqAn will be enabled or disabled, depending on whether the dependencies could be found.

.. caution::

    Optional dependencies of SeqAn have to be searched **before** the SeqAn module is searched!

Currently, the following dependencies enable optional features:

``ZLIB``
  zlib compression library

``BZip2``
  libbz2 compression library

``OpenMP``
  OpenMP language extensions to C/C++

An example of where you only want ZLIB and OpenMP support, but not BZip2, would look like this:

.. code-block:: cmake

    find_package (ZLIB)
    find_package (OpenMP)
    find_package (SeqAn)

From within CMake you can check the variables ``ZLIB_FOUND`` or ``OpenMP_FOUND`` to see the results of these dependency searches, but you can also use the following macros from within your source code to escape certain optional code paths:

``SEQAN_HAS_ZLIB``
  ``TRUE`` if zlib was found.

``SEQAN_HAS_BZIP2``
  ``TRUE`` if libbz2 was found.

``SEQAN_HAS_OPENMP``
  ``TRUE`` if OpenMP was found.

CMake build variables
^^^^^^^^^^^^^^^^^^^^^

As can be seen from the example above, the following variables need to be passed to ``include_directories()``, ``target_link_directories()``, and ``add_definitions()`` in your ``CMakeLists.txt``:

``SEQAN_INCLUDE_DIRS``
  A list of include directories.

``SEQAN_LIBRARIES``
  A list of libraries to link against.

``SEQAN_DEFINITIONS``
  A list of definitions to be passed to the compiler.

Required additions to C++ compiler flags are in the following variable:

``SEQAN_CXX_FLAGS``
  C++ Compiler flags to add.
