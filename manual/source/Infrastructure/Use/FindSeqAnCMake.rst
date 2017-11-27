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

Cmake supports two different modes to load settings from an external project: The **module** and the **config** mode. 
Please read the `cmake documentation <https://cmake.org/cmake/help/v3.0/command/find_package.html>`_ to learn more about this feature.

Install SeqAn from package maintainer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The recommended way for SeqAn 2.3 or newer is to use the config mode. 
If you installed/updated SeqAn from one of the downstream package maintainer listed in :ref:`Getting Started with SeqAn <infra-use-install>`, then a file called ``seqan-config.cmake`` was installed in a system path that is automatically searched by the cmake system (see the cmake documentation for `find_package <https://cmake.org/cmake/help/v3.0/command/find_package.html>`_).
If everything was done with default settings, than you can simply build your project like:

.. code-block:: console
   
   # cmake ../../my_project

Install SeqAn into user defined prefix or clone from GitHub
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In case you obtained SeqAn from a git clone, or installed SeqAn into a user defined location, you need to specify the install location by setting the ``CMAKE_PREFIX_PATH`` in your cmake call.
In addition you also have to specify the ``SEQAN_INCLUDE_PATH`` variable to find the SeqAn headers. 
Assume you have cloned SeqAn into ``~/devel/seqan``, then your setup could look as the following:

.. code-block:: console
   
   # cmake ../../my_project \
      -DCMAKE_PREFIX_PATH="$HOME/devel/seqan/util/cmake" \
      -DSEQAN_INCLUDE_PATH="$HOME/devel/seqan/include"

Backwards compatibility
~~~~~~~~~~~~~~~~~~~~~~~

Before SeqAn 2.3 we used the module mode to setup SeqAn as an external project.
To allow backwards compatibility we added a redirect from the ``FindSeqAn.cmake`` to ``seqan-config.cmake`` in our sources.
In this case configuing your project with the old approach using the ``CMAKE_MODULE_PATH`` variable, will still work:

.. code-block:: console
   
   # cmake ../../my_project \
      -DCMAKE_MODULE_PATH="$HOME/devel/seqan/util/cmake" \
      -DSEQAN_INCLUDE_PATH="$HOME/devel/seqan/include"

.. tip::

    Depending on your setup you might need to manually choose a more modern compiler! Please read :ref:`this page <infra-use-cmake-build-dirs>` for more information on configuring CMake builds. Don't forget to clean your CMake build directory after changing the compiler!

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

Checking for newer Versions of SeqAn (optional)
-----------------------------------------------

The argument parser has a new feature to check for updates for the SeqAn library or for an application.
This can be a very helpful reminder to stay up to date since SeqAn evolves rapidly to resolve issues or to supply new functionality.
If none of the following options are selected the version update feature is activated by default.

  =================================  ==========================================
            Cmake Option                                Description
  =================================  ==========================================
  ``-DSEQAN_VERSION_CHECK_OPT_IN``   Turn update feature on but make it opt-in.

  ``-DSEQAN_DISABLE_VERSION_CHECK``  Turn update feature off.
  =================================  ==========================================

.. note::

    This does only affect applications or scipts that use the SeqAn :ref:`Argument Parser <tutorial-getting-started-parsing-command-line-arguments>`!

Details of the SeqAn Module
---------------------------

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

``_OPENMP``
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

  .. caution::

    Please note that these variables include whatever has been added by the dependencies mentioned above so **do not add** e.g. ``${OpenMP_CXX_FLAGS}`` yourself!

Intel Compiler specifics
^^^^^^^^^^^^^^^^^^^^^^^^

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
``icpc -gxx-name=g++-5.4.0 -gcc-name=gcc-5.4.0 …``. Use
``cmake -DCMAKE_CXX_FLAGS="-gxx-name=g++-5.4.0 -gcc-name=gcc-5.4.0" …``
to propagate those options through cmake.

You may have to add the path of the library to ```$LD_LIBRARY_PATH`` for the
linker.

Static builds
^^^^^^^^^^^^^

If you want to build your app statically, please do not use gcc-4.9 or make sure you add the ``-static`` flag **before** calling ``find_package (SeqAn)``. Otherwise a broken binary will be built that crashes immediately.
