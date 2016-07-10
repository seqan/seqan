.. sidebar:: ToC

    .. contents::

.. _infra-use-cmake-build-dirs:

CMake Build directories in detail
=================================

Motivation
----------

Why would you need more than one build directory or more than one IDE project file?
This is very useful

* if you want to use the same set of source files from multiple version of the same IDE (e.g. two Visual Studio versions),
* if you want to have both debug builds (for debugging) and release builds (for performance tests) in parallel,
* if you have your source files stored on a shared network location and want to have build files on two computer and/or operating systems, or
* if you want to build the sources with two different compilers or compiler versions at the same time (e.g. to see whether you can figure out compiler errors better from the messages by another compiler).

The overall idea is very simple: you create one build directory for each variant and call CMake in each of it using different settings.

.. tip::

    A nice side-effect of separating source and build directories is also that you can just delete you build directory and recreate it if you feel that something went wrong configuring your build.

CMake Parameters
----------------

A central question with CMake is the choice of the so called generator. Enter ``cmake -G`` to get a list of the supported ones. The most common generators are the **Unix Makefiles** which are default on Linux/Mac/BSD. But there are also specific generators for IDEs, such as Visual Studio, XCode or CodeBlocks.

For most of the IDEs further choices like "Release or Debug" are available from the graphical user interface of the IDE, whereas, for the Unix Makefile generator, we can specify the *build types* using a command line option.
Also, the compiler program (and version) can be switched using a command line option.

Examples
--------

We assume that your project source is at ``~/devel/my_project``.

Unix Makefiles
^^^^^^^^^^^^^^

Different compilers:

.. code-block:: console

    # mkdir -p ~/devel/my_project-build/release_gcc5
    # cd ~/devel/my_project-build/release_gcc5
    # cmake ../../my_project -DCMAKE_CXX_COMPILER=g++-5
    [...]
    # mkdir -p ~/devel/my_project-build/release_clang37
    # cd ~/devel/my_project-build/release_clang37
    # cmake ../../my_project -DCMAKE_CXX_COMPILER=clang++-3.7
    [...]

Please note that the above only works if your compiler is in your PATH. You can instead also specify a full path like ``-DCMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-4.9``.

Debug and release builds:

.. code-block:: console

    # mkdir -p ~/devel/my_project-build/release
    # cd ~/devel/my_project-build/release
    # cmake ../../my_project
    [...]
    # mkdir -p ~/devel/my_project-build/debug
    # cd ~/devel/my_project-build/debug
    # cmake ../../my_project -DCMAKE_BUILD_TYPE=Debug
    [...]

Of course the above can also be combined to have ``debug_clang37`` et cetera.

Visual Studio
^^^^^^^^^^^^^

Different versions (please note that versions older than 2015 are not supported any longer):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: console

    # mkdir -p ~/devel/my_project-build/vs2015
    # cd ~/devel/my_project-build/vs2015
    # cmake ../../my_project -G "Visual Studio 14 2015"
    [...]
    # mkdir -p ~/devel/my_project-build/vs2013
    # cd ~/devel/my_project-build/vs2013
    # cmake ../../my_project -G "Visual Studio 12 2013"
    [...]


32Bit and 64Bit:
~~~~~~~~~~~~~~~~

.. code-block:: console

    # mkdir -p ~/devel/my_project-build/vs2015_32
    # cd ~/devel/my_project-build/vs2015_32
    # cmake ../../my_project -G "Visual Studio 14 2015"
    [...]
    # mkdir -p ~/devel/my_project-build/vs2015_64
    # cd ~/devel/my_project-build/vs2015_64
    # cmake ../../my_project -G "Visual Studio 14 2015 Win64"
    [...]

.. caution::

    **64Bit builds on Windows**

    You almost always want 64Bit builds when using SeqAn, so don't forget to specify a generator that ends in "Win64". It is not the default, even on 64Bit Windows installations.

Different Compilers:
~~~~~~~~~~~~~~~~~~~~

`Intel Compiler 2016
<https://software.intel.com/en-us/articles/intel-parallel-studio-xe-2016-release-notes>`_:

.. code-block:: console

    # mkdir -p ~/devel/my_project-build/intel_32
    # cd ~/devel/my_project-build/intel_32
    # cmake ../../my_project -G "Visual Studio 14 2015" -T "Intel C++ Compiler 16.0"
    [...]

    # mkdir -p ~/devel/my_project-build/intel_64
    # cd ~/devel/my_project-build/intel_64
    # cmake ../../my_project -G "Visual Studio 14 2015 Win64" -T "Intel C++ Compiler 16.0"
    [...]

`Clang/C2 3.7 or 3.8
<https://blogs.msdn.microsoft.com/vcblog/2015/12/04/clang-with-microsoft-codegen-in-vs-2015-update-1/>`_
(requires CMake ≥ 3.6):

.. code-block:: console

    # mkdir -p ~/devel/my_project-build/clang_c2_32
    # cd ~/devel/my_project-build/clang_c2_32
    # cmake ../../my_project -G "Visual Studio 14 2015" -T "v140_clang_3_7"
    [...]

    # mkdir -p ~/devel/my_project-build/clang_c2_64
    # cd ~/devel/my_project-build/clang_c2_64
    # cmake ../../my_project -G "Visual Studio 14 2015 Win64" -T "v140_clang_3_7"
    [...]

.. note::

    If Clang/C2 3.8 is installed, the tool-chain name in Visual Studio 14 is
    still "v140_clang_3_7" even though the name says otherwise.

.. caution::

    Clang/C2 is currently experimental and shouldn't be used in production.

XCode
^^^^^

.. code-block:: console

    # mkdir -p ~/devel/my_project-build/xcode
    # cd ~/devel/my_project-build/xcode
    # cmake ../../my_project -G "Xcode"
    [...]
