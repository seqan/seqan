.. sidebar:: ToC

   .. contents::


.. _how-to-use-parallel-build-directories:

Using Parallel Build Directories
--------------------------------

Motivation
~~~~~~~~~~

Why would you need more than one build directory or more than one IDE project file?
This is very useful

* if you want use the same set of source files from multiple version of the same IDE (e.g. two Visual Studio versions),
* if you want to have both debug builds (for debugging) and release builds (for performance tests) in paralell,
* if you have your source files stored on a shared network location and want to have build files on two computer and/or operating systems, or
* if you want to build the sources with two different compilers or compiler versions at the same time (e.g. to see whether you can figure out compiler errors better from the messages by another compiler).

This How-To also serves as a collection of CMake command lines for copy-and-paste.

The Overall Idea
~~~~~~~~~~~~~~~~

The overall idea is very simple: you create one build directory for each variant and call CMake in each of it using different settings.

If you want to have different IDE project files then you use different *CMake generators*.
In most IDEs, there is an option to select debug or release builds.
For the CMake Makefile generator, however, we can select the *build types* using a command line option.
Also, the compiler program (and version) can be switched using a command line option.

Generating Parallel IDE Project Files (Visual Studio, Xcode etc.)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You will only be able to generate files for Xcode when on Mac Os X, for Visual Studio when on Windows and so on.

The following section assumes that you have a subdirectory of your SeqAn checkout called ``build``.
We will create subdirectories for each IDE project we create.

**Creating Directories**

For example, if we have installed Visual studio 8, 9, and 10 and want to create project files for each, we might use the following commands to create the directories:

.. code-block:: console

    > cd seqan-trunk\build
    > mkdir vs8
    > mkdir vs9
    > mkdir vs10

For XCode on Mac, we could do the following:

.. code-block:: console

    $ cd seqan-trunk/build
    $ mkdir xcode

Note that you can choose any directory name.
You have to take care that no such directory exists before.
Previously generated project files can break the generation process!

**Generating Project Files**

We can now use CMake to generate the project fiels specifying a generator with the command line parameter ``-G``.

Let us generate the Visual Studio projects in the directories we mentioned above:

.. code-block:: console

    > cd vs8
    > cmake -G "Visual Studio 8 2005" ..\..\..
    > cd ..\vs9
    > cmake -G "Visual Studio 9 2008" ..\..\..
    > cd ..\vs10
    > cmake -G "Visual Studio 10" ..\..\..

Click *more...* to see the commands for 64 bit builds.

.. container:: foldable

   .. code-block:: console

      > cd vs8
      > cmake -G "Visual Studio 8 2005 Win64" ..\..\..
      > cd ..\vs9
      > cmake -G "Visual Studio 9 2008 Win64" ..\..\..
      > cd ..\vs10
      > cmake -G "Visual Studio 10 Win64" ..\..\..

On Mac Os X, we can generate XCode build files as follows:

.. code-block:: console

    # cd xcode
    # cmake -G "Xcode" ../../..

How To Make Release/Debug builds
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using the Makefile generator, it is useful to have multiple build
types. CMake supports the following build types:

**Different Build Types**

Debug
  No optimization, with debug symbols.

Release
 Optimization, without debug symbols.

RelWithDebInfo
  Optimization, with debug symbols.
  Useful for profiling.

MinSizeRel
  Size-optimized release binary without debug symbols.

You can select the build type with a command line parameter to ``cmake``, e.g. ``-DCMAKE_BUILD_TYPE=Debug`` or ``-DCMAKE_BUILD_TYPE=Release``.

**Picking A Compiler**

You can pick a C++ compiler using the command line parameter to ``cmake``, e.g. ``-DCMAKE_CXX_COMPILER=g++-4.1`` or ``-DCMAKE_CXX_COMPILER=clang++``.

**Creating Directories**

Let's create a build directory with the system's default compiler both in debug and release mode.
Also, we create one directory for the Clang compiler in debug mode.

.. code-block:: console

    # cd seqan-trunk/build
    # mkdir Debug
    # mkdir Release
    # mkdir Debug-clang

Note that you should use fresh directories.
Previously generated Makefiles can break the generation process!

**Generating Project Files**

.. code-block:: console

    # cd Debug
    # cmake ../..
    # cd ../Release
    # cmake -DCMAKE_BUILD_TYPE=Release ../..
    # cd ../Debug-clang
    # cmake -DCMAKE_CXX_COMPILER=clang++

Note that when using clang, you actually have to use ``clang++`` and not ``clang`` (although ``clang++`` usually only is a symlink to ``clang``).
If you use ``clang`` then all C++ features will be disabled and you will get configuration errors.
