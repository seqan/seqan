.. sidebar:: ToC

   .. contents::


.. _tutorial-getting-started-linux-makefiles:

Getting Started With SeqAn On Linux Using Makefiles
---------------------------------------------------

This tutorial explains how to get started with SeqAn on Linux using Makefiles.

We assume that you are using the Debian or a Debian-like Linux distributions such as Ubuntu.
The only difference to other distributions is the name of the packages and the package management system in the Prerequisites section.
It should be very simple for you to tailor these instructions to your requirements.

Prerequisites
~~~~~~~~~~~~~

Use the following command line to install the required dependencies: the `Git`__ client, the `GNU C++ compiler <http://gcc.gnu.org/>`_, `CMake <http://cmake.org>`_ for the build system and the `Python <http://python.org>`_ script interpreter for running helper scripts.

.. __: http://git-scm.com/

.. code-block:: console

    ~ # sudo apt-get install git g++ cmake python

The following command line installs optional dependencies:
developer versions of `zlib <http://zlib.org>`_ and `libbzip2 <http://bzip.org>`_ (for compressed I/O support) and the `Boost <http://boost.org>`_ library (required by a few apps).

.. code-block:: console

    ~ # sudo apt-get install zlib1g-dev libbz2-dev libboost-dev

Install
~~~~~~~

.. important::
	
	In the following we describe the easiest way to get up and running with SeqAn.
	This is especially recommended for novel users working through the tutorials in the beginning.
	If you are planning to contribute to SeqAn at any point, you need to read the :ref:`infrastructure-seqan-git-workflow` instructions first. 
	This manual will guide you through the SeqAn workflow required to submit bug-fixes and new features.

Go to the directory you want to keep your SeqAn install in (e.g. ``Development`` in your home folder).

.. code-block:: console

    ~ # cd $HOME/Development

Then, use git to retrieve the current SeqAn source-base:

.. code-block:: console

    # Development # git clone https://github.com/seqan/seqan.git seqan-src

You can now find the whole tree with the SeqAn library and applications in ``$HOME/Development/seqan-src``.

.. tip::

    By default git creates a local branch pointing to the stable master branch.
    This branch is only updated when hot fixes are applied or a new release is published.
    
    If you want to have access to regular updates and new features you can switch to the ``develop`` branch of SeqAn:
    
    .. code-block:: console

		# Development # cd seqan-src
		# Development/seqan-src # git checkout -b develop origin/develop
	
    For more help on git, please read the documentation ``git help`` and consult the homepage `Git`__.

.. __: http://git-scm.com/

.. warning::

    Note that the state of develop is not guaranteed to be stable at any time.


A First Build
~~~~~~~~~~~~~

Next, we will use CMake to create Makefiles for building the applications, demo programs (short: demos), and tests.
For this, we create a separate folder ``seqan-build`` on the same level as the folder ``seqan-src``.

.. code-block:: console

    # Development # mkdir seqan-build

When using Makefiles, we have to create separate Makefiles for debug builds (including debug symbols with no optimization) and release builds (debug symbols are stripped, optimization is high).
Thus, we create a subdirectory for each build type.

We start with debug builds since this is best for learning:
debug symbols are enabled and assertions are active.

.. warning ::

    Compiling **debug mode yields very slow binaries** since optimizations are disabled.
    Compile your programs in release mode if you want to run them on large data sets.

    The reason for disabling optimizations in debug mode is that the compiler performs less inlining and does not optimize variables away.
    This way, debugging your programs in a debugger becomes much easier.

.. code-block:: console

    # Development # mkdir seqan-build/debug
    # Development # cd seqan-build/debug

The resulting directory structure will look as follows.

::

       ~/Development
         ├─ seqan-src                 source directory
         └─ seqan-build
            └─ debug                  build directory with debug symbols

Within the **build directory** ``debug``, we call CMake to generate Makefiles in *Debug* mode.

.. code-block:: console

    # debug # cmake ../../seqan-src -DCMAKE_BUILD_TYPE=Debug

We can then build one application, for example RazerS 2:

.. code-block:: console

    # debug # make razers2

Optionally, we could also use "``make``\ " instead of "``make razers2``\ ".
However, this builds all demos, tests and applications, which **can take a long time and is not really necessary**.

Hello World!
~~~~~~~~~~~~

Now it is time to write your first little application in SeqAn.
Go to the demos folder in the ``seqan-src`` directory and create a new folder with the same name as your username.
In this tutorial we use ``seqan_dev`` as the username.
Create a new cpp file called ``hello_seqan.cpp``

.. code-block:: console
	
    # debug # cd ../../seqan-src/demos
    # demos # mkdir seqan_dev; cd seqan_dev
    # seqan_dev # echo "" > hello_seqan.cpp

Now, we go back into the build directory and call CMake again to make it detect the new source file.

.. code-block:: console

    # seqan_dev # cd ../../../seqan-build/debug
    # debug # cmake .

.. tip::

    When and where do you have to call CMake?

    CMake is a cross-platform tool for creating and updating build files (IDE projects or Makefiles).
    When you first create the build files, you can configure things such as the build mode or the type of the project files.

    Whenever you add a new application, a demo or a test or whenever you make changes to ``CMakeLists.txt`` you need to call CMake again.
    Since CMake remembers the settings you chose the first time you called CMake in a file named ``CMakeCache.txt``, all you have to do is to switch to your ``debug`` or ``release`` build directory and call "``cmake .``" in there.

    .. code-block:: console

       ~ # cd $HOME/Development/seqan-build/debug
       # debug # cmake .

    Do not try to call "``cmake .``" from within the ``seqan-src`` directory **but only from your build directory**.

.. raw:: html

   </pre>

Open the file ``demos/seqan_dev/hello_seqan.cpp`` (in your ``seqan-src`` directory) with a text editor and replace its contents with the following:

.. code-block:: cpp

    #include <iostream>
    #include <seqan/sequence.h>  // CharString, ...
    #include <seqan/stream.h>    // to stream a CharString into cout

    int main(int, char const **)
    {
        std::cout << "Hello World!" << std::endl;
        seqan::CharString mySeqAnString = "Hello SeqAn!";
        std::cout << mySeqAnString << std::endl;
        return 1;
    }

Afterwards, you can simply compile and run your application:

.. code-block:: console

    # debug # make demo_seqan_dev_hello_seqan
    # debug # ./bin/demo_seqan_dev_hello_seqan

On completion, you should see the following output:

.. code-block:: console

    Hello World!
    Hello SeqAn!

Congratulations, you have successfully created your first application within the SeqAn build system with Makefiles!

Further Steps
~~~~~~~~~~~~~

As a next step, we suggest the following:

* :ref:`Continue with the Tutorials <tutorial>`
* For the tutorial, using the SeqAn build system is great!
  If you later want to use SeqAn as a library, have a look at :ref:`build-manual-integration-with-your-own-build-system`.
* If you plan to contribute to SeqAn, please read the following document: :ref:`infrastructure-seqan-git-workflow`.
