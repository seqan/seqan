.. sidebar:: ToC

   .. contents::


.. _tutorial-getting-started-mac-makefiles:

Getting Started With SeqAn On Mac Os X Using Makefiles
------------------------------------------------------

This tutorials explains how to get started with SeqAn on Mac Os X using Makefiles.

We assume that you want to use `MacPorts <http://www.macports.org/>`_ for installing some dependencies (MacPorts is a package management system that easily allows you to install Unix software on Os X).
Of course, if you want to use a different way for installing the dependencies (e.g. using Homebrew) then you are free to do so.

Prerequisites
~~~~~~~~~~~~~

First, you have to install the Apple `Xcode SDK <https://developer.apple.com/downloads/index.action>`_ (Apple ID needed).

.. warning::

    Please choose Xcode SDK version 4.2 or lower, because the current version has some compatibility problems with the SeqAn build system.

After installing the Xcode SDK, please install MacPorts following `these instructions <http://www.macports.org/install.php>`_.
To check that the MacPorts install was successful, enter the following on your shell.
If the ``port`` program is found then you can go on.

.. code-block:: console

    # port info

Next, install `CMake <http://cmake.org>`_ and `Subversion <http://subversion.apache.org/>`_ using the ``port`` command.

.. code-block:: console

    # port install cmake
    # port install subversion

Install
~~~~~~~

Now, go to the directory you want to keep your SeqAn install in (e.g.  ``Development`` in your home folder).

.. code-block:: console

    ~ # cd $HOME/Development

Then, use Subversion to retrieve the current SeqAn trunk:

.. code-block:: console

    Development # svn co https://github.com/seqan/seqan/branches/master seqan-trunk

You can now find the whole tree with the SeqAn library and applications in ``seqan-trunk``.

A First Build
~~~~~~~~~~~~~

Next, we will use CMake to create Makefiles for building the applications, demo programs (short: demos), and tests.
For this, we create a separate folder ``seqan-trunk-build`` on the same level as the folder ``seqan-trunk``.

.. code-block:: console

    Development # mkdir seqan-trunk-build

When using Makefiles, we have to create separate Makefiles for debug builds (including debug symbols with no optimization) and release builds (debug symbols are stripped, optimization is high).
Thus, we create a subdirectory for each build type.
We start with debug builds since this is best for learning: Debug symbols are enabled and assertions are active

.. warning::

    Compiling ''debug mode yields very slow binaries''' since optimizations are disabled.
    Compile your programs in release mode if you want to run them on large data sets.

    The reason for disabling optimizations in debug mode is that the compiler performs less inlining and does not optimize variables away.
    This way, debugging your programs in a debugger becomes much easier.

.. code-block:: console

    Development # mkdir seqan-trunk-build/debug
    Development # cd seqan-trunk-build/debug

The resulting directory structure will look as follows.

::

       ~/Development
         +-- seqan-trunk        source directory
         `-- seqan-trunk-build
             `-- debug          build directory with debug symbols

Within the **build directory** ``debug``, we use CMake to generate Makefiles in *Debug* mode.

.. code-block:: console

    debug # cmake ../../seqan-trunk -DCMAKE_BUILD_TYPE=Debug

We can then build one application, for example RazerS 2:

.. code-block:: console

    debug # make razers2

Optionally, we could also use "``make``" instead of "``make razers2``". However, this **can take a long time and is not really necessary**.

Hello World!
~~~~~~~~~~~~

Now, let us create a **sandbox** for you.
This sandbox will be your local workspace and you might want to have it versionized on your own Subversion repository at a later point.
All of your development will happen in your sandbox.

We go back to the source directory and then use the SeqAn code generator to create a new sandbox.

.. code-block:: console

    debug # cd ../../seqan-trunk
    seqan-trunk # ./util/bin/skel.py repository sandbox/my_sandbox

Now that you have your own working space, we create a new application ``first_app``.

.. code-block:: console

    seqan-trunk # ./util/bin/skel.py app first_app sandbox/my_sandbox

Details about the code generator are explained in :ref:`how-to-use-the-code-generator`.

Now, we go back into the build directory and call CMake again to make it detect the added app.

.. code-block:: console

    seqan-trunk # cd ../seqan-trunk-build/debug
    debug # cmake .

.. tip::

    When and where do you have to call CMake?

    CMake is a cross-platform tool for creating and updating build files (IDE projects or Makefiles).
    When you first create the build files, you can configure things such as the build mode or the type of the project files.

    Whenever you add a new application, a demo or a test or whenever you make changes to ``CMakeLists.txt`` you need to call CMake again.
    Since CMake remembers the settings you chose the first time you called CMake in a file named ``CMakeCache.txt``, all you have to do is to switch to your ``debug`` or ``release`` build directory and call "``cmake .``" in there.

    .. code-block: console

       ~ # cd $HOME/Development/seqan-trunk-build/debug
       debug # cmake .

    Do not try to call "``cmake .``" from within the ``seqan-trunk`` directory **but only from your build directory**.

The step above creates the starting point for a real-world application, including an argument parser and several other things that are a bit too complicated to fit into the Getting Started tutorial.
Therefore, we will replace the program of the app *first_app* with a very simple example program.

Open the file ``sandbox/my_sandbox/apps/first_app/first_app.cpp`` (in your ``seqan-trunk`` directory) with a text editor and replace its contents with the following:

.. code-block:: cpp

    #include <iostream>
    #include <seqan/sequence.h>  // CharString, ...
    #include <seqan/file.h>      // to stream a CharString into cout

    int main(int, char const **)
    {
        std::cout << "Hello World!" << std::endl;
        seqan::CharString mySeqAnString = "Hello SeqAn!";
        std::cout << mySeqAnString << std::endl;
        return 1;
    }

Afterwards, you can simply compile and run your application:

.. code-block:: console

    debug # make first_app
    debug # ./bin/first_app

On completion, you should see the following output:

.. code-block:: console

    Hello World!
    Hello SeqAn!

Congratulations, you have successfully created your first application within the SeqAn build system with Makefiles!

Further Steps
~~~~~~~~~~~~~

As a next step, we suggest the following:

* :ref:`Continue with the Tutorials <tutorial>`
* Look around in the files in ``sandbox/my_sandbox/apps/first_app`` or the demos in ``core/demos`` and ``extras/demos``.
* For the tutorial, using the SeqAn build system is great!
  If you later want to use SeqAn as a library, have a look at :ref:`build-manual-integration-with-your-own-build-system`.

