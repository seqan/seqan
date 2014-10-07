.. _tutorial-getting-started:

Getting Started
---------------

This chapter gives you the necessary steps to get started with SeqAn:

-  Necessary Prerequisites
-  Installing SeqAn from Subversion
-  Creating a first build.
-  Creating your own first application.

Use the following links to select your target operating system and IDE/build system.
The bold items show the recommended build system for the given platforms.

Linux using
  * :ref:`Makefiles <tutorial-getting-started-linux-makefiles>`
  * :ref:`Eclipse <tutorial-getting-started-linux-eclipse>`

Mac Os X
  * :ref:`Makefiles <tutorial-getting-started-mac-makefiles>`
  * :ref:`Xcode <tutorial-getting-started-mac-xcode>`

Windows
  * :ref:`Visual Studio <tutorial-getting-started-windows-visual-studio>`

Click "more" for details on the supported development platforms.

.. container:: foldable

   .. note ::
      **In-Depth Information:** Supported OS, Build Systems, and Compilers

      The content of this box is meant as additional information.
      You do not need to understand it to use SeqAn or follow the tutorials.

      There are three degrees of freedom when selecting a SeqAn development platform.
      The degrees of freedom are:

      #. The **operating system**.
         We support Linux, Mac Os X and Windows.
      #. The **build system**.
         This is partially orthogonal to the operating system, although each build system is only available on some platforms (e.g. Visual Studio is only supported on Windows).
         We use CMake to generate the actual build files and the build system maps to "CMake generators".
         A CMake generator creates either build files for a build system (e.g. GNU Make) or a project file for an IDE (e.g. for Visual Studio 2008).
      #. The **compiler**.
         This is partially orthogonal to the operating system and build system, although only some combinations of each are possible.
         For example, Visual Studio projects of a particular version can only use the Visual Studio compiler of the same version.

    The SeqAn team offers support for the following operating systems, build systems, and compilers:

    * **Operating System:** Linux, Mac Os X, Windows.
    * **Build System:** Makefiles, Visual Studio projects, XCode projects, Eclipse CDT projects.
    * **Compilers:** GNU g++ from version 4.1, LLVM/Clang from version 3.0, Visual C++ from Version 8.

    We are told that SeqAn also works on FreeBSD.
    It should work with all `generators available in CMake <http://www.cmake.org/cmake/help/v2.8.8/cmake.html#section_Generators>`_ that work with the supported compilers (e.g. the CodeBlocks generator will probably work as long as you use it on a operating system with a supported compiler, although we cannot offer any support for CodeBlocks).

Relevant How-Tos
~~~~~~~~~~~~~~~~

Although slightly more advanced than "getting started", the following
How-Tos apply to setting up your build environment:

* :ref:`how-to-use-parallel-build-directories`
* :ref:`how-to-install-contribs-on-windows`
* :ref:`build-manual-integration-with-your-own-build-system`

.. toctree::
   :hidden:
   :maxdepth: 2

   GettingStarted/LinuxMakefiles
   GettingStarted/LinuxEclipse
   GettingStarted/MacMakefiles
   GettingStarted/MacXcode
   GettingStarted/WindowsVisualStudio
