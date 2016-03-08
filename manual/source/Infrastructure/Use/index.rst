.. _infra-use:

User Guide
==========

Before you can start using the SeqAn library you have to make sure that it is installed and setup correctly:

  * :ref:`Installing SeqAn <infra-use-install>`
  * :ref:`Installing Dependencies <infra-use-install-contribs-windows>`
Once this is done you have two choices for using SeqAn:

  * :ref:`Using SeqAn in CMake-based projects <infra-use-cmake>`
  * :ref:`Integration with your own Build System <infra-use-custom>`

We highly recommend using the CMake-Module as it correctly handles many settings you otherwise need to set manually. It also provides the best cross-platform compatibility.

If you choose to use CMake, the following document contains more information on build configurations and using multiple build directories in parallel:

  * :ref:`CMake Build Directories <infra-use-cmake-build-dirs>`


.. toctree::
    :glob:
    :titlesonly:
    :hidden:

    Install.rst
    FindSeqAnCMake.rst
    CustomBuildSystem.rst
    CMakeBuildDirs.rst
    InstallContribsWindows.rst

