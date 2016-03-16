.. _infra-use:

User Guide
==========

Before you can start using the SeqAn library you have to make sure that it is installed and set up correctly. Optionally also install the dependencies:

  * :ref:`Installing SeqAn <infra-use-install>`
  * :ref:`Installing Dependencies <infra-use-install-dependencies>`

Once this is done you have two choices for using SeqAn:

  * :ref:`Using SeqAn in CMake-based projects <infra-use-cmake>`
  * :ref:`Integration with your own Build System <infra-use-custom>`

We highly recommend using the CMake-Module as it correctly handles many settings you otherwise need to set manually. It also provides the best cross-platform compatibility.

If you choose to use CMake, the following document contains more information on build configurations and using multiple build directories in parallel:

  * :ref:`CMake Build Directories <infra-use-cmake-build-dirs>`


Some notes on using this manual
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You will often see something like this in the manual:

.. code-block:: console

   # mkdir -p /tmp/mytempdir

This is called a terminal and you can find it on all major operating systems, sometimes also called "Konsole" or "Shell" (although that is something different strictly speaking).

You are expected to enter whatever comes right of the ``#`` (sometimes also the ``$``) and then press ``RETURN``. Knowing your way around the Terminal will make things easier, but it should be possible to just copy'n'paste.

.. important::

    On Windows you are expected to be using the **PowerShell**, and not the legacy command prompt. Our tutorials will only work with the PowerShell!

.. toctree::
    :glob:
    :titlesonly:
    :hidden:

    Install.rst
    InstallDependencies.rst
    FindSeqAnCMake.rst
    CustomBuildSystem.rst
    CMakeBuildDirs.rst
