.. sidebar:: ToC

   .. contents::


.. _infrastructure-build-system:

The CMake-Based Build System
----------------------------

We use `CMake <http://www.cmake.org>`_ for building the SeqAn demos, applications and tests.
This wiki page explains our usage of CMake, the variables we define and how to extend the build scripts for new demos, apps etc.
`CMake's documentation (v2.8) <http://www.cmake.org/cmake/help/cmake-2-8-docs.html>`_ supplements this document.
The `documentation of CTest (v2.8) <http://www.cmake.org/cmake/help/ctest-2-8-docs.html>`_ could also be of interest.

.. todo:: Document usage of "dummy targets" for IDE project file generation.

Directory Layout
~~~~~~~~~~~~~~~~

The CMake files live in *projects/library/cmake*:

.. code-block:: console

    $ cd projects/library/cmake
    $ tree
    .
    |-- CMakeLists.txt
    |-- apps
    |   `-- CMakeLists.txt
    |-- demos
    |   `-- CMakeLists.txt
    |-- seqan
    |   `-- CMakeLists.txt
    `-- tests
        |-- CMakeLists.txt
        `-- CTestConfig.cmake

Target Structure
~~~~~~~~~~~~~~~~

There is a target for each program to be built.

Additionally, there is a target called *Seqan* that represents the library.
When using the GCC, we need to build the generated forward headers.
In this case, *Seqan* also depends on the generated forward headers and the *CMakeList.txt* files define how to generate these generated forwards.

External Dependencies
~~~~~~~~~~~~~~~~~~~~~

SeqAn is a C++ header library and thus does not need to build itself.
However, some applications have dependencies on external libraries, such as `Boost <http://www.boost.org>`_ or `Threading Building Blocks <http://www.threadingbuildingblocks.org/>`_.

The policy is to install these external dependencies on your system and let CMake find them.

The policy for missing dependencies is not to build the programs that depend on them and print an error message.


Adding New Programs
~~~~~~~~~~~~~~~~~~~

The process of adding a new demo, test or app is really simple: create a new directory *my_app* under *projects/library/app*, ''my_test" under *projects/test* or *my_demo* under *projects/demos*.
Within this directory create a new file *my_app.cpp*, *my_test.cpp* or *my_demo.cpp* and write your program.

Go to *projects/library/cmake* in your shell and execute ``cmake .`` again.
The new target will appear in your IDE.
If you use Makefiles then you can now type ``make my_app``, ``make my_test`` or ``make my_demo``.

Multiple Build Types
~~~~~~~~~~~~~~~~~~~~

You can call CMake in different directories to be able to build Debug, Release etc. binaries without having to re-cmake.
The process is described `here in the CMake wiki <http://www.vtk.org/Wiki/CMake_FAQ#How_can_I_build_multiple_modes_without_switching_.3F>`_.
