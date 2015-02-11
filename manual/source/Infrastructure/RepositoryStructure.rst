.. sidebar:: ToC

   .. contents::


.. _infrastructure-repository-structure:

SeqAn Repository Structure
--------------------------

This article describes the SeqAn repository structure.
After reading it, you will have knowledge about the repository structure and the reasons for the design decisions.

Note that this article describes the structure of the Subversion repository, not the structure of the release version of SeqAn which you can download as a ZIP archive.

Overview
~~~~~~~~

The main repository structure is shown in the following picture.

::

    seqan
      |-- CMakeLists.txt      CMake script file.
      |
      |-- LICENSE             Top-Level Information Files
      |-- README.rst
      |
      |-- apps              Apps area
      |
      |-- demos             Demos area
      |
      |-- dox               Documentation System
      |
      |-- include/seqan     Core Library
      |
      |-- manual            Tutorials
      |
      |-- tests             Unit Tests for Library Modules
      |                
      `-- util				Miscellaneous and Utility Code	

* The repository root contains some **information files** such as the ``LICENSE`` and ``README.rst``.
  The file names should speak for themselves.
* The folder ``apps`` contains the applications that are shipped together with the library.
  Each application directory contains the source files for one or more binaries, documentation, example files, and app tests.
  More information is available in the section `Application Structure`_.
* The folder ``demos`` contains **demo programs**.
* The folder ``dox`` contains the documentation building system, which builds the documentation from the API documentation.
* The folders ``include/seqan`` contains main library code, with all modules.
  This is described in more detail in the section `Library Modules`_.
* The folders ``manual`` contains the tutorials and the building system to build the tutorials.
* The folder ``tests`` contains the unit tests for the library modules. 
  For each library module, there is a directory below ``tests`` with the same name that contains the tests for this module.
  Simpler modules have one tests executable, whereas there might be multiple tests executables for larger modules.
  For example, the module ``index`` has multiple test programs ``test_index_qgram``, ``test_index_shapes`` etc.
  Writing tests is explained in detail in the article :ref:`how-to-write-tests`.
* The folder ``util`` contains **miscellaneous files** and **utility code**.  

Application Structure
~~~~~~~~~~~~~~~~~~~~~

Each application directory contains one ``CMakeLists.txt`` file and the files for compiling one binary.
Usually, apps have tests, too.
In this case, there is a subdirectory ``tests``.
Writing application tests is covered in detail in the article :ref:`how-to-write-app-tests`.

The general structure of an app is as follows:

::

    seqan/apps/razers
      |-- CMakeLists.txt      CMake script file
      |
      |-- README              Documentation and License Files
      |-- LICENSE
      |
      |-- example             Small Example Files
      |     |-- genome.fa
      |     |-- reads.fa
      |     `-- ...
      |
      |-- razers.cpp          Source Files for Executables
      |-- razers.h
      |-- ...
      |
      `-- tests               App Tests Files

Library Modules
~~~~~~~~~~~~~~~

The library modules area looks as follows:

::

    include/
      |-- seqan/
      |     |-- basic/                       Library Module basic
      |     |     |-- aggregate_concept.h
      |     |     |-- debug_test_system.h
      |     |     `-- ...
      |     |-- basic.h
      |     |
      |     |-- sequence/                    Library Module sequence
      |     |-- sequence.h
      |     |
      |     `-- ...                          Other Library Modules

On the top level, there is the folder ``seqan`` that contains the
library modules. Inside the folder ``seqan``, there is one directory and
one header for each module.

The folder ``<module-name>`` contains the headers for the module module-name.
The header ``<module-name>.h`` includes the headers from the module module-name.
Including the header makes the code in the module available.

Documentation System
~~~~~~~~~~~~~~~~~~~~

The folder ``dox`` is used for building the documentation with the old and the new documentation system.
You can build them by going into the directory and then calling ``./dox_only.sh``.
This will build the documentation into the sub directory ``html``:

.. code-block:: console

   seqan # cd dox
   dox # ./dox_only.sh
   ...
