.. sidebar:: ToC

    .. contents::

.. _internal-infrastructure-repository-structure:

SeqAn Repository Structure
==========================

This article describes the SeqAn repository structure.
After reading it, you will have knowledge about the repository structure and the reasons for the design decisions.

Note that this article describes the structure of the Git repository and the "full sources", not the "library package" available from our downloads.

Overview
--------

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
      |-- manual            Manuals
      |
      |-- tests             Unit Tests for Library Modules
      |
      `-- util              Miscellaneous and Utility Code

* The repository root contains some **information files** such as the ``LICENSE`` and ``README.rst``.
  The file names should speak for themselves.
* The folder ``apps`` contains the applications that are shipped together with the library.
  Each application directory contains the source files for one or more binaries, documentation, example files, and app tests.
  More information is available in the section `Application Structure`_.
* The folder ``demos`` contains **demo programs**.
* The folder ``dox`` contains the documentation building system, which builds the API documentation from the source code, see also section `API Documentation`_.
* The folders ``include/seqan`` contains main library code, with all modules.
  This is described in more detail in the section `Library Modules`_.
* The folders ``manual`` contains the tutorials and manuals you are reading right now, see section `Manual`_.
* The folder ``tests`` contains the unit tests for the library modules. 
  For each library module, there is a directory below ``tests`` with the same name that contains the tests for this module.
  Simpler modules have one tests executable, whereas there might be multiple tests executables for larger modules.
  For example, the module ``index`` has multiple test programs ``test_index_qgram``, ``test_index_shapes`` etc.
  Writing tests is explained in detail in the article :ref:`how-to-recipes-write-tests`.
* The folder ``util`` contains **miscellaneous files** and **utility code**.  

Application Structure
---------------------

Each application directory contains one ``CMakeLists.txt`` file and the files for compiling one binary.
Usually, apps have tests, too.
In this case, there is a subdirectory ``tests``.
Writing application tests is covered in detail in the article :ref:`how-to-recipes-write-app-tests`.

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
---------------

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


API Documentation
-----------------

The SeqAn API documentation is created using a customly-written system called *dox*.
You can find out more about the syntax in :ref:`internal-style-guide-dox-api-docs`.

You can build the documentation in the `dox` folder:

.. code-block:: console

   dox # ./dox_only.sh

This will build the documentation into the sub directory ``html``.


Manual
------

The SeqAn manual is created using the `Sphinx <http://sphinx-doc.org/>`_ documentation system.

Follow these instructions to setup a local sphinx environment and build the manual:

.. code-block:: console

    $ virtualenv ~/seqan-manual-env
    $ source ~/seqan-manual-env/bin/activate
    (seqan-manual-env) $ cd ~/seqan/manual
    (seqan-manual-env) $ pip install -r requirements.txt
    (seqan-manual-env) $ make html

Note that you have to first build the dox documentation since plugins for generating the ``:dox:`` links rely on the generated search index for checks.
In order to get correct dox-links within the generated manuals, you have to specify the correct branch version.
If you are working on the develop branch there is nothing to do, since ``'develop'`` is set by default.
But if you are working on another branch, for example ``master``, you can set the correct branch by calling

.. code-block:: console

    (seqan-manual-env) $ export READTHEDOCS_VERSION='master'

before you call ``make html`` as described in the previous step.
This will generate the correct links to the master's version of the dox, i.e., ``http://docs.seqan.de/seqan/master/``
