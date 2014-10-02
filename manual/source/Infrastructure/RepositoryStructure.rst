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
      |-- GETTING_STARTED
      |-- LICENSE             Top-Level Information Files
      |-- README
      |
      |-- core                Core Area
      |
      |-- extras              Extras Area
      |
      |-- sandbox             Sandboxes Area
      |
      |-- docs                Documentation System
      |-- docs2
      |
      |-- misc                Miscellaneous and Utility Code
      `-- util

* The repository root contains some **information files** such as the ``LICENSE``, ``README``, and ``GETTING_STARTED`` files.
  The file names should speak for themselves.
* The folder ``core`` contains the **core area** with apps, tests, and library modules that are (1) stable and (2) of general interest.
  Furthermore, it contains demos for the library modules in the SeqAn core.
* The folder ``extras`` contains the **extras area** with apps, tests, and library modules that are either (1) not stable enough yet or (2) of special interest only, as well as demos for the SeqAn extras library modules.
* The folder ``sandbox`` contains the **sandbox area**.
  Users can create their own **user areas** inside this folder as described in the section `Sandboxes`_.
* The folders ``docs`` and ``docs2`` contain the scripts for the **documentation system**.
  At the moment, there are two concurrent documentation systems.
  In the midterm future, we aim to replace this by one new documentation system.
* The folders ``misc`` and ``util`` contain **miscellaneous files** and **utility code**.
  For example the :ref:`Code Generator <how-to-use-the-code-generator>` Python scripts are located here as well as SeqAn logo image files and CMake modules.

Core Area
~~~~~~~~~

The core area is structured as follows. Note that we generally refer to such a structure as a **repository** below.

::

    seqan/core
      |-- CMakeLists.txt      CMake script file
      |
      |-- apps                Applications
      |
      |-- demos               Demo Programs
      |
      |-- include             SeqAn Library Modules
      |
      `-- tests               Tests for Library Modules

* The ``apps`` directory contains **applications**.
  Each application directory contains the source files for one or more binaries, documentation, example files, and app tests.
  More information is available in the section `Application Structure`_.
* The ``demos`` directory contains **demo programs**.
  The ``CMakeLists.txt`` file in this directory is written such that each file ending in ``.cpp`` is compiled into an executable with default SeqAn flag options.
* The ``include`` directory contains **library modules**.
  This is described in more detail in the section `Library Modules`_.
* The ``tests`` directory contains **tests for library modules**.
  For each library module, there is a directory below ``tests`` with the same name that contains the tests for this module.
  Simpler modules have one tests executable, whereas there might be multiple tests executables for larger modules.
  For example, the module ``index`` has multiple test programs ``test_index_qgram``, ``test_index_shapes`` etc.
  Writing tests is explained in detail in the article :ref:`how-to-write-tests`.

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

    seqan/core/include
      |-- seqan
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

Extras Area
~~~~~~~~~~~

The extras area has the same "repository" structure as the core area.
The main difference is that it contains code that is not stable enough or not of general interest.

Sandboxes
~~~~~~~~~

The sandbox area is a location where users can place their own repositories (i.e. directory trees having the same structure as the core and extras area) into.
Currently, the sandboxes are also stored in the SeqAn SVN repository but that will change in the near future.
Sandboxes can be generated using the ``skel.py`` :ref:`Code Generator <how-to-use-the-code-generator>`.

The following example shows how to create a user sandbox in the sandboxes area in an already existing Subversion repository.
We assume that ``https://svn.example.com/trunk`` is an empty directory in a Subversion repository.

.. code-block:: console

   seqan # cd sandbox
   sandbox # svn co https://svn.example.com/trunk sandbox_example
   ...
   sandbox # cd ..
   seqan # ./util/bin/skel.py --force repository sandbox/sandbox_example

Next, we can create an application from a simple template in this sandbox:

.. code-block:: console

   seqan # ./util/bin/skel.py app first_app sandbox/sandbox_example

Finally, commit this new sandbox into your Subversion repository:

.. code-block:: console

   seqan # cd sandbox/sandbox_example
   seqan # svn add *
   ...
   seqan # svn commit -m "Initial sandbox structure with one app."
   ...

Note that for the Subversion repository containing sandboxes, we recommend the following layout.
Using the classic SVN ``trunk``, ``tags``, ``branches`` structure allows for tagging releases or points of returns.
Furthermore, you can create folders parallel to those for documentation (for example a folder ``slides`` parallel to ``trunk``) without polluting your repository structure:

::

    Subversion repository root
      |-- trunk
      |     |-- CMakeLists.txt
      |     |-- apps
      |     |-- demos
      |     |-- include
      |     `-- tests
      |-- tags
      `-- branches

Documentation System
~~~~~~~~~~~~~~~~~~~~

The folders ``docs`` and ``docs2`` are used for building the documentation with the old and the new documentation system.
You can build them by going into the directory and then calling ``./make.sh``.
This will build the documentation into the sub directory ``html``:

.. code-block:: console

   seqan # cd docs
   docs # ./make.sh
   ...
   seqan # cd ../docs2
   docs2 # ./make.sh
   ...

If you want to include documentation for code from your sandbox then you can pass the path to the library (or library module) in your sandbox as a parameter to ``./make.sh``:

.. code-block:: console

   docs2 # ./make.sh ../sandbox/sandbox_example/include
   ...
   docs2 # ./make.sh ../sandbox/sandbox_example/include/jus_this_module
   ...

