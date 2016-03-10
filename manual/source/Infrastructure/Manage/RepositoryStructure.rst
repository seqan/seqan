.. sidebar:: ToC

    .. contents::

.. _infra-manage-repo:

The SeqAn Repository
====================

This article describes the SeqAn repository structure and how to work with full SeqAn sources.

Getting Started
---------------

We assume that you have read :ref:`infra-use-install` and have cloned or unzipped the **full SeqAn sources** to ``~/devel/seqan`` (not the "library sources" described in other places).

SeqAn supports the usual CMake build types and we recommend :ref:`using multiple build directories <infra-use-cmake-build-dirs>`. Start like this:

.. code-block:: console

    # mkdir -p ~/devel/seqan-build/release
    # cd ~/devel/seqan-build/release
    # cmake ../../seqan -DCMAKE_BUILD_TYPE=Release


In addition to ``CMAKE_BUILD_TYPE`` there is also the ``SEQAN_BUILD_SYSTEM`` parameter which can be one of

#. ``DEVELOP`` -- all build targets (apps, demos, tests) and documentation (dox, manual) are created (the default).
#. ``SEQAN_RELEASE_LIBRARY`` -- only dox and library targets are created.
#. ``SEQAN_RELEASE_APPS`` -- all app targets are created, but nothing else.
#. ``APP:$APPNAME`` -- only a single app target is created for the chosen app.

All build systems other than ``DEVELOP`` are only relevant to :ref:`packaging releases <infra-manage-deploy>`.

As usual, calling ``make $TARGET`` will build a single target and just ``make`` will build all targets. On Windows, run ``cmake --build . --target $TARGET`` or just ``cmake --build`` instead.


Overview
--------

The main repository structure is shown in the following picture.

::

    seqan
      |-- CMakeLists.txt    CMake script file.
      |
      |-- LICENSE           Top-Level Information Files
      |-- README.rst
      |
      |-- apps              Applications
      |
      |-- demos             Demos
      |
      |-- dox               API documentation system
      |
      |-- include/seqan     SeqAn header ("the library")
      |
      |-- manual            Manuals
      |
      |-- tests             Unit tests for library modules
      |
      `-- util              Miscellaneous and Utility Code

The repository root contains some **information files** such as the ``LICENSE`` and ``README.rst``.
The other folders are as follows:

apps
----

The apps folder contains many applications, and
each application directory contains one ``CMakeLists.txt`` file and the files for compiling at least one binary.
Usually, apps have tests, too.
In this case, there is a subdirectory ``tests``.
Writing application tests is covered in detail in the article :ref:`infra-manage-write-app-tests`.

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


Note that some applications have binary names (make targets) that are not identical to the app-name, e.g. yara has ``yara_mapper`` and ``yara_indexer``.


demos
-----

The demos are short programs and code snippets that are used in the dox or the manual.
They serve as small examples and also functions as additional unit tests.


dox
---

The SeqAn API documentation is created using a customly-written system called *dox*.
It is very similar to doxygen, you can find out more about the syntax in :ref:`infra-contribute-dox`.

You can build the documentation in the `dox` subfolder of the *source folder*:

.. code-block:: console

   ~   # cd ~/devel/seqan/dox
   dox # ./dox_only.sh

This will build the documentation into the sub directory ``html``.


include/seqan
---------------

This is the actual library consisting of multiple modules:

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

.. note:: Header only library

   Remember that SeqAn is a template library that consists entirely of headers.
   No build steps are required for building the library and no shared objects will be created.

manual
------

The SeqAn manual is created using the `Sphinx <http://sphinx-doc.org/>`_ documentation system.

Follow these instructions to set up a local sphinx environment and build the manual:

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

tests
-----

The folder ``tests`` contains the unit tests for the library modules.

For each library module, there is a directory below ``tests`` with the same name that contains the tests for this module.
Simpler modules have one tests executable, whereas there might be multiple tests executables for larger modules.
For example, the module ``index`` has multiple test programs ``test_index_qgram``, ``test_index_shapes`` etc.
Writing tests is explained in detail in the article :ref:`infra-manage-write-unit-tests`.

