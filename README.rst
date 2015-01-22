.. image:: https://readthedocs.org/projects/seqan/badge/?version=master
   :target: https://seqan.readthedocs.org/en/master?badge=master
   :alt: Documentation Status

SeqAn - The Library for Sequence Analysis
=========================================

What Is SeqAn?
--------------

SeqAn is an open source C++ library of efficient algorithms and data structures for the analysis of sequences with the focus on biological data.
Our library applies a unique generic design that guarantees high performance, generality, extensibility, and integration with other libraries.
SeqAn is easy to use and simplifies the development of new software tools with a minimal loss of performance.

License
-------

The SeqAn library itself, the tests and demos are licensed under the very permissing 3-clause BSD License.
The licenses for the applications themselves can be found in the LICENSE files.

Supported Platforms
-------------------

* Visual C++ 10 (2010), 11 (2012), 12 (2013)
* Clang >= 3.3 (e.g. Xcode >= 5)
* GCC >= 4.7 (e.g. Debian stable/wheezy)

Prerequisites
-------------

* Supported C++ compiler
* CMake (http://cmake.org/)

For the Impatient
-----------------

Assuming that you have checked out the repository already, all prerequisites are installed and you are on Linux or Mac OS X.

::

    # mkdir build/Debug
    # cd build/Debug
    # cmake ../.. -DCMAKE_BUILD_TYPE=Debug
    # make test_basic
    # ./tests/basic/test_basic
    ... the tests for module basic will run ...

Documentation Resources
-----------------------

* `Getting Started <http://seqan.readthedocs.org/en/master/Tutorial/GettingStarted.html>`_
* `Manual <http://seqan.readthedocs.org/en/master>`_
* `Tutorial <http://seqan.readthedocs.org/en/master/Tutorial.html>`_
* `How-Tos <http://seqan.readthedocs.org/en/master/HowTo.html>`_
* `API Documentation (stable) <http://docs.seqan.de/seqan/master/>`_


Contact
=======

* `Mailing List <https://lists.fu-berlin.de/listinfo/seqan-dev#subscribe>`_
* `GitHub Project (issues, source code) <https://github.com/seqan/seqan>`_
