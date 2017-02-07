.. image:: https://readthedocs.io/projects/seqan/badge/?version=develop
   :target: https://seqan.readthedocs.io/en/develop?badge=develop
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

Prerequisites
-------------------

Linux, Mac OSX, FreeBSD:
  * GCC ≥ 4.9 [GCC ≥ 5 recommended]
  * Clang/LLVM ≥ 3.5
  * Intel Compiler ≥ 16.0.2
Windows:
  * Visual C++ ≥ 14.0 / Visual Studio ≥ 2015
  * Intel Compiler ≥ 16.0.3 / Visual Studio ≥ 2015u2
  * Clang/C2 ≥ 3.8.0 / Visual Studio ≥ 2015u3 [experimental, requires CMake ≥ 3.6]

Architecture support:
  * 32bit and 64bit Intel/AMD officially supported
  * 64bit is highly recommended and some parts of SeqAn are optimized for ``POPCNT``, ``SSE4`` and/or ``AVX2`` instruction sets available on newer CPUs
  * many other architectures, including ARM, PowerPC and Sparc64 work, but receive less testing and optimizations

To build tests, demos, and official SeqAn applications you also need:
  * CMake ≥ 3.0 [CMake ≥ 3.4 recommended]

Some of the official applications might have additional requirements or work only on a subset of platforms.

Documentation Resources
-----------------------

* `Getting Started <http://seqan.readthedocs.io/en/master/Tutorial/GettingStarted>`_
* `Manual <http://seqan.readthedocs.io/en/master>`_
* `Tutorial <http://seqan.readthedocs.io/en/master/index.html#tutorials>`_
* `How-Tos <http://seqan.readthedocs.io/en/master/Tutorial/HowTo>`_
* `API Documentation (stable) <http://docs.seqan.de/seqan/master/>`_

Contact
=======

* `Mailing List <https://lists.fu-berlin.de/listinfo/seqan-dev#subscribe>`_
* `GitHub Project (issues, source code) <https://github.com/seqan/seqan>`_
