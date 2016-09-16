.. sidebar:: ToC

    .. contents::

.. _infra-manage-deploy:

Library and App releases
========================

There are three different "packaging targets":

#. a source package of the SeqAn library (``-DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_LIBRARY``)
#. a package containing all apps (``-DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_APPS``)
#. a package containing a single SeqAn app (``-DSEQAN_BUILD_SYSTEM=APP:$appname``)

We assume that you have read :ref:`infra-use-install` and have cloned or unzipped the **full SeqAn sources** to ``~/devel/seqan`` (not the "library sources" described in other places).

The instructions for all packaging targets are the same (replace ``$pack_target`` with the above string):

.. code-block:: console

    ~ # mkdir -p ~/devel/seqan-build/deploy
    ~ # cd ~/devel/seqan-build/deploy
    deploy # cmake ../../seqan -DSEQAN_BUILD_SYSTEM=$pack_target -DSEQAN_STATIC_APPS=1 -DSEQAN_ARCH_SSE4=1
    deploy # make package

On Windows, replace the last command with

.. code-block:: console

    deploy # cmake --build . --target PACKAGE

Depending on the platform this might create a ZIP-file, a tarball and/or a platform specific installer.

Official Packages
-----------------

We provide (1) a source package of SeqAn library; and for each official application (3) single binary packages for different operating systems and architectures.

.. note::

    Especially when creating packages, make sure that the cmake generator and/or compiler are the ones you want!

GNU/Linux, macOS & BSD
^^^^^^^^^^^^^^^^^^^^^^

* The binary packages should be built on the **oldest supported kernel** and with the **oldest supported GCC** compiler.
* The CMake version on the building system should be at least 3.1.
* Builds should be static (``-DSEQAN_STATIC_APPS=1``).
* There should be a 32Bit package, built on a 32Bit system or cross-compiled (``-DCMAKE_CXX_FLAGS="-m32"``).
* There should be a 64Bit package.
* There should be an optimized 64Bit build (``-DSEQAN_ARCH_SSE4=1``).
* For applications where it makes sense, a further optimized build *can* be provided (``-DSEQAN_ARCH_AVX2=1``)

Windows
^^^^^^^

* The binary packages should be built with the latest **Intel C++ Compiler** for performance and compatibility reasons (see :ref:`here <infra-use-cmake-build-dirs>`).
* There should be a 32Bit package, built on a 32Bit system or cross-compiled (see :ref:`here <infra-use-cmake-build-dirs>`).
* There should be a 64Bit package.

Downstream Packaging
--------------------

These are some guidelines for creating SeqAn packages for operating system specific packaging
systems, like *apt* (Debian/Ubuntu) or *rpm* (Fedora/RedHat/CentOS/SUSE).

Library Package
^^^^^^^^^^^^^^^

We recommend that downstream package maintainers provide one package named **seqan** that contains only the header-library and the api-docs and that is built from our *library packages* available here: http://packages.seqan.de

They have the advantage of not requiring any build steps, simply copy the ``include`` and ``share`` directories to the desired locations.

Application Package(s)
^^^^^^^^^^^^^^^^^^^^^^

Beyond that package maintainers have the choice to create either a single package called **seqan-apps** that contains all the applications *or* a seperate package per application (with the respective name of that app). Based on the above instructions this should be fairly easy to accomplish.

