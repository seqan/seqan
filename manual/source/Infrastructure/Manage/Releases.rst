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


This will create a ZIP-file and on unix also a tarball (``.tar.xz``) of the package.
For the official packages distributed by us we also set ``-DSEQAN_STATIC_APPS=1 -DSEQAN_ARCH_SSE4=1`` which means that  binaries contain optimizations and are built statically, i.e. they will not depend on external libraries (OpenMP, Zlib...). For downstream packages you will likely not want to set these options.

.. note::

    Especially when creating packages, make sure that the cmake generator and/or compiler are the ones you want!

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

