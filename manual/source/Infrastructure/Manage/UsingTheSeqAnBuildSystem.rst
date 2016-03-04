.. sidebar:: ToC

    .. contents::

.. _infra-manage-deploy:

Library and App releases
========================

There are three different "packaging targets":

#. a source package of the SeqAn library (``-DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_LIBRARY``)
#. a package containing all apps (``-DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_APPS``)
#. a package containing single SeqAn app (``-DSEQAN_BUILD_SYSTEM=APP:$appname``)

We assume that you have read :ref:`infra-use-install`, are using Linux/Mac/Unix and have cloned or unzipped the full SeqAn sources to ``~/devel/seqan`` (not the "library sources" described in other places).

The instructions for all packaging targets are the same (replace ``$pack_target`` with the above string):

.. code-block:: console

    ~ # mkdir -p ~/devel/seqan-build/deploy
    ~ # cd ~/devel/seqan-build/deploy
    deploy # cmake ../../seqan -DSEQAN_BUILD_SYSTEM=$pack_target -DSEQAN_OFFICIAL_PKGS=1
    deploy # make package

This will create a ZIP-file and on unix also a tarball (`.tar.xz`) of the package.
``SEQAN_OFFICIAL_PKGS`` enables some optimizations on our official packages and also makes them static, i.e. they will not depend on external libraries (OpenMP, Zlib...).

TODO: Remove all the cruft until Windows?

Packaging Library Releases
--------------------------

Packaging the library and documentation is quite simple. Note that we
have to build the documentation using ``make dox`` before calling
``make package`` because of a `bug in
CMake <http://public.kitware.com/Bug/view.php?id=8438>`_ that prevents
us from doing it automatically.

The version is automatically detected from the constants in the
``seqan/version.h`` header. There also is a marker variable that marks
whether the checked out repository version has a version number or
whether it is a pre-release of the next version.

.. code-block:: console

    ~ # mkdir -p ~/devel/seqan-build/release_library
    ~ # cd ~/devel/seqan-build/release_library
    release_library # cmake ../../seqan-src -DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_LIBRARY
    release_library # make package

On Linux, this will build three archives:

.. code-block:: console

    release_library # ls -l seqan-library-pre1.4.0-Linux.*
    -rw-rw-r-- 1 USER GROUP 2357465 Nov 20 13:57 seqan-library-pre1.4.0-Linux.tar.xz
    -rw-rw-r-- 1 USER GROUP 5953550 Nov 20 13:57 seqan-library-pre1.4.0-Linux.zip

Let us look at the contents of one (they all contain the same files):

TODO(h-2): replace by tar tf

.. code-block:: console

    release_library # dpkg --contents seqan-library-pre1.4.0-Linux.deb
    drwxrwxr-x root/root         0 2012-11-20 13:57 ./usr/
    drwxrwxr-x root/root         0 2012-11-20 13:57 ./usr/share/
    drwxrwxr-x root/root         0 2012-11-20 13:57 ./usr/share/seqan/
    drwxrwxr-x root/root         0 2012-11-20 13:57 ./usr/share/seqan/docs/
    drwxr-xr-x root/root         0 2012-11-20 13:57 ./usr/share/seqan/docs/html/
    -rw-r--r-- root/root      2012 2012-11-20 13:50 ./usr/share/seqan/docs/html/FUNCTION.prefix_Sum.html
    -rw-r--r-- root/root     24116 2012-11-20 13:50 ./usr/share/seqan/docs/html/SPEC_Super_Max_Repeats_Fast+_Iterator.html
    -rw-r--r-- root/root      1270 2012-11-20 13:50 ./usr/share/seqan/docs/html/MEMVAR_Triple_23i3.html
    ...
    -rw-r--r-- root/root      2940 2012-11-06 13:28 ./usr/share/doc/seqan/README
    -rw-r--r-- root/root      1517 2012-11-06 13:28 ./usr/share/doc/seqan/LICENSE
    drwxrwxr-x root/root         0 2012-11-20 13:57 ./usr/include/
    drwxrwxr-x root/root         0 2012-11-20 13:57 ./usr/include/seqan/
    drwxrwxr-x root/root         0 2012-11-20 13:57 ./usr/include/seqan/statistics/
    -rw-r--r-- root/root     24044 2012-11-06 13:28 ./usr/include/seqan/statistics/statistics_markov_model.h
    -rw-r--r-- root/root     15533 2012-11-06 13:28 ./usr/include/seqan/statistics/statistics_base.h
    ...

Packaging All Apps
------------------

It is simple to create a SeqAn Apps release:

.. code-block:: console

    ~ # git clone https://github.com/seqan/seqan seqan-src
    ~ # mkdir -p seqan-build/release_apps
    ~ # cd release_apps
    release_apps # cmake ../../seqan-src -DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_APPS
    release_apps # make package
    release_apps # ls -l seqan-apps-pre1.4.0-Linux*
    -rw-rw-r-- 1 USER GROUP 532 Nov 20 14:22 seqan-apps-pre1.4.0-Linux.deb
    -rw-rw-r-- 1 USER GROUP  42 Nov 20 14:22 seqan-apps-pre1.4.0-Linux.tar.bz2
    -rw-rw-r-- 1 USER GROUP  22 Nov 20 14:22 seqan-apps-pre1.4.0-Linux.zip

The contents of the archives is as follows:

.. code-block:: console

    release_library # dpkg --contents seqan-apps-pre1.4.0-Linux.deb
     dpkg --contents seqan-apps-pre1.4.0-Linux.deb
    drwxrwxr-x root/root         0 2012-11-20 14:30 ./usr/
    drwxrwxr-x root/root         0 2012-11-20 14:30 ./usr/bin/
    -rwxr-xr-x root/root   2253741 2012-11-20 14:27 ./usr/bin/masai_mapper
    -rwxr-xr-x root/root    191351 2012-11-20 14:24 ./usr/bin/tree_recon
    -rwxr-xr-x root/root    349878 2012-11-20 14:26 ./usr/bin/param_chooser
    ...
    drwxrwxr-x root/root         0 2012-11-20 14:30 ./usr/share/
    drwxrwxr-x root/root         0 2012-11-20 14:30 ./usr/share/doc/
    drwxrwxr-x root/root         0 2012-11-20 14:30 ./usr/share/doc/tree_recon/
    drwxrwxr-x root/root         0 2012-11-20 14:30 ./usr/share/doc/tree_recon/example/
    -rw-r--r-- root/root       475 2012-11-20 13:32 ./usr/share/doc/tree_recon/example/example.dist
    -rw-r--r-- root/root        20 2012-11-20 13:32 ./usr/share/doc/tree_recon/README
    -rw-r--r-- root/root       843 2012-11-20 13:32 ./usr/share/doc/tree_recon/LICENSE
    ...
    drwxrwxr-x root/root         0 2012-11-20 14:30 ./usr/share/doc/razers3/
    drwxrwxr-x root/root         0 2012-11-20 14:30 ./usr/share/doc/razers3/example/
    -rw-r--r-- root/root       105 2012-11-06 13:28 ./usr/share/doc/razers3/example/reads2.fa
    -rw-r--r-- root/root       985 2012-11-06 13:28 ./usr/share/doc/razers3/example/genome.fa
    -rw-r--r-- root/root       105 2012-11-06 13:28 ./usr/share/doc/razers3/example/reads.fa
    -rw-r--r-- root/root     23338 2012-11-06 13:28 ./usr/share/doc/razers3/README
    -rw-r--r-- root/root      1044 2012-11-20 13:32 ./usr/share/doc/razers3/LICENSE

Packaging Individual Apps
-------------------------

The release manager would check out an app in a specific revision, e.g.
through a tag or the current master version:

.. code-block:: console

    ~ # git clone -b yara-v0.9.2 https://github.com/seqan/seqan yara-v0.9.2
    ~ # mkdir yara-v0.9.2-build
    ~ # cd yara-0.v9.2-build
    yara-0.9.2-build # cmake ../yara-v0.9.2 -DSEQAN_BUILD_SYSTEM=APP:yara \
                          -DSEQAN_APP_VERSION=0.9.2
    yara-0.9.2-build # make package

On Unix, this will create a Tarball (``.tar.bz2``) and a ZIP file with
the binaries, documentation, and example files:

.. code-block:: console

    yara-0.9.2-build # ls -l yara-0.9.2-Linux-x86_64.*
    -rw-rw-r-- 1 USER GROUP  918587 Jan 16 18:15 yara-0.9.2-Linux-x86_64.tar.bz2
    -rw-rw-r-- 1 USER GROUP 1238990 Jan 16 18:15 yara-0.9.2-Linux-x86_64.zip

The packages have the following structure:

.. code-block:: console

    yara-0.9.2-build # tar tjf yara-0.9.2-Linux-x86_64.tar.bz2
    yara-0.9.2-Linux-x86_64/bin/yara_mapper
    yara-0.9.2-Linux-x86_64/bin/yara_indexer
    yara-0.9.2-Linux-x86_64/LICENSE
    yara-0.9.2-Linux-x86_64/README.rst




Windows Notes
-------------

TODO fuer cpockrandt: aktualisieren

The descriptions above apply to Linux/Mac systems.
On Windows we can use the GitHub client which can be downloaded `here <https://windows.github.com>`_.
Following the installation instructions will install a GitHub GUI client to manage your repository and a command line tool called ``Git Shell`` which emulates a unix like system so we can use the same commands as described before.

However, the main difference is that when building with the Visual Studio tools, one does not use ``make`` for building applications.
When developing, users can simply open the generated Visual Studio ``*.sln`` solution files and then use Visual Studio for building the applications.
When packaging, users can use the ``msbuild`` command as described below.

As an example, we adapt the description of creating an SeqAn application release on Windows.
The next steps are typed into the Command Prompt (``Start > All Programs > GitHub, Inc > Git Shell``).

.. code-block:: console

    ~ # git clone https://github.com/seqan/seqan seqan-src
    ~ # mkdir seqan-build
    ~ # cd seqan-build
    seqan-build # cmake ../seqan-src -DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_APPS

You can then open the generated ``seqan.sln`` file in ``C:\seqan-build`` with Visual Studio and build the packages from there.

Alternatively, ``msbuild`` can be used.
This program is only available when using the Visual Studio Command Prompt.
For Visual Studio 2010, you can start it through the start menu as follows:
``Start > Programs > Microsoft Visual Studio 2010 > Visual Studio Tools > Visual Studio Command Prompt 2010``.
For other Visual Studio versions, the path is similar.
If you want 64 bit builds then you have to start ``Visual Studio x86 Win64 Command Prompt (2010)``.

.. code-block:: console

    C:\> cd seqan-build
    C:\seqan-build> msbuild /p:Configuration=Release PACKAGE.vcxproj

This will create a ZIP file with the app build of the seqan apps.

Note that you could also input the first part of commands from this example into the Visual Studio Command Prompt.


Downstream Packaging
--------------------

These are some guidelines for creating SeqAn packages for operating system specific paackaging
systems, like *apt* (Debian/Ubuntu), , in GNU/Linux distr
