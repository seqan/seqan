.. sidebar:: ToC

   .. contents::


.. _build-manual-using-the-seqan-build-system:

Using the SeqAn Build System
----------------------------

We describe the SeqAn build system from three perspectives:

* The **app user** who just wants to compile a couple of SeqAn applications from the SeqAn SVN repository.
* The **SeqAn release manager** who wants to create SeqAn releases.
* The **SeqAn developer** who wants to write his own applications using the SeqAn build system.

But first, we will give a short overview of the repository and how
versioning applications and the whole project works.

Repository Structure and Versioning
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The SeqAn library including all applications, demos and documentation are hosted under the repository https://github.com/seqan/seqan.
For a more detailed view over the repository structure please read this :ref:`document <infrastructure-repository-structure>`.

The SeqAn workflow is a mix between the **git-workflow** and the **forking-workflow** as described :ref:`here <infrastructure-seqan-git-workflow>`.
Please read the document and all linked sources carefully before you start developing tools and modules for SeqAn.

Note that there is no separation between apps and the library.
This means that all apps are released together with a new **library release** and are updated along with the chaneges in the ``develop`` branch.
A **library release** is achieved by tagging the corresponding commit to the ``master`` with the new version number, e.g. ``seqan-v2.0.0`` as described in the git-workflow.

Independently of this, an **app release** can be performed by tagging the new version with an increased version number, e.g. as ``yara-v0.9.0`` for the app Yara in version 0.9.0.
The tagged commit can either point to the ``master`` or to the ``develop`` branch, depending where it was applied.

Note that tags are final and a new tag has to be created if any code is to be changed.

User Perspective
~~~~~~~~~~~~~~~~

The user can clone either the ``master`` or the ``develop`` branch or any tagged version (e.g. ``yara-v0.9.0`` or ``seqan-v1.4.2``) to his local computer. 
The user could then proceed as the developer (see below) but there are dedicated modes in the SeqAn build system for easier installation. 
A user might also want to install the library to an include folder. 
We will look at both use cases.

User App Installation
^^^^^^^^^^^^^^^^^^^^^

Note that we assume Unixoid systems in this document and only refer to makefile based build systems.
The easiest way to install an application is described in the :ref:`getting started tutorials using linux makefiles <tutorial-getting-started-linux-makefiles>`.
By default the binaries are deployed in the bin folder of the build directory, e.g., ``${HOME}/Development/build-seqan/release/bin``.

However, it will be more convenient for the user to build the app and then install it, for example to a certain directory like ``~/local/bin/app``:
Here is an example for the application Razers 3.

.. code-block:: console

    ~ # git clone https://github.com/seqan/seqan seqan-src
    ~ # mkdir -p seqan-build/release-razers3
    ~ # cd seqan-build/release-razers3
    release-razers3 # cmake ../../seqan-src -DCMAKE_INSTALL_PREFIX=~/local/bin/razers3 \
                      -DSEQAN_BUILD_SYSTEM=APP:razers3
    release-razers3 # make install

.. hint::
    
    
	The user can of course install any tagged version by using the command
	
	.. code-block:: console
	    
	    # git clone -b <tag> https://github.com/seqan/seqan tag-src

After executing this, the user will find the following structure in ``~/local/bin/razers3``, including the example files and documentation.

.. code-block:: console

    razers3 # tree ~/local/bin/razers3
    /home/${USER}/local/bin/razers3/
    ├── bin
    │   └── razers3
    ├── example
    │   ├── genome.fa
    │   ├── reads2.fa
    │   └── reads.fa
    ├── LICENSE
    └── README

User Library Installation
^^^^^^^^^^^^^^^^^^^^^^^^^

The user could also want to install the library headers only.
The checkout step is the same as above, but he has to create a new build directory and execute CMake with different parameters. 
The library will be installed to ``~/local/seqan``.

.. code-block:: console

    ~ # git clone https://github.com/seqan/seqan seqan-src
    ~ # mkdir -p seqan-build/library_only
    ~ # cd seqan-build/library_only
    library_only # cmake ../../seqan-src -DCMAKE_INSTALL_PREFIX=~/local/seqan \
                     -DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_LIBRARY
    library_only # make dox
    library_only # make install

The user can now find the SeqAn library in ``~/local/seqan/include``:

.. code-block:: console

    library_only # tree ~/local/seqan/
    /home/${USER}/local/seqan/
    ├── include
    │   └── seqan
    │       ├── align
    │       │   ├── align_base.h
    │       │   ├── align_cols.h
    │       │   ├── align_config.h
    │       │   ├── align_iterator_base.h
    │       │   ├── alignment_algorithm_interface.h
    │       │   ├── alignment_algorithm_tags.h
    │       │   ├── alignment_operations.h
    │       │   ├── align_metafunctions.h
    │       │   ├── align_traceback.h
    │       │   ├── gap_anchor.h
    ...
    │       ├── system.h
    │       └── version.h
    └── share
        └── doc
            └── seqan
                ├── LICENSE
                └── README

SeqAn Release Manager Perspective
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The SeqAn release manager wants to create release packages of (1)
individual apps from the SeqAn repository, (2) create a SeqAn library
release that includes the library and documentation, and (3) create a
SeqAn apps release that contains the built apps. The manager wants to
build the binary packages for different platforms, e.g. 32 bit and 64
bit Linux and Windows, Mac Os X, etc.

We will give examples for Unixoid operating systems.

Note that the packaging described below can be automatized. App and
project releases can simply be tagged in the Subversion repository. A
script that runs nightly can then pick up new tags from the GitHub
repository and create binary packages for them. This can also automatize
nightly builds on different platforms without much work for the release
manager.

Packaging Individual Apps
^^^^^^^^^^^^^^^^^^^^^^^^^

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
    

Packaging Library Releases
^^^^^^^^^^^^^^^^^^^^^^^^^^

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

    ~ # git clone https://github.com/seqan/seqan seqan-src
    ~ # mkdir -p seqan-build/release_library
    ~ # cd seqan-build/release_library
    release_library # cmake ../../seqan-src -DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_LIBRARY
    release_library # make dox
    release_library # make package

On Linux, this will build three archives:

.. code-block:: console

    release_library # ls -l seqan-library-pre1.4.0-Linux.*
    -rw-rw-r-- 1 USER GROUP 3367876 Nov 20 13:57 seqan-library-pre1.4.0-Linux.deb
    -rw-rw-r-- 1 USER GROUP 2357465 Nov 20 13:57 seqan-library-pre1.4.0-Linux.tar.bz2
    -rw-rw-r-- 1 USER GROUP 5953550 Nov 20 13:57 seqan-library-pre1.4.0-Linux.zip

Let us look at the contents of one (they all contain the same files):

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
    drwxrwxr-x root/root         0 2012-11-20 13:57 ./usr/include/seqan/random/
    -rw-r--r-- root/root     15590 2012-11-06 13:28 ./usr/include/seqan/random/ext_MersenneTwister.h
    -rw-r--r-- root/root      4767 2012-11-06 13:28 ./usr/include/seqan/random/random_rng_functor.h
    -rw-r--r-- root/root      5810 2012-11-06 13:28 ./usr/include/seqan/random/random_uniform.h
    -rw-r--r-- root/root      4796 2012-11-06 13:28 ./usr/include/seqan/random/random_normal.h
    -rw-r--r-- root/root      3879 2012-11-06 13:28 ./usr/include/seqan/random/random_shuffle.h
    ...

Packaging All Apps
^^^^^^^^^^^^^^^^^^

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

Nightly Builds
^^^^^^^^^^^^^^

It is also possible to create nightly builds of the library, all apps,
or individual apps. Simply define the CMake variable
``SEQAN_NIGHTLY_RELEASE`` to ``TRUE`` on the command line. In the
following examples, we skip the checkout step and simply show the CMake
and build steps:

One App
^^^^^^^

.. code-block:: console

    masai-build # cmake ../yara-v0.9.2 -DSEQAN_BUILD_SYSTEM=APP:yara \
                          -DSEQAN_NIGHTLY_RELEASE=TRUE
    masai-build # make package
    masai-build # ls -l yara-20121120-Linux-x86_64.*
    -rw-rw-r-- 1 USER GROUP  918587 Nov 20 14:11 yara-20121120-Linux-x86_64.tar.bz2
    -rw-rw-r-- 1 USER GROUP 1238990 Nov 20 14:11 yara-20121120-Linux-x86_64.zip
    masai-build # tar tjf masai-20121120-Linux-x86_64.tar.bz2
    masai-20121120-Linux-x86_64/bin/masai_mapper
    masai-20121120-Linux-x86_64/bin/masai_indexer
    masai-20121120-Linux-x86_64/README
    masai-20121120-Linux-x86_64/LICENSE

All Apps
^^^^^^^^

.. code-block:: console

    release_apps # cmake ../../seqan-src -DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_APPS \
                     -DSEQAN_NIGHTLY_RELEASE=TRUE
    release_apps # make package
    release_apps #  ls -l seqan-apps-20121120-*
    -rw-rw-r-- 1 USER GROUP 10232442 Nov 20 14:37 seqan-apps-20121120-Linux.deb
    -rw-rw-r-- 1 USER GROUP  8847407 Nov 20 14:37 seqan-apps-20121120-Linux.tar.bz2
    -rw-rw-r-- 1 USER GROUP 10266596 Nov 20 14:37 seqan-apps-20121120-Linux.zip

Library Only
^^^^^^^^^^^^

.. code-block:: console

    release_library # cmake ../../seqan-src -DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_LIBRARY \
                        -DSEQAN_NIGHTLY_RELEASE=TRUE
    release_library # make dox
    release_library # make package
    release_library # ls -l seqan-library-20121120-*
    -rw-rw-r-- 1 USER GROUP 3368034 Nov 20 14:07 seqan-library-20121120-Linux.deb
    -rw-rw-r-- 1 USER GROUP 2356769 Nov 20 14:07 seqan-library-20121120-Linux.tar.bz2
    -rw-rw-r-- 1 USER GROUP 5955755 Nov 20 14:06 seqan-library-20121120-Linux.zip

SeqAn Developer Perspective
~~~~~~~~~~~~~~~~~~~~~~~~~~~

SeqAn developers want to develop their own applications using SeqAn.
When they want to use the SeqAn build system, they can follow these 
instructions to (1) fork the application template from github, 
(2) setup their apps, and (3) create releases of the applications.

Getting the Template
^^^^^^^^^^^^^^^^^^^^

Getting the application template can be achieved by forking the project ``https://github.com/seqan/APP_TEMPLATE.git``. 
This repository contains a template structure for the application containing all necessary files and a starting point from which to begin the development.
One of the files already present is the template repository is the ``CMakeLists.txt`` file.
Since you will have to adjust the file to your project, let us have a look at the file in detail.
You can look up details in the `CMake documentation <http://www.cmake.org/cmake/help/v2.8.8/cmake.html>`_ in case that some CMake functions are not clear to you.

The file starts out with a header describing where the file lives and what it is for.
This is useful when having many ``CMakeLists.txt`` files open and you want to quickly identifyin the file in the current window.

.. code-block:: cpp

   # ===========================================================================
   #                  SeqAn - The Library for Sequence Analysis
   # ===========================================================================
   # File: src/CMakeLists.txt
   #
   # CMakeLists.txt file for my_app.
   # ===========================================================================

   cmake_minimum_required (VERSION 2.8.2)
   project (src_my_app)
   message (STATUS "Configuring src/my_app")

Then comes the section that searches for the app's dependencies.
By default, the app only depends on the package SeqAn.
By setting the variable ``SEQAN_FIND_DEPENDENCIES``, we can configure which dependencies the call to ``find_package (SeqAn REQUIRED)`` will try to find.
See the :ref:`build-manual-using-the-find-seqan-cmake-module` for more details.

.. code-block:: cmake

    # ----------------------------------------------------------------------------
    # Dependencies
    # ----------------------------------------------------------------------------

    # Search SeqAn and select dependencies.
    set (SEQAN_FIND_DEPENDENCIES NONE)
    find_package (SeqAn REQUIRED)

The call to ``find_package (SeqAn REQUIRED)`` will then set the
following variables that we will then use below to add the correct
parameters to the compiler and linker.

* ``SEQAN_INCLUDE_DIRS``: Required include directories for the headers.
  Pass to ``include_directories()``
* ``SEQAN_DEFINITIONS``: Additional precompiler macros to pass to the
  compiler. Pass to ``add_definitions()``
* ``SEQAN_CXX_FLAGS``: Additional C++ compiler flags. Extend
  ``CMAKE_CXX_FLAGS`` by this list.
* ``SEQAN_LIBRARIES``: The libraries to link against. Pass to
  ``target_link_libraries()`` for each target.

We then need one ``add_executable()`` call for each program executable
that we want to build. We also need to link the libraries into the
program.

.. code-block:: cmake

    # ----------------------------------------------------------------------------
    # Build Setup
    # ----------------------------------------------------------------------------

    # Add CXX flags found by find_package(SeqAn).
    set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS})

    # Add include directories.
    include_directories (${SEQAN_INCLUDE_DIRS})

    # Add definitions set by find_package(SeqAn).
    add_definitions (${SEQAN_DEFINITIONS})

    # Update the list of file names below if you add source files to your application.
    add_executable (dfi dfi.cpp)

    # Add dependencies found by find_package(SeqAn).
    target_link_libraries (dfi ${SEQAN_LIBRARIES})

We then configure the app for installation. Note that this is a distinct
step than configuring CPack for packaging. The following controls which
files to copy when calling ``make install``. CPack will use the result
of ``make install`` for creating its packages.

We first call ``seqan_setup_install_vars()`` (to set the variable
``SEQAN_PREFIX_SHARE_DOC``. This is required for installing
documentation and example files to ``share/${PROGRAM_NAME}`` when
building multiple apps and directly to the current directory ``.`` when
building only one app.

The macro ``seqan_setup_install_vars`` is specific to the SeqAn build
system.

The ``util/skel.py`` script will create files ``LICENSE`` and ``README``
for you. If you want to include additional files then you should use one
of the given ``install()`` calls. Install documentation to
``${SEQAN_PREFIX_SHARE_DOC}`` and examples to
``${SEQAN_PREFIX_SHARE_DOC}/example``.

.. code-block:: console

    # ----------------------------------------------------------------------------
    # Installation
    # ----------------------------------------------------------------------------

    # Set variables for installing, depending on the selected build type.
    if (NOT SEQAN_PREFIX_SHARE_DOC)
      seqan_setup_install_vars (dfi)
    endif (NOT SEQAN_PREFIX_SHARE_DOC)

    # Install dfi in ${PREFIX}/bin directory
    install (TARGETS dfi
             DESTINATION bin)

    # Install non-binary files for the package to "." for app builds and
    # ${PREFIX}/share/doc/dfi for SeqAn release builds.
    install (FILES LICENSE
                   README
             DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
    #install (FILES example/example.txt
    #         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

Then, we can use the macro ``seqan_add_app_test()`` from the SeqAn build system to register app tests.
If you want to add an app test for your program then simply uncomment the ``seqan_add_app_test()`` call and follow the instructions in :ref:`how-to-write-app-tests` to write such an app tests.

.. code-block:: console

    # ----------------------------------------------------------------------------
    # App Test
    # ----------------------------------------------------------------------------

    #seqan_add_app_test(dfi)

Finally, we configure the application packaging system for building individual apps.

.. code-block:: console

    # ----------------------------------------------------------------------------
    # CPack Install
    # ----------------------------------------------------------------------------

    if (SEQAN_BUILD_SYSTEM STREQUAL "APP:my_app")
      set (CPACK_PACKAGE_NAME "my_app")
      set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "My App - Catch Summary")
      set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.name@example.com>")
      set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

      seqan_configure_cpack_app(my_app "My App")
    endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:my_app")

.. hint::

    If you use the markdown feature for your ``README`` with the file ending ``*.rst``, then you need to explicitly tell **CPack**, which the correct README file is.
    You can do this by adding the following line to the **CPack Install** section.
    
    .. code-block:: console
        
        set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.rst")
        
    Also make sure to replace all occurrences of ``README`` with ``README.rst`` in the **INSTALLATION** section.

Building Apps
^^^^^^^^^^^^^

Simply use CMake to generate project files for the whole SeqAn
repository. Let us say that we want to build the app
``my_app``:

.. code-block:: console

    ~ # mkdir -p seqan-build/release
    ~ # cd seqan-build/release
    release # cmake ../../seqan-src
    release # make my_app

Note that the default build type is the release mode. 
The binaries will be built with optimization and without debug symbols. 
To build apps with debug symbols and without optimization with Makefiles, use the CMake paraemter ``-DCMAKE_BUILD_TYPE=Debug``.
When using IDE files such as for Xcode, you can select the optimization state from within the IDE.

.. code-block:: console

    Release # cd ../..
    ~ # mkdir -p seqan-build/debug
    ~ # cd seqan-build/debug
    debug # cmake ../../seqan-src
    debug # make my_app

Windows Notes
~~~~~~~~~~~~~

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