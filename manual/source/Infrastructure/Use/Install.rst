.. sidebar:: ToC

   .. contents::


.. _infra-use-install:

Installing SeqAn
================

There are different ways to install SeqAn, we recommend to try these in the given order:

#. Native package management of the operating system.
#. Unpacking the library package from http://packages.seqan.de
#. Using the full sources from our github repository.

If possible, use the first option. If SeqAn is not available for your operating system, or if it is outdated, use the second option.

Use the third option if you want to use the master or develop branch which might contain bug-fixes and new features.

Native package management
-------------------------

SeqAn is available natively on the following platforms.

.. tip::

    Before you install, please make sure that the version supplied is not completely out of date (a difference in the third digit is ok, but if the difference is bigger use the `Library Package`_ below).
    The current version of SeqAn is always shown on the `SeqAn-Homepage <http://www.seqan.de/>`__ and the version available on your platform is usually displayed in the info-link below.

* **GNU/Linux**
    * Arch in the `Arch User Repository <https://aur.archlinux.org/packages/?O=0&K=seqan>`__
    * Debian as "seqan-dev": ``apt install seqan-dev`` [`info <https://packages.debian.org/stable/seqan-dev>`__ | `contact <mailto:debian-med-packaging()lists.alioth.debian.org>`__]
    * Fedora as "seqan-devel": ``yum install seqan-devel`` [`info <https://apps.fedoraproject.org/packages/seqan-devel>`__ | `contact <mailto:sagitter()fedoraproject.org>`__]
    * Ubuntu as "seqan-dev": ``apt install seqan-dev`` [`info <http://packages.ubuntu.com/xenial/seqan-dev>`__ | `contact <mailto:ubuntu-motu()lists.ubuntu.com>`__]
* **Mac OS X**
    * Homebrew as "seqan" in `science <http://brew.sh/homebrew-science/>`__ : ``brew install homebrew/science/seqan`` [`info <http://braumeister.org/repos/Homebrew/homebrew-science/formula/seqan>`__ | `contact <mailto:tim()tim-smith.us>`__]
    * Macports as "seqan": ``port install seqan`` [`info <https://trac.macports.org/browser/trunk/dports/science/seqan/Portfile>`__ | `contact <mailto:rene.rahn()fu-berlin.de>`__]

* **BSD**
    * FreeBSD as "seqan": ``pkg install seqan`` [`info <http://freshports.org/biology/seqan>`__ | `contact <mailto:h2+fbsdports()fsfe.org>`__]

You should execute the above commands in a terminal as the ``root`` user or prefix them with ``sudo``. If you have problems installing the package on your operating system, or it is outdated, please write to the contact shown above (and replace ``()`` in the e-mail-address with ``@``).


Library Package
---------------

First you need to download the most recent "library package" from http://packages.seqan.de and extract its contents. Now copy the `include` and `share` folders to their target location. This could be one of the following:

* ``/usr/local`` so they are available system-wide and automatically found by your program [requires root or sudo]
* ``/opt/seqan`` available system-wide and easy to remove again [requires root or sudo]
* ``~/devel/seqan`` some place in your home directory [does not require root or sudo]

In any case it is important to remember where you installed it to.

Full Sources
------------

Make sure that you have git installed. For the operating systems mentioned above it can usually be achieved by using the respective command with `git` as package name.

For Windows there is Git client and shell available `here <https://windows.github.com/>`__.

Next create the required folders and clone our master branch:

.. code-block:: console

    ~ # mkdir -p ~/devel
    ~ # cd ~/devel
    ~ # git clone https://github.com/seqan/seqan.git seqan


You can update this branch at a later point by running ``git pull`` in ``~/devel/seqan`` .
