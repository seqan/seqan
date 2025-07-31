.. sidebar:: ToC

   .. contents::


.. _infra-use-install:

Installing SeqAn
================

There are different ways to install SeqAn, we recommend to try these in the given order:

#. Native package management of the operating system.
#. Using the full sources from our github repository.

If possible, use the first option. If SeqAn is not available for your operating system, or if it is outdated, use the second option.

Native package management
-------------------------

SeqAn is available natively on the following platforms.

.. tip::

    Before you install, please make sure that the version supplied is not completely out of date (a difference of 0.1.* is okay, but if the difference is bigger use the `Library Package`_ below).
    The current version of SeqAn is always shown on the `GitHub repository <https://github.com/seqan/seqan/releases/latest>`__ and the version available on your platform is usually displayed in the info-link below.

.. |br| raw:: html

    <br/>

+-------------------------+----------------+-----------------------------------------+------------------------------------------------------------------------------------------------------------+
| Operating System        | Package Name   | Command                                 | links                                                                                                      |
+============+============+================+=========================================+============================================================================================================+
| **L** |br| | Arch       | seqan (AUR)    |  *depends*                              | `info <https://aur.archlinux.org/packages/seqan/>`__                                                       |
| **I** |br| +------------+----------------+-----------------------------------------+------------------------------------------------------------------------------------------------------------+
| **N** |br| | Debian     | libseqan2-dev  | ``apt install libseqan2-dev``           | `info <https://packages.debian.org/search?keywords=libseqan2-dev>`__                                       |
| **U** |br| +------------+----------------+-----------------------------------------+------------------------------------------------------------------------------------------------------------+
| **X** |br| | Ubuntu     | libseqan2-dev  | ``apt install libseqan2-dev``           | `info <https://packages.ubuntu.com/search?keywords=libseqan2-dev&searchon=names&suite=all&section=all>`__  |
+------------+------------+----------------+-----------------------------------------+------------------------------------------------------------------------------------------------------------+
| **M** |br| | Homebrew   | seqan          | ``brew install brewsci/bio/seqan@2``    | `info <https://github.com/brewsci/homebrew-bio/blob/develop/Formula/seqan%402.rb>`__                       |
| **A** |br| +------------+----------------+-----------------------------------------+------------------------------------------------------------------------------------------------------------+
| **C** |br| | MacPorts   | seqan          | ``port install seqan2``                 | `info <https://github.com/macports/macports-ports/blob/master/science/seqan2/Portfile>`__                  |
+------------+------------+----------------+-----------------------------------------+------------------------------------------------------------------------------------------------------------+
| **B** |br| | FreeBSD    | seqan          | ``pkg install seqan``                   | `info <https://freshports.org/biology/seqan>`__                                                            |
| **S** |br| +------------+----------------+-----------------------------------------+------------------------------------------------------------------------------------------------------------+
| **D** |br| |            |                |                                         |                                                                                                            |
+------------+------------+----------------+-----------------------------------------+------------------------------------------------------------------------------------------------------------+

You should execute the above commands in a terminal as the ``root`` user or prefix them with ``sudo``. If you have problems installing the package on your operating system, or it is outdated, please open an issue on https://github.com/seqan/seqan.


Library Package
---------------

First you need to download the most recent "library package" from https://github.com/seqan/seqan/releases and extract its contents. Now copy the `include` and `share` folders to their target location. This could be one of the following:

* ``/usr/local`` so they are available system-wide and automatically found by your program [requires root or sudo]
* ``/opt/seqan`` available system-wide and easy to remove again [requires root or sudo]
* ``~/devel/seqan`` some place in your home directory [does not require root or sudo]

In any case it is important to remember where you installed it to.

Full Sources
------------

Make sure that you have git installed. For the operating systems mentioned above it can usually be achieved by using the respective command with `git` as package name.

For Windows there is Git client and shell available `here <https://windows.github.com/>`__.

Next create the required folders and clone our main branch:

.. code-block:: console

    ~ # mkdir -p ~/devel
    ~ # cd ~/devel
    ~ # git clone https://github.com/seqan/seqan.git seqan


You can update this branch at a later point by running ``git pull`` in ``~/devel/seqan`` .
