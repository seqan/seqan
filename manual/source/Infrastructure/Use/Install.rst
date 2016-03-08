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

SeqAn is available natively on:

* **GNU/Linux**
    * Debian as "seqan-dev": ``apt install seqan-dev``
* **Mac OS X**
    * Macports as "seqan": ``port install seqan``
    * Homebrew as TODO
* **BSD**
    * FreeBSD as "seqan": ``pkg install seqan``

You should execute the above commands in a terminal as the ``root`` user or prefix them with ``sudo``.


Library Package
---------------

First you need to download the most recent "library package" from http://packages.seqan.de and extract its contents. Now copy the `include` and `share` folders to their target location. This could be one the following:

* ``/usr/local`` so they are available system-wide and automatically foundby your program [requires root or sudo]
* ``/opt/seqan`` available system-wide and easy to remove again [requires root or sudo]
* ``~/devel/seqan`` some place in your home directory [does not require root or sudo]

In any case it is important to remember where you installed it to.

Full Sources
------------

Make sure that you have git installed. For the operating systems mentioned above it can usually be achieved by using the respective command with `git` as package name.

For Windows there is Git client and shell available `here <https://windows.github.com/>`_.

Next create the required folders and clone our master branch:

.. code-block:: console

    ~ # mkdir -p ~/devel
    ~ # cd ~/devel
    ~ # git clone https://github.com/seqan/seqan.git seqan


You can update this branch at a later point by running ``git pull`` in ``~/devel/seqan`` .
