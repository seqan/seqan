.. sidebar:: ToC

    .. contents::

.. _infra-use-install-dependencies:

Installing Dependencies
=======================

SeqAn can optionally make use of ZLIB and BZip2. This is relevant mostly for Input/Output.
Depending on your operating system you may need to install extra packages of these libraries or their headers.

GNU/Linux
---------

It depends on your distribution whether these packages are installed by default or not.

On Debian, Ubuntu, Mint and similar distributions:

.. code-block:: console

    # sudo apt install zlib1g-dev libbz2-dev

Mac and BSD
-----------

Nothing needs to be done, the libraries and their headers are pre-installed.

Windows
-------

The downloadable contribs contain precompiled library binaries (zlib, libbz2) for Windows by the supported compilers.
The contribs come in 32 bit and 64 bit variants.

* `Download contribs for 32 bit builds <http://ftp.seqan.de/contribs/seqan-contrib-D20160115-x86.zip>`_.
* `Download contribs for 64 bit builds <http://ftp.seqan.de/contribs/seqan-contrib-D20160115-x64.zip>`_.

You can install both variants in parallel if you want to do both 32 bit and 64 bit builds.
Previous contribs packages are available `here <http://ftp.seqan.de/contribs/>`__ and `here <http://svn.mi.fu-berlin.de/seqan-contrib/>`__ you can find the code for building contribs with new VS versions, for example.

Now, extract the downloaded ZIP file either to ``C:\Program Files`` or ``C:\``.

**After downloading the 64 bit variant**, you should now have a folder named ``C:\Program Files\seqan-contrib-D20160115-x64`` or a folder named ``C:\seqan-contrib-D20130710-x64``.

**After downloading the 32 bit variant**, you should now have a folder named ``C:\Program Files\seqan-contrib-D20160115-x86`` or a folder named ``C:\seqan-contrib-D20130710-x86``.

