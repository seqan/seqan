.. sidebar:: ToC

   .. contents::


.. _how-to-install-contribs-on-windows:

Installing Contribs On Windows
------------------------------

Download Contribs
~~~~~~~~~~~~~~~~~

The downloadable contribs contain precompiled library binaries (zlib, libbz2) for Windows by the supported compilers.
The contribs come in 32 bit and 64 bit variants.

* `Download contribs for 32 bit builds <http://ftp.seqan.de/contribs/seqan-contrib-D20130710-x86.zip>`_.
* `Download contribs for 64 bit builds <http://ftp.seqan.de/contribs/seqan-contrib-D20130710-x64.zip>`_.

You can install both variants in parallel if you want to do both 32 bit and 64 bit builds.

Extract Contribs
~~~~~~~~~~~~~~~~

Now, extract the downloaded ZIP file either to ``C:\Program Files`` or ``C:\``.

**After downloading the 64 bit variant**, you should now have a folder named ``C:\Program Files\seqan-contrib-D20130710-x64`` or a folder named ``C:\seqan-contrib-D20130710-x64``.

**After downloading the 32 bit variant**, you should now have a folder named ``C:\Program Files\seqan-contrib-D20130710-x86`` or a folder named ``C:\seqan-contrib-D20130710-x86``.

Re-run CMake
~~~~~~~~~~~~

You now have to re-run CMake to find the libraries.
You also have to remove the CMake Cache so it finds the new libraries.
You might also need to update your SeqAn Checkout.

The following assumes that your checkout is in ``c:\seqan-trunk`` and your build directory is ``c:\seqan-build\vs10``.

.. code-block:: console

    > cd c:\seqan-trunk
    > svn update .
    > cd c:\seqan-build\vs10
    > del CMakeCache.txt
    > cmake c:\seqan-trunk -G "Visual Studio 2010"
